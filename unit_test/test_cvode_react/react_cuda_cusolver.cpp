#include <cvode/cvode.h>                  /* prototypes for CVODE fcts., consts.      */
#include <nvector/nvector_cuda.h>         /* access to CUDA N_Vector                  */
#include <cusolver/cvode_cusolver_spqr.h>    /* access to cuSolver interface             */
#include <sundials/sundials_types.h>      /* definition of realtype                   */
#include <sundials/sundials_math.h>       /* contains the macros ABS, SUNSQR, and EXP */
#include <AMReX_MultiFab.H>
#include "test_react.H"
#include "test_react_F.H"
#include <iostream>

using namespace amrex;

void do_react(const int* lo, const int* hi,
	      amrex::Real* state, const int* s_lo, const int* s_hi,
	      const int ncomp, const amrex::Real dt)
{
  const int size_x = hi[0]-lo[0]+1;
  const int size_y = hi[1]-lo[1]+1;
  const int size_z = hi[2]-lo[2]+1;
  const int size_state = size_x * size_y * size_z;

  int idx_spec, idx_spec_old, idx_dens, idx_temp, idx_omegadot, idx_dens_hnuc;

  get_species_index(&idx_spec);
  get_species_old_index(&idx_spec_old);
  get_density_index(&idx_dens);
  get_temperature_index(&idx_temp);
  get_omegadot_index(&idx_omegadot);
  get_density_hnuc_index(&idx_dens_hnuc);

  int neqs, nspec_not_evolved, size_flat;

  get_number_equations(&neqs, &nspec_not_evolved);
  size_flat = neqs * size_state;
  int nspec_evolve;
  sk_get_nspec_evolve(&nspec_evolve);
  int size_rpar_per_cell;
  sk_get_num_rpar_comps(&size_rpar_per_cell);
  const int size_rpar = size_rpar_per_cell * size_state;

  UserData user_data;
  cudaMallocManaged(&user_data, sizeof(struct CVodeUserData));
  cudaMallocManaged(&user_data->rpar,
		    size_rpar * sizeof(amrex::Real));
  cudaMallocManaged(&user_data->dense_jacobians,
		    size_state * neqs * neqs * sizeof(amrex::Real));
  user_data->num_cells = size_state;
  user_data->num_eqs_per_cell = neqs;
  user_data->num_rpar_per_cell = size_rpar_per_cell;

  initialize_rpar_indices(user_data, nspec_not_evolved, neqs);
  zero_rpar_data(user_data, size_rpar);

  realtype reltol=1.0e-6, time=0.0e0, tout;

  realtype abstol_values[size_flat];
  amrex::Real* state_y;
  amrex::Real* state_final_y;
  cudaMallocManaged(&state_y, size_flat * sizeof(amrex::Real));
  
#ifdef PRINT_DEBUG
  cudaMallocManaged(&state_final_y, size_flat * sizeof(amrex::Real));
#else
  state_final_y = state_y;
#endif

  N_Vector y = NULL, yout=NULL;
  N_Vector abstol = NULL;
  
  cuSolver_method LinearSolverMethod = QR;
  int jac_number_nonzero;

  sk_get_csr_jac_nnz(&jac_number_nonzero);
  
  int csr_row_count[neqs+1];
  int csr_col_index[jac_number_nonzero];

  sk_get_csr_jac_rowcols(&csr_row_count[0], &csr_col_index[0]);
  
  void* cvode_mem = NULL;
  int flag;

  // Create NVectors
  y = N_VNew_Cuda(size_flat);
  yout = N_VNew_Cuda(size_flat);
  abstol = N_VNew_Cuda(size_flat);

  // Initialize y, abstol from flattened state
  int nzone = 0;
  for (int i=lo[0]; i<=hi[0]; i++) {
    for (int j=lo[1]; j<=hi[1]; j++) {
      for (int k=lo[2]; k<=hi[2]; k++) {
	// Put mass fractions into integration vector
	int scomp = 0;
        for (int n=idx_spec_old; n<idx_spec_old+nspec_evolve; n++) {
          get_state(state, s_lo, s_hi, ncomp, i, j, k, n, &state_y[nzone*neqs + scomp]);
	  abstol_values[nzone*neqs + scomp] = 1.0e-6;
	  scomp++;
        }

	// Temperature absolute tolerance
	abstol_values[nzone*neqs + nspec_evolve] = 1.0e-6;
	
	// Put temperature into integration vector
	get_state(state, s_lo, s_hi, ncomp, i, j, k, idx_temp,
		  &state_y[nzone*neqs + nspec_evolve]);

	// Energy absolute tolerance
	abstol_values[nzone*neqs + nspec_evolve + 1] = 1.0e-6;

	// Initialize energy to 0, we'll get it from the EOS
	state_y[nzone*neqs + nspec_evolve + 1] = 0.0e0;

	// Put density in user data
	get_state(state, s_lo, s_hi, ncomp, i, j, k, idx_dens,
		  &user_data->rpar[nzone*user_data->num_rpar_per_cell + user_data->irp_dens]);

	// Put unevolved mass fractions in user data
	scomp=0;
	for (int n=idx_spec_old+nspec_evolve; n<idx_spec_old+nspec_evolve+nspec_not_evolved; n++) {
	  get_state(state, s_lo, s_hi, ncomp, i, j, k, n,
		    &user_data->rpar[nzone*user_data->num_rpar_per_cell + user_data->irp_xn_not_evolved + scomp]);
	  scomp++;
	}

	// Set zone size dx in rpar(irp_dx)
	user_data->rpar[nzone*user_data->num_rpar_per_cell + user_data->irp_dx] = 1.0;

	// Set time offset to the simulation time of 0
	user_data->rpar[nzone*user_data->num_rpar_per_cell + user_data->irp_t0] = 0.0;

	nzone++;
      }
    }
  }

  // Prepare cell data for integration (X->Y, normalization, call EOS)
  initialize_system(&state_y[0], user_data);
  
  // Set CVODE initial values and absolute tolerances
  set_nvector_cuda(y, &state_y[0], size_flat);
  set_nvector_cuda(abstol, &abstol_values[0], size_flat);

  // Initialize CVODE
  cvode_mem = CVodeCreate(CV_BDF);
  flag = CVodeSetUserData(cvode_mem, static_cast<void*>(user_data));
  if (flag != CV_SUCCESS) amrex::Abort("Failed to set user data");
  flag = CVodeInit(cvode_mem, fun_rhs, time, y);
  if (flag != CV_SUCCESS) amrex::Abort("Failed to initialize CVode");  
  flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (flag != CV_SUCCESS) amrex::Abort("Failed to set tolerances");  
  flag = CVodeSetMaxNumSteps(cvode_mem, 150000);
  if (flag != CV_SUCCESS) amrex::Abort("Failed to set max steps");  

  // Initialize cuSolver Linear Solver
  flag = cv_cuSolver_SetLinearSolver(cvode_mem, LinearSolverMethod);
  flag = cv_cuSolver_CSR_SetSizes(cvode_mem, neqs, jac_number_nonzero, size_state);
  flag = cv_cuSolver_SetJacFun(cvode_mem, &fun_csr_jac);
  flag = cv_cuSolver_SystemInitialize(cvode_mem, &csr_row_count[0], &csr_col_index[0]);

  // Do Integration
  time = time + static_cast<realtype>(dt);
  flag = CVode(cvode_mem, time, yout, &tout, CV_NORMAL);

  cudaError_t cuda_status = cudaSuccess;
  cuda_status = cudaGetLastError();
#ifdef PRINT_DEBUG  
  std::cout << "After CVode, got CUDA last error of: " << cudaGetErrorString(cuda_status) << std::endl;
#endif  
  assert(cuda_status == cudaSuccess);


  if (flag != CV_SUCCESS) amrex::Abort("Failed integration");

#ifdef PRINT_DEBUG
  std::cout << "Desired end time = " << time << std::endl;
  std::cout << "Integrated to tout = " << tout << std::endl;
  std::cout << "Time difference = " << time-tout << std::endl;
#endif

  // Get Final State
  get_nvector_cuda(yout, &state_final_y[0], size_flat);

  // Finalize cell data to save
  finalize_system(&state_final_y[0], user_data);

  // Save Final State
  nzone = 0;
  for (int i=lo[0]; i<=hi[0]; i++) {
    for (int j=lo[1]; j<=hi[1]; j++) {
      for (int k=lo[2]; k<=hi[2]; k++) {
	// Put evolved mass fractions into state
	int scomp = 0;
        for (int n=idx_spec; n<idx_spec+nspec_evolve; n++) {
          set_state(state, s_lo, s_hi, ncomp, i, j, k, n, state_final_y[nzone*neqs + scomp]);
	  scomp++;
        }

	// Put unevolved mass fractions into state
	scomp=0;
	for (int n=idx_spec+nspec_evolve; n<idx_spec+nspec_evolve+nspec_not_evolved; n++) {
	  set_state(state, s_lo, s_hi, ncomp, i, j, k, n,
		    user_data->rpar[nzone*user_data->num_rpar_per_cell + user_data->irp_xn_not_evolved + scomp]);
	  scomp++;
	}

	// Put omegadot into state
	amrex::Real xn_start, xn_final, wscratch;
	scomp = 0;
	int n_spec_old;
        for (int n=idx_omegadot; n<idx_omegadot+nspec_evolve+nspec_not_evolved; n++) {
	  n_spec_old = n-idx_omegadot+idx_spec_old;
          get_state(state, s_lo, s_hi, ncomp, i, j, k, n_spec_old, &xn_start);
	  if (scomp < nspec_evolve) {
	    xn_final = state_final_y[nzone*neqs + scomp];
	  } else {
	    xn_final = user_data->rpar[nzone*user_data->num_rpar_per_cell + user_data->irp_xn_not_evolved + scomp - nspec_evolve];
	  }
	  wscratch = (xn_final - xn_start)/dt;
	  set_state(state, s_lo, s_hi, ncomp, i, j, k, n, wscratch);
	  scomp++;
        }

	// Set rho*Hnuc
	get_state(state, s_lo, s_hi, ncomp, i, j, k, idx_dens, &wscratch);
	wscratch = wscratch * state_final_y[nzone*neqs + nspec_evolve + 1]/dt;
	set_state(state, s_lo, s_hi, ncomp, i, j, k, idx_dens_hnuc, wscratch);

	// Print initial and final vectors
#ifdef PRINT_DEBUG
	std::cout << "initial integration vector = ";
	for (int ii = 0; ii < neqs; ii++)
	  std::cout << state_y[nzone*neqs + ii] << " ";
	std::cout << std::endl;	
	std::cout << "final integration vector = ";
	for (int ii = 0; ii < neqs; ii++)
	  std::cout << state_final_y[nzone*neqs + ii] << " ";
	std::cout << std::endl;
#endif

	nzone++;
      }
    }
  }

  // Free Memory
  cudaFree(user_data->dense_jacobians);
  cudaFree(user_data->rpar);
  cudaFree(user_data);
  cudaFree(state_y);
#ifdef PRINT_DEBUG
  cudaFree(state_final_y);
#endif
  N_VDestroy(y);
  N_VDestroy(yout);
  N_VDestroy(abstol);
  CVodeFree(&cvode_mem);

}


void initialize_rpar_indices(UserData user_data, const int nspec_not_evolved,
			     const int num_eqs_per_cell)
{
  int i = 0;
  user_data->irp_dens = i; i++;
  user_data->irp_cv = i; i++;
  user_data->irp_cp = i; i++;
  user_data->irp_xn_not_evolved = i; i+=nspec_not_evolved;
  user_data->irp_abar = i; i++;
  user_data->irp_zbar = i; i++;
  user_data->irp_eta = i; i++;
  user_data->irp_ye = i; i++;
  user_data->irp_cs = i; i++;
  user_data->irp_dx = i; i++;
  user_data->irp_t_sound = i; i++;
  user_data->irp_y_init = i; i+=num_eqs_per_cell;
  user_data->irp_self_heat = i; i++;
  user_data->irp_Told = i; i++;
  user_data->irp_dcvdt = i; i++;
  user_data->irp_dcpdt = i; i++;
  user_data->irp_t0 = i; i++;
  user_data->irp_energy_offset = i;
}


void zero_rpar_data(UserData user_data, const int size_rpar)
{
  for (int i=0; i<size_rpar; i++) {
    user_data->rpar[i] = 0.0;
  }
}


void initialize_system(realtype* y, UserData udata)
{
  cudaError_t cuda_status = cudaSuccess;
  cuda_status = cudaGetLastError();
#ifdef PRINT_DEBUG
  std::cout << "In initialize_system, got CUDA last error of: " << cudaGetErrorString(cuda_status) << std::endl;
#endif
  assert(cuda_status == cudaSuccess);  

  int numThreads = std::min(32, udata->num_cells);
  int numBlocks = static_cast<int>(ceil(((double) udata->num_cells)/((double) numThreads)));
  initialize_cell<<<numBlocks, numThreads>>>(y, udata);

  cuda_status = cudaDeviceSynchronize();
#ifdef PRINT_DEBUG
  std::cout << "In initialize_system, got cudaDeviceSynchronize error of: " << cudaGetErrorString(cuda_status) << std::endl;
#endif
  assert(cuda_status == cudaSuccess);
}


__global__ static void initialize_cell(realtype* y, UserData udata)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < udata->num_cells) {
    int offset = tid * udata->num_eqs_per_cell;
    int rpar_offset = tid * udata->num_rpar_per_cell;
    sk_initialize_cell_device(&y[offset], &udata->rpar[rpar_offset]);
  }
}


void finalize_system(realtype* y, UserData udata)
{
  cudaError_t cuda_status = cudaSuccess;
  cuda_status = cudaGetLastError();
#ifdef PRINT_DEBUG
  std::cout << "In finalize_system, got CUDA last error of: " << cudaGetErrorString(cuda_status) << std::endl;
#endif
  assert(cuda_status == cudaSuccess);  
  
  int numThreads = std::min(32, udata->num_cells);
  int numBlocks = static_cast<int>(ceil(((double) udata->num_cells)/((double) numThreads)));
  finalize_cell<<<numBlocks, numThreads>>>(y, udata);

  cuda_status = cudaDeviceSynchronize();
#ifdef PRINT_DEBUG
  std::cout << "In finalize_system, got cudaDeviceSynchronize error of: " << cudaGetErrorString(cuda_status) << std::endl;
#endif
  assert(cuda_status == cudaSuccess);
}


__global__ static void finalize_cell(realtype* y, UserData udata)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < udata->num_cells) {
    int offset = tid * udata->num_eqs_per_cell;
    int rpar_offset = tid * udata->num_rpar_per_cell;

    sk_finalize_cell_device(&y[offset], &udata->rpar[rpar_offset]);
  }
}


static void set_nvector_cuda(N_Vector vec, realtype* data, sunindextype size)
{
  realtype* vec_host_ptr = N_VGetHostArrayPointer_Cuda(vec);
  for (sunindextype i = 0; i < size; i++) {
    vec_host_ptr[i] = data[i];
  }
  N_VCopyToDevice_Cuda(vec);
}


static void get_nvector_cuda(N_Vector vec, realtype* data, sunindextype size)
{
  N_VCopyFromDevice_Cuda(vec);  
  realtype* vec_host_ptr = N_VGetHostArrayPointer_Cuda(vec);
  for (sunindextype i = 0; i < size; i++) {
    data[i] = vec_host_ptr[i];
  }
}


static int fun_rhs(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  cudaError_t cuda_status = cudaSuccess;
  cuda_status = cudaGetLastError();
#ifdef PRINT_DEBUG
  std::cout << "In fun_rhs, got CUDA last error of: " << cudaGetErrorString(cuda_status) << std::endl;
#endif
  assert(cuda_status == cudaSuccess);
  realtype* ydot_d = N_VGetDeviceArrayPointer_Cuda(ydot);
  realtype* y_d = N_VGetDeviceArrayPointer_Cuda(y);
  UserData udata = static_cast<CVodeUserData*>(user_data);
  int numThreads = std::min(32, udata->num_cells);
  int numBlocks = static_cast<int>(ceil(((double) udata->num_cells)/((double) numThreads)));
  fun_rhs_kernel<<<numBlocks, numThreads>>>(t, y_d, ydot_d,
					    user_data);

  cuda_status = cudaDeviceSynchronize();
#ifdef PRINT_DEBUG
  std::cout << "In fun_rhs, got cudaDeviceSynchronize error of: " << cudaGetErrorString(cuda_status) << std::endl;
#endif
  assert(cuda_status == cudaSuccess);

#ifdef PRINT_DEBUG
  // debugging rhs
  realtype ydot_host[udata->num_eqs_per_cell*2];

  // copy first and second system rhs to host
  cuda_status = cudaMemcpy(&ydot_host[0], ydot_d, sizeof(realtype) * udata->num_eqs_per_cell * 2, cudaMemcpyDeviceToHost);
  assert(cuda_status == cudaSuccess);
  
  // print out the first system rhs
  std::cout << "first system rhs = " << std::endl;
  for (int i = 0; i < udata->num_eqs_per_cell; i++)
    std::cout << ydot_host[i] << " ";
  std::cout << std::endl;

  // print out the second system rhs
  std::cout << "second system rhs = " << std::endl;
  for (int i = udata->num_eqs_per_cell; i < 2*udata->num_eqs_per_cell; i++)
    std::cout << ydot_host[i] << " ";
  std::cout << std::endl;
#endif
  
  return 0;
}


__global__ static void fun_rhs_kernel(realtype t, realtype* y, realtype* ydot,
				      void *user_data)
{
  UserData udata = static_cast<CVodeUserData*>(user_data);
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < udata->num_cells) {
    int offset = tid * udata->num_eqs_per_cell;
    int rpar_offset = tid * udata->num_rpar_per_cell;

    sk_f_rhs_device(&t, &y[offset], &ydot[offset], &udata->rpar[rpar_offset]);
  }
}


int fun_csr_jac(realtype t, N_Vector y, N_Vector fy,
		CV_cuSolver_csr_sys csr_sys, void* user_data)
{
  cudaError_t cuda_status = cudaSuccess;
  cuda_status = cudaGetLastError();

#ifdef PRINT_DEBUG
  std::cout << "In fun_csr_jac, got CUDA Last Error of: ";
  std::cout << cudaGetErrorString(cuda_status) << std::endl;
#endif

  assert(cuda_status == cudaSuccess);
  realtype* y_d   = N_VGetDeviceArrayPointer_Cuda(y);
  realtype* fy_d  = N_VGetDeviceArrayPointer_Cuda(fy);
  UserData udata = static_cast<CVodeUserData*>(user_data);

#ifdef PRINT_DEBUG
  std::cout << "num_cells = " << udata->num_cells << std::endl;
  std::cout << "size_per_subsystem = " << csr_sys->size_per_subsystem << std::endl;
  std::cout << "csr_number_nonzero = " << csr_sys->csr_number_nonzero << std::endl;
  std::cout << "number_subsystems = " << csr_sys->number_subsystems << std::endl;
#endif

  int numThreads = std::min(32, udata->num_cells);
  int numBlocks = static_cast<int>(ceil(((double) udata->num_cells)/((double) numThreads)));
  fun_csr_jac_kernel<<<numBlocks, numThreads>>>(t, y_d, fy_d, csr_sys->d_csr_values,
						user_data,
						csr_sys->size_per_subsystem,
						csr_sys->csr_number_nonzero,
						csr_sys->number_subsystems);
  cuda_status = cudaDeviceSynchronize();

#ifdef PRINT_DEBUG
  std::cout << "Got CUDA Synchronize return message of: ";
  std::cout << cudaGetErrorString(cuda_status) << std::endl;
#endif

#ifdef PRINT_DEBUG
  // debugging jac
  realtype jac_host[csr_sys->csr_number_nonzero * 2];

  // copy first and second system rhs to host
  cuda_status = cudaMemcpy(&jac_host[0], csr_sys->d_csr_values, sizeof(realtype) * csr_sys->csr_number_nonzero * 2, cudaMemcpyDeviceToHost);
  assert(cuda_status == cudaSuccess);
  
  // print out the first system jac
  std::cout << "first system jac = " << std::endl;
  for (int i = 0; i < csr_sys->csr_number_nonzero; i++)
    std::cout << jac_host[i] << " ";
  std::cout << std::endl;

  // print out the second system jac
  std::cout << "second system jac = " << std::endl;
  for (int i = csr_sys->csr_number_nonzero; i < 2*csr_sys->csr_number_nonzero; i++)
    std::cout << jac_host[i] << " ";
  std::cout << std::endl;
#endif
  
  assert(cuda_status == cudaSuccess);
  return 0;
}


__global__ static void fun_csr_jac_kernel(realtype t, realtype* y, realtype* fy,
					  realtype* csr_jac, void* user_data,
					  const int size, const int nnz, const int nbatched)
{
  UserData udata = static_cast<CVodeUserData*>(user_data);
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < udata->num_cells) {
    int offset = tid * udata->num_eqs_per_cell;
    int jac_offset = tid * udata->num_eqs_per_cell * udata->num_eqs_per_cell;
    int rpar_offset = tid * udata->num_rpar_per_cell;

    int csr_jac_offset = tid * nnz;
    realtype* csr_jac_cell = csr_jac + csr_jac_offset;

    sk_dense_jac_device(&t, &y[offset], &udata->dense_jacobians[jac_offset],
			&udata->rpar[rpar_offset]);

    sk_fill_csr_jac_device(&udata->dense_jacobians[jac_offset], csr_jac_cell);
  }  
}
