#include <cvode/cvode.h>                  /* prototypes for CVODE fcts., consts.      */
#include <nvector/nvector_cuda.h>         /* access to CUDA N_Vector                  */
#include <sunlinsol/sunlinsol_spgmr.h>    /* access to SPGMR SUNLinearSolver          */
#include <cvode/cvode_spils.h>            /* access to CVSpils interface              */
#include <sundials/sundials_types.h>      /* definition of realtype                   */
#include <sundials/sundials_math.h>       /* contains the macros ABS, SUNSQR, and EXP */
#include <AMReX_MultiFab.H>
#include "test_react.H"
#include "test_react_F.H"
#include <iostream>

using namespace amrex;

void do_react(const int* lo, const int* hi,
	      amrex::Real* state, const int* s_lo, const int* s_hi,
	      const int ncomp, const amrex::Real dt,
	      long* n_rhs, long* n_jac, long* n_linsetup)
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

  CVodeUserData* user_data;
  cudaMallocManaged(&user_data, sizeof(CVodeUserData));

  realtype reltol=1.0e-6, time=0.0e0, tout;

  realtype abstol_values[size_flat];
  amrex::Real* state_y;
  cudaMallocManaged(&state_y, size_flat * sizeof(amrex::Real));

  N_Vector y = NULL, yout=NULL;
  N_Vector abstol = NULL;
  SUNLinearSolver Linsol = NULL;

  int jac_number_nonzero;

  sk_get_sparse_jac_nnz(&jac_number_nonzero);

  int csr_row_count[neqs+1];
  int csr_col_index[jac_number_nonzero];

  sk_get_csr_jac_rowcols(&csr_row_count[0], &csr_col_index[0]);

  new (user_data) CVodeUserData(size_flat, size_state, neqs,
				size_rpar_per_cell, jac_number_nonzero,
				nspec_not_evolved);

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
	  abstol_values[nzone*neqs + scomp] = 1.0e-12;
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

  // Initialize Linear Solver
  Linsol = SUNSPGMR(y, PREC_NONE, 0);
  flag = CVSpilsSetLinearSolver(cvode_mem, Linsol);
  if (flag != CV_SUCCESS) amrex::Abort("Failed to set linear solver");  
  flag = CVSpilsSetJacTimes(cvode_mem, NULL, fun_jac_times_vec);
  if (flag != CV_SUCCESS) amrex::Abort("Failed to set jac*vec function");

  // Do Integration
  time = time + static_cast<realtype>(dt);
  flag = CVode(cvode_mem, time, yout, &tout, CV_NORMAL);
  if (flag != CV_SUCCESS) amrex::Abort("Failed integration");

  flag = CVodeGetNumRhsEvals(cvode_mem, n_rhs);
  flag = CVSpilsGetNumJtimesEvals(cvode_mem, n_jac);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, n_linsetup);
  
  // Get Final State
  get_nvector_cuda(yout, &state_y[0], size_flat);

  // Finalize cell data to save
  finalize_system(&state_y[0], user_data);

  // Save Final State
  nzone = 0;
  for (int i=lo[0]; i<=hi[0]; i++) {
    for (int j=lo[1]; j<=hi[1]; j++) {
      for (int k=lo[2]; k<=hi[2]; k++) {
	// Put evolved mass fractions into state
	int scomp = 0;
        for (int n=idx_spec; n<idx_spec+nspec_evolve; n++) {
          set_state(state, s_lo, s_hi, ncomp, i, j, k, n, state_y[nzone*neqs + scomp]);
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
	    xn_final = state_y[nzone*neqs + scomp];
	  } else {
	    xn_final = user_data->rpar[nzone*user_data->num_rpar_per_cell + user_data->irp_xn_not_evolved + scomp - nspec_evolve];
	  }
	  wscratch = (xn_final - xn_start)/dt;
	  set_state(state, s_lo, s_hi, ncomp, i, j, k, n, wscratch);
	  scomp++;
        }

	// Set rho*Hnuc
	get_state(state, s_lo, s_hi, ncomp, i, j, k, idx_dens, &wscratch);
	wscratch = wscratch * state_y[nzone*neqs + nspec_evolve + 1]/dt;
	set_state(state, s_lo, s_hi, ncomp, i, j, k, idx_dens_hnuc, wscratch);

	nzone++;
      }
    }
  }

  // Free Memory
  user_data->~CVodeUserData();
  cudaFree(user_data);
  cudaFree(state_y);
  N_VDestroy(y);
  N_VDestroy(yout);
  N_VDestroy(abstol);
  CVodeFree(&cvode_mem);
  SUNLinSolFree(Linsol);

}


void initialize_system(realtype* y, CVodeUserData* udata)
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


__global__ static void initialize_cell(realtype* y, CVodeUserData* udata)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < udata->num_cells) {
    int offset = tid * udata->num_eqs_per_cell;
    int rpar_offset = tid * udata->num_rpar_per_cell;
    sk_initialize_cell_device(&y[offset], &udata->rpar[rpar_offset]);
  }
}


void finalize_system(realtype* y, CVodeUserData* udata)
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


__global__ static void finalize_cell(realtype* y, CVodeUserData* udata)
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
  CVodeUserData* udata = static_cast<CVodeUserData*>(user_data);
  int numThreads = std::min(32, udata->num_cells);
  int numBlocks = static_cast<int>(ceil(((double) udata->num_cells)/((double) numThreads)));
  fun_rhs_kernel<<<numBlocks, numThreads>>>(t, y_d, ydot_d, udata);

  cuda_status = cudaDeviceSynchronize();
#ifdef PRINT_DEBUG
  std::cout << "In fun_rhs, got cudaDeviceSynchronize error of: " << cudaGetErrorString(cuda_status) << std::endl;
#endif
  assert(cuda_status == cudaSuccess);
  
  return 0;
}


__global__ static void fun_rhs_kernel(realtype t, realtype* y, realtype* ydot,
				      CVodeUserData* udata)
{
  int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < udata->num_cells) {
    int offset = tid * udata->num_eqs_per_cell;
    int rpar_offset = tid * udata->num_rpar_per_cell;

    sk_f_rhs_device(&t, &y[offset], &ydot[offset], &udata->rpar[rpar_offset]);
  }
}


static int fun_jac_times_vec(N_Vector v, N_Vector Jv, realtype t,
			     N_Vector y, N_Vector fy,
			     void *user_data, N_Vector tmp)
{
  cudaError_t cuda_status = cudaSuccess;
  cuda_status = cudaGetLastError();
#ifdef PRINT_DEBUG
  std::cout << "In fun_jac_times_vec, got CUDA last error of: " << cudaGetErrorString(cuda_status) << std::endl;
#endif
  assert(cuda_status == cudaSuccess);  
  
  realtype* v_d   = N_VGetDeviceArrayPointer_Cuda(v);
  realtype* Jv_d  = N_VGetDeviceArrayPointer_Cuda(Jv);
  realtype* y_d   = N_VGetDeviceArrayPointer_Cuda(y);
  realtype* fy_d  = N_VGetDeviceArrayPointer_Cuda(fy);
  realtype* tmp_d = N_VGetDeviceArrayPointer_Cuda(tmp);

  CVodeUserData* udata = static_cast<CVodeUserData*>(user_data);
  int numThreads = std::min(32, udata->num_cells);
  int numBlocks = static_cast<int>(ceil(((double) udata->num_cells)/((double) numThreads)));

  // allocate space for the sparse jacobian evaluation
  realtype* csr_jac_d = NULL;
  cuda_status = cudaMalloc((void**) &csr_jac_d, udata->num_cells * udata->num_sparse_jac_nonzero * sizeof(realtype));
#ifdef PRINT_DEBUG
  std::cout << "In fun_jac_times_vec, got cudaDeviceSynchronize error of: " << cudaGetErrorString(cuda_status) << std::endl;
#endif
  assert(cuda_status == cudaSuccess);
  
  fun_jtv_kernel<<<numBlocks, numThreads>>>(v_d, Jv_d, t,
					    y_d, fy_d,
					    udata, tmp_d, csr_jac_d);
  
  cuda_status = cudaDeviceSynchronize();
#ifdef PRINT_DEBUG
  std::cout << "In fun_jac_times_vec, got cudaDeviceSynchronize error of: " << cudaGetErrorString(cuda_status) << std::endl;
#endif
  assert(cuda_status == cudaSuccess);

  // free space for the sparse jacobian evaluation
  cuda_status = cudaFree(csr_jac_d);
#ifdef PRINT_DEBUG
  std::cout << "In fun_jac_times_vec, got cudaDeviceSynchronize error of: " << cudaGetErrorString(cuda_status) << std::endl;
#endif
  assert(cuda_status == cudaSuccess);

  return 0;
}


__global__ static void fun_jtv_kernel(realtype* v, realtype* Jv, realtype t,
				      realtype* y, realtype* fy,
				      CVodeUserData* udata, realtype* tmp, realtype* csr_jac)
{
  const int tid = blockIdx.x * blockDim.x + threadIdx.x;
  if (tid < udata->num_cells) {
    const int offset = tid * udata->num_eqs_per_cell;
    const int csr_jac_offset = tid * udata->num_sparse_jac_nonzero;
    const int rpar_offset = tid * udata->num_rpar_per_cell;

    sk_analytic_jac_device(&t, &y[offset], &csr_jac[csr_jac_offset],
			   &udata->rpar[rpar_offset]);

    sk_jac_times_vec_device(&csr_jac[csr_jac_offset], &v[offset], &Jv[offset]);
  }
}
