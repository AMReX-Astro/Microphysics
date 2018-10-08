#include <cvode/cvode.h>                  /* prototypes for CVODE fcts., consts.      */
#ifdef AMREX_USE_CUDA
#include <nvector/nvector_cuda.h>         /* access to CUDA N_Vector                  */
#include <sunlinsol/sunlinsol_spgmr.h>    /* access to SPGMR SUNLinearSolver          */
#include <cvode/cvode_spils.h>            /* access to CVSpils interface              */
#else
#include "nvector/nvector_serial.h"         /* access to serial N_Vector                */
#include "sunmatrix/sunmatrix_dense.h"        /* access to dense SUNMatrix                */
#include "sunlinsol/sunlinsol_dense.h"        /* access to dense SUNLinearSolver          */
#include "cvode/cvode_direct.h"           /* access to CVDls interface                */
#endif
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

  user_data = static_cast<UserData>(std::malloc(sizeof(struct CVodeUserData)));
  user_data->rpar = static_cast<amrex::Real*>(std::malloc(size_rpar * sizeof(amrex::Real)));

  user_data->num_eqs = size_state * neqs;
  user_data->num_cells = size_state;
  user_data->num_eqs_per_cell = neqs;
  user_data->num_rpar_per_cell = size_rpar_per_cell;

  initialize_rpar_indices(user_data, nspec_not_evolved, neqs);
  zero_rpar_data(user_data, size_rpar);

  realtype reltol=1.0e-6, time=0.0e0, tout;

  realtype abstol_values[size_flat];
  amrex::Real* state_y;
  state_y = static_cast<amrex::Real*>(std::malloc(size_flat * sizeof(amrex::Real)));

  N_Vector y = NULL, yout=NULL;
  N_Vector abstol = NULL;
  SUNMatrix Amat = NULL;
  SUNLinearSolver Linsol = NULL;
  void* cvode_mem = NULL;
  int flag;

  // Create NVectors
  y = N_VNew_Serial(size_flat);
  yout = N_VNew_Serial(size_flat);
  abstol = N_VNew_Serial(size_flat);

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

	// Energy absolute tolerance
	abstol_values[nzone*neqs + nspec_evolve + 1] = 1.0e-6;
	
	// Put temperature into integration vector
	get_state(state, s_lo, s_hi, ncomp, i, j, k, idx_temp,
		  &state_y[nzone*neqs + nspec_evolve]);

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
  set_nvector_serial(y, &state_y[0], size_flat);
  set_nvector_serial(abstol, &abstol_values[0], size_flat);

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
  Amat = SUNDenseMatrix(size_flat, size_flat);
  Linsol = SUNDenseLinearSolver(y, Amat);
  flag = CVDlsSetLinearSolver(cvode_mem, Linsol, Amat);
  if (flag != CV_SUCCESS) amrex::Abort("Failed to set linear solver");  
  flag = CVDlsSetJacFn(cvode_mem, fun_jac);
  if (flag != CV_SUCCESS) amrex::Abort("Failed to set jac function");

  // Do Integration
  time = time + static_cast<realtype>(dt);
  flag = CVode(cvode_mem, time, yout, &tout, CV_NORMAL);
  if (flag != CV_SUCCESS) amrex::Abort("Failed integration");

  // Get Final State
  get_nvector_serial(yout, &state_y[0], size_flat);

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
  N_VDestroy(y);
  N_VDestroy(yout);
  N_VDestroy(abstol);
  CVodeFree(&cvode_mem);
  SUNMatDestroy(Amat);
  SUNLinSolFree(Linsol);
  std::free(user_data->rpar);
  std::free(user_data);
  std::free(state_y);

}

/*
void do_rhs(const int* lo, const int* hi,
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
  amrex::Real* state_rhs;  
  cudaMallocManaged(&state_y, size_flat * sizeof(amrex::Real));
  cudaMallocManaged(&state_rhs, size_flat * sizeof(amrex::Real));  

  // Initialize y, abstol from flattened state
  int nzone = 0;
  for (int i=lo[0]; i<=hi[0]; i++) {
    for (int j=lo[1]; j<=hi[1]; j++) {
      for (int k=lo[2]; k<=hi[2]; k++) {
	// Put mass fractions into integration vector
	int scomp = 0;
        for (int n=idx_spec_old; n<idx_spec_old+nspec_evolve; n++) {
          get_state(state, s_lo, s_hi, &ncomp, &i, &j, &k, &n, &state_y[nzone*neqs + scomp]);
	  scomp++;
        }

	// Put temperature into integration vector
	get_state(state, s_lo, s_hi, &ncomp, &i, &j, &k, &idx_temp,
		  &state_y[nzone*neqs + scomp]);
	scomp++;

	// Initialize energy to 0, we'll get it from the EOS
	state_y[nzone*neqs + scomp] = 0.0e0;

	// Put density in user data
	get_state(state, s_lo, s_hi, &ncomp, &i, &j, &k, &idx_dens,
		  &user_data->rpar[nzone*user_data->num_rpar_per_cell + user_data->irp_dens]);

	// Put unevolved mass fractions in user data
	scomp=0;
	for (int n=idx_spec_old+nspec_evolve; n<idx_spec_old+nspec_evolve+nspec_not_evolved; n++) {
	  get_state(state, s_lo, s_hi, &ncomp, &i, &j, &k, &n,
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
  
  // Evaluate RHS
  rhs_only(&state_y[0], &state_rhs[0], user_data);

  // Save RHS
  nzone = 0;
  for (int i=lo[0]; i<=hi[0]; i++) {
    for (int j=lo[1]; j<=hi[1]; j++) {
      for (int k=lo[2]; k<=hi[2]; k++) {
	// Put RHS for evolved mass fractions into state
	int scomp = 0;
        for (int n=idx_spec; n<idx_spec+nspec_evolve; n++) {
          set_state(state, s_lo, s_hi, &ncomp, &i, &j, &k, &n, &state_rhs[nzone*neqs + scomp]);
	  scomp++;
        }

	// Put RHS for unevolved mass fractions into state
	amrex::Real actual_zero = 0.0;
	scomp=0;
	for (int n=idx_spec+nspec_evolve; n<idx_spec+nspec_evolve+nspec_not_evolved; n++) {
	  set_state(state, s_lo, s_hi, &ncomp, &i, &j, &k, &n, &actual_zero);
	  scomp++;
	}

	// Put 0 for omegadot into state
	amrex::Real xn_start, wscratch;
	scomp = 0;
        for (int n=idx_omegadot; n<idx_omegadot+nspec_evolve; n++) {
	  set_state(state, s_lo, s_hi, &ncomp, &i, &j, &k, &n, &actual_zero);
	  scomp++;
        }
	for (int n=idx_omegadot+nspec_evolve; n<idx_omegadot+nspec_evolve+nspec_not_evolved; n++) {
	  set_state(state, s_lo, s_hi, &ncomp, &i, &j, &k, &n, &actual_zero);
	}

	// Set rho*Hnuc to the RHS for energy
	set_state(state, s_lo, s_hi, &ncomp, &i, &j, &k, &idx_dens_hnuc, &state_rhs[nzone*neqs + nspec_evolve + 1]);

	nzone++;
      }
    }
  }

}
*/

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
  for (int i=0; i<udata->num_cells; i++) {
    initialize_cell(y, udata, i);
  }
}


static void initialize_cell(realtype* y, UserData udata, int cell_id)
{
  int offset = cell_id * udata->num_eqs_per_cell;
  int rpar_offset = cell_id * udata->num_rpar_per_cell;
  sk_initialize_cell(&y[offset], &udata->rpar[rpar_offset]);
}


void finalize_system(realtype* y, UserData udata)
{
  for (int i=0; i<udata->num_cells; i++) {
    finalize_cell(y, udata, i);
  }
}


static void finalize_cell(realtype* y, UserData udata, int cell_id)
{
  int offset = cell_id * udata->num_eqs_per_cell;
  int rpar_offset = cell_id * udata->num_rpar_per_cell;

  sk_finalize_cell(&y[offset], &udata->rpar[rpar_offset]);
}


static void set_nvector_serial(N_Vector vec, realtype* data, sunindextype size)
{
  for (sunindextype i = 0; i < size; i++) {
    NV_Ith_S(vec, i) = data[i];
  }
}


static void get_nvector_serial(N_Vector vec, realtype* data, sunindextype size)
{
  for (sunindextype i = 0; i < size; i++) {
    data[i] = NV_Ith_S(vec, i);
  }
}


static int rhs_only(amrex::Real* y, amrex::Real* ydot, UserData udata)
{
  amrex::Real t = 0.0;
  for (int i=0; i<udata->num_cells; i++) {
    fun_rhs_kernel(t, y, ydot, udata, i);
  }
  return 0;
}


static int fun_rhs(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  UserData udata = static_cast<CVodeUserData*>(user_data);
  amrex::Real* state_y = static_cast<amrex::Real*>(std::malloc(udata->num_cells * udata->num_eqs_per_cell * sizeof(amrex::Real)));
  amrex::Real* state_ydot = static_cast<amrex::Real*>(std::malloc(udata->num_cells * udata->num_eqs_per_cell * sizeof(amrex::Real)));
  get_nvector_serial(y, state_y, udata->num_cells * udata->num_eqs_per_cell);
  for (int i=0; i<udata->num_cells; i++) {
    fun_rhs_kernel(t, state_y, state_ydot, user_data, i);
  }
  set_nvector_serial(ydot, state_ydot, udata->num_cells * udata->num_eqs_per_cell);
  std::free(state_y);
  std::free(state_ydot);
  return 0;
}


static void fun_rhs_kernel(realtype t, realtype* y, realtype* ydot,
				      void *user_data, int cell_id)
{
  UserData udata = static_cast<CVodeUserData*>(user_data);
  int offset = cell_id * udata->num_eqs_per_cell;
  int rpar_offset = cell_id * udata->num_rpar_per_cell;

  sk_f_rhs(&t, &y[offset], &ydot[offset], &udata->rpar[rpar_offset]);
}


static int fun_jac(realtype tn, N_Vector y, N_Vector fy, SUNMatrix J,
                   void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  UserData udata = static_cast<CVodeUserData*>(user_data);  
  realtype* Jdata;
  amrex::Real* state_y = static_cast<amrex::Real*>(std::malloc(udata->num_cells * udata->num_eqs_per_cell * sizeof(amrex::Real)));
  get_nvector_serial(y, state_y, udata->num_cells * udata->num_eqs_per_cell);
  Jdata = SUNDenseMatrix_Data(J);

  sk_full_jac(state_y, Jdata, udata->rpar,
	      &udata->num_eqs, &udata->num_cells,
	      &udata->num_eqs_per_cell, &udata->num_rpar_per_cell);

  std::free(state_y);

  return 0;
}

