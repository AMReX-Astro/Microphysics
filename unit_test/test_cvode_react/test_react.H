#ifndef TEST_EOS_H
#define TEST_EOS_H

#include "cvode/cvode.h"                  /* prototypes for CVODE fcts., consts.      */
#ifdef AMREX_USE_CUDA
#include "nvector/nvector_cuda.h"         /* access to CUDA N_Vector                */
#ifdef CVODE_CUSOLVER
#include "cusolver/cvode_cusolver_spqr.h"    /* access to cuSolver interface             */
#endif
#else
#include "nvector/nvector_serial.h"       /* access to serial N_Vector                */
#include "sunmatrix/sunmatrix_dense.h"    /* access to dense SUNMatrix                */
#endif
#include "sundials/sundials_types.h"      /* definition of realtype                   */

typedef struct CVodeUserData {
  int num_eqs;
  int num_cells;
  int num_eqs_per_cell;
  int num_rpar_per_cell;
  int num_sparse_jac_nonzero;

  int irp_dens;
  int irp_cv;
  int irp_cp;
  int irp_xn_not_evolved;
  int irp_abar;
  int irp_zbar;
  int irp_eta;
  int irp_ye;
  int irp_cs;
  int irp_dx;
  int irp_t_sound;
  int irp_y_init;
  int irp_self_heat;
  int irp_Told;
  int irp_dcvdt;
  int irp_dcpdt;
  int irp_t0;
  int irp_energy_offset;

  amrex::Real* rpar;

}* UserData;

void main_main();

void do_react(const int* lo, const int* hi,
	      amrex::Real* state, const int* s_lo, const int* s_hi,
	      const int ncomp, const amrex::Real dt);

void initialize_rpar_indices(UserData user_data, const int nspec_not_evolved,
			     const int num_eqs_per_cell);

void zero_rpar_data(UserData user_data, const int size_rpar);

void initialize_system(realtype* y, UserData udata);

void finalize_system(realtype* y, UserData udata);

static int fun_rhs(realtype t, N_Vector y, N_Vector ydot, void *user_data);

#ifdef AMREX_USE_CUDA
static void set_nvector_cuda(N_Vector vec, realtype* data, sunindextype size);

static void get_nvector_cuda(N_Vector vec, realtype* data, sunindextype size);

__global__ static void initialize_cell(realtype* y, UserData udata);

__global__ static void finalize_cell(realtype* y, UserData udata);

__global__ static void fun_rhs_kernel(realtype t, realtype* y, realtype* ydot,
				      void *user_data);

#ifdef CVODE_CUSOLVER

int fun_csr_jac(realtype t, N_Vector y, N_Vector fy,
		CV_cuSolver_csr_sys csr_sys, void* user_data);

__global__ static void fun_csr_jac_kernel(realtype t, realtype* y, realtype* fy,
					  realtype* csr_jac, void* user_data,
					  const int size, const int nnz, const int nbatched);

#else

static int fun_jac_times_vec(N_Vector v, N_Vector Jv, realtype t,
			     N_Vector y, N_Vector fy,
			     void *user_data, N_Vector tmp);

__global__ static void fun_jtv_kernel(realtype* v, realtype* Jv, realtype t,
				      realtype* y, realtype* fy,
				      void* user_data, realtype* tmp, realtype* csr_jac);

#endif

#else
static void set_nvector_serial(N_Vector vec, realtype* data, sunindextype size);

static void get_nvector_serial(N_Vector vec, realtype* data, sunindextype size);

static int fun_jac(realtype tn, N_Vector y, N_Vector fy, SUNMatrix J,
                   void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static void initialize_cell(realtype* y, UserData udata, int cell_id);

static void finalize_cell(realtype* y, UserData udata, int cell_id);

static void fun_rhs_kernel(realtype t, realtype* y, realtype* ydot,
			   void *user_data, int cell_id);
#endif

#endif
