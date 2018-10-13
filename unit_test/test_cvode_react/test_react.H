#ifndef TEST_EOS_H
#define TEST_EOS_H

#include "cvode/cvode.h"                  /* prototypes for CVODE fcts., consts.      */
#ifdef AMREX_USE_CUDA
#include <AMReX_Device.H>
#include "nvector/nvector_cuda.h"         /* access to CUDA N_Vector                */
#ifdef CVODE_CUSOLVER
#include "cusolver/cvode_cusolver_spqr.h"    /* access to cuSolver interface             */
#endif
#else
#include "nvector/nvector_serial.h"       /* access to serial N_Vector                */
#include "sunmatrix/sunmatrix_dense.h"    /* access to dense SUNMatrix                */
#endif
#include "sundials/sundials_types.h"      /* definition of realtype                   */

class CVodeUserData {

public:

  const int num_cells;
  const int num_eqs;
  const int num_eqs_per_cell;
  const int num_rpar;  
  const int num_rpar_per_cell;
  const int num_sparse_jac_nonzero;
  const int num_spec_not_evolved;

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

  CVodeUserData() : num_cells(0),
		    num_eqs_per_cell(0),
		    num_eqs(0),
		    num_rpar_per_cell(0),
		    num_rpar(0),
		    num_sparse_jac_nonzero(0),
		    num_spec_not_evolved(0)
  {
    rpar = NULL;
  }

  CVodeUserData(int n_eqs, int n_cells, int n_eqs_per_cell,
		int n_rpar_per_cell, int n_sparse_jac_nonzero,
		int n_spec_not_evolved) : num_cells(n_cells),
					  num_eqs_per_cell(n_eqs_per_cell),					  
					  num_eqs(n_eqs_per_cell*n_cells),
					  num_rpar_per_cell(n_rpar_per_cell),					  
					  num_rpar(n_rpar_per_cell*n_cells),
					  num_sparse_jac_nonzero(n_sparse_jac_nonzero),
					  num_spec_not_evolved(n_spec_not_evolved)
  {
    cudaError_t cuda_status = cudaSuccess;
    cuda_status = cudaMallocManaged(&rpar, num_rpar * sizeof(amrex::Real));
    assert(cuda_status == cudaSuccess);

    const int device = amrex::Device::deviceId();
    amrex::Device::mem_advise_set_preferred(rpar, num_rpar * sizeof(amrex::Real), device);

    initialize_rpar_indices();

    zero_rpar_data();
  }

  ~CVodeUserData()
  {
    cudaError_t cuda_status = cudaSuccess;    

    if (rpar != NULL) {
      cuda_status = cudaFree(rpar);
      assert(cuda_status == cudaSuccess);
    }
  }

  void initialize_rpar_indices()
  {
    int i = 0;
    irp_dens = i; i++;
    irp_cv = i; i++;
    irp_cp = i; i++;
    irp_xn_not_evolved = i; i+=num_spec_not_evolved;
    irp_abar = i; i++;
    irp_zbar = i; i++;
    irp_eta = i; i++;
    irp_ye = i; i++;
    irp_cs = i; i++;
    irp_dx = i; i++;
    irp_t_sound = i; i++;
    irp_y_init = i; i+=num_eqs_per_cell;
    irp_self_heat = i; i++;
    irp_Told = i; i++;
    irp_dcvdt = i; i++;
    irp_dcpdt = i; i++;
    irp_t0 = i; i++;
    irp_energy_offset = i;
  }

  void zero_rpar_data()
  {
    for (int i = 0; i < num_rpar; i++) {
      rpar[i] = 0.0;
    }
  }
};

void main_main();

void do_react(const int* lo, const int* hi,
	      amrex::Real* state, const int* s_lo, const int* s_hi,
	      const int ncomp, const amrex::Real dt,
	      long* n_rhs, long* n_linsetup);

void initialize_system(realtype* y, CVodeUserData* udata);

void finalize_system(realtype* y, CVodeUserData* udata);

static int fun_rhs(realtype t, N_Vector y, N_Vector ydot, void *user_data);

#ifdef AMREX_USE_CUDA
static void set_nvector_cuda(N_Vector vec, realtype* data, sunindextype size);

static void get_nvector_cuda(N_Vector vec, realtype* data, sunindextype size);

__global__ static void initialize_cell(realtype* y, CVodeUserData* udata);

__global__ static void finalize_cell(realtype* y, CVodeUserData* udata);

__global__ static void fun_rhs_kernel(realtype t, realtype* y, realtype* ydot,
				      CVodeUserData* udata);

#ifdef CVODE_CUSOLVER

int fun_csr_jac(realtype t, N_Vector y, N_Vector fy,
		CV_cuSolver_csr_sys csr_sys, void* user_data);

__global__ static void fun_csr_jac_kernel(realtype t, realtype* y, realtype* fy,
					  realtype* csr_jac, CVodeUserData* udata);
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

static void initialize_cell(realtype* y, CVodeUserData* udata, int cell_id);

static void finalize_cell(realtype* y, CVodeUserData* udata, int cell_id);

static void fun_rhs_kernel(realtype t, realtype* y, realtype* ydot,
			   void *user_data, int cell_id);
#endif

#endif
