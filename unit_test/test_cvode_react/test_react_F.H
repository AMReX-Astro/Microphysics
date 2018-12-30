#ifndef TEST_EOS_F_H_
#define TEST_EOS_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif
  void init_variables();

  void get_ncomp(int* ncomp);

  void get_number_equations(int* neqs, int* nspec_not_evolved);

  void get_species_index(int* index);

  void get_species_old_index(int* index);

  void get_temperature_index(int* index);

  void get_density_index(int* index);

  void get_omegadot_index(int* index);

  void get_density_hnuc_index(int* index);

  void get_name_len(int* nlen);

  void get_var_name(char* cstring[], int* idx);

  void init_unit_test(const int* name, const int* namlen);

  void init_state(const int* lo, const int* hi,
                  amrex::Real* state, const int* s_lo, const int* s_hi,
                  const int* npts);

  void get_state(amrex::Real* state, const int* s_lo, const int* s_hi, const int ncomp,
                 const int i, const int j, const int k, const int c,
                 amrex::Real* f);

  void set_state(amrex::Real* state, const int* s_lo, const int* s_hi, const int ncomp,
                 const int i, const int j, const int k, const int c,
                 const amrex::Real f);

  void sk_get_csr_jac_rowcols(int* csr_row_count, int* csr_col_index);

  void sk_get_sparse_jac_nnz(int* nnz);

  void sk_get_store_jacobian(int* sjac);

  void sk_get_num_steps_save_jacobian(int* nstepssavejac);

  void sk_renormalize_species(amrex::Real* y, amrex::Real* rpar);

  void sk_convert_to_molar(amrex::Real* y);

  void sk_get_num_rpar_comps(int* number_rpar_comps);

  void sk_get_nspec_evolve(int* number_species_evolved);

#ifdef AMREX_USE_CUDA
__device__  void sk_initialize_cell_device(amrex::Real* y, amrex::Real* rpar);

__device__  void sk_finalize_cell_device(amrex::Real* y, amrex::Real* rpar);

__device__  void sk_f_rhs_device(const amrex::Real* t, const amrex::Real* y,
				 amrex::Real* ydot, amrex::Real* rpar);

__device__  void sk_analytic_jac_device(const amrex::Real* t, const amrex::Real* y,
					amrex::Real* jac_mat, amrex::Real* rpar);

#ifdef CVODE_CUSOLVER

__device__  void sk_fill_csr_jac_device(const amrex::Real* jac_mat,
					amrex::Real* csr_jac_cell);

#else

__device__  void sk_jac_times_vec_device(const amrex::Real* jac_mat, const amrex::Real* vec,
					 amrex::Real* jtv);

#endif

#else
void sk_initialize_cell(amrex::Real* y, amrex::Real* rpar);

void sk_finalize_cell(amrex::Real* y, amrex::Real* rpar);

void sk_f_rhs(const amrex::Real* t, const amrex::Real* y,
	      amrex::Real* ydot, amrex::Real* rpar);

void sk_analytic_jac(const amrex::Real* t, const amrex::Real* y,
		     amrex::Real* jac_mat, amrex::Real* rpar);

void sk_full_jac(amrex::Real* y, amrex::Real* jac_mat, amrex::Real* rpar,
		 const int* neq_total, const int* ncells,
		 const int* neq_per_cell, const int* nrpar_per_cell);
#endif

#ifdef __cplusplus
}
#endif

#endif
