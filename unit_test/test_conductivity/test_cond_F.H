#ifndef TEST_COND_F_H
#define TEST_COND_F_H

#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  void init_variables_F();

  void get_ncomp(int* ncomp);

  void get_name_len(int* nlen);

  void get_var_name(char* cstring[], int* idx);

  void get_cond_len(int* cond_len);

  void get_cond_name(char* cond_string[]);

  void init_unit_test(const int* name, const int* namlen); 

  void do_conductivity(const int* lo, const int* hi,
                       const Real dlogrho, const Real dlogT, const Real dmetal,
                       amrex::Real* state, const int* s_lo, const int* s_hi);

#ifdef __cplusplus
}
#endif

#endif
