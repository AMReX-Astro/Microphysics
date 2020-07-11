#ifndef TEST_EOS_F_H_
#define TEST_EOS_F_H_

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

  void get_eos_len(int* eos_len);

  void get_eos_name(char* eos_string[]);

  void init_unit_test(const int* name, const int* namlen); 

  void do_eos(const int* lo, const int* hi,
              amrex::Real dlogrho, amrex::Real dlogT, amrex::Real dmetal,
              amrex::Real* state, const int* s_lo, const int* s_hi);

#ifdef __cplusplus
}
#endif

#endif
