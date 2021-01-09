#ifndef TEST_EOS_F_H_
#define TEST_EOS_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif
  void init_variables_F();

  void init_runtime_parameters(const int* name, const int* namlen); 

  void do_eos(const int* lo, const int* hi,
              amrex::Real dlogrho, amrex::Real dlogT, amrex::Real dmetal,
              amrex::Real* state, const int* s_lo, const int* s_hi);

#ifdef __cplusplus
}
#endif

#endif
