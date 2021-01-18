#ifndef TEST_APROX_RATES_F_H_
#define TEST_APROX_RATES_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif
  void init_variables_F();

  void do_rates(const int* lo, const int* hi,
                const Real dlogrho, const Real dlogT, const Real dmetal,
                amrex::Real* state, const int* s_lo, const int* s_hi);


#ifdef __cplusplus
}
#endif

#endif
