#ifndef TEST_SCREENING_F_H
#define TEST_SCREENING_F_H

#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  void init_variables_F();

  void do_screening(const int* lo, const int* hi,
                    const Real dlogrho, const Real dlogT, const Real dmetal,
                    amrex::Real* state, const int* s_lo, const int* s_hi);

#ifdef __cplusplus
}
#endif

#endif
