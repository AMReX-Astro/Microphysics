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

  void do_eos(const int* lo, const int* hi,
              amrex::Real* state, const int* s_lo, const int* s_hi,
              const int* npts);

#ifdef __cplusplus
}
#endif

#endif
