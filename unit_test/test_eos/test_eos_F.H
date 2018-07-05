#ifndef TEST_EOS_F_H_
#define TEST_EOS_F_H_

#include <AMReX_BLFort.H>

#ifdef __cpluscplus
extern "C"
{
#endif
  void init_variables();

  void get_ncomp(int* ncomp);

  void do_eos(const int* lo, const int* hi,
              BL_FORT_FAB_ARG_3D(state));

#ifdef __cpluscplus
}
#endif

#endif
