#ifndef TEST_EOS_F_H_
#define TEST_EOS_F_H_

#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif

  void init_unit_test(const int* name, const int* namlen); 

  void do_f90_parameters();

#ifdef __cplusplus
}
#endif

#endif
