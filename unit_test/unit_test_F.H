#ifndef UNIT_TEST_F_H
#define UNIT_TEST_F_H

#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif
  void init_unit_test(const int* name, const int* namlen);

#ifdef __cplusplus
}
#endif

#endif
