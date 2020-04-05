#ifndef TEST_APROX_RATES_F_H_
#define TEST_APROX_RATES_F_H_

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
