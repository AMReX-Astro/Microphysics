#ifndef UNIT_TEST_F_H
#define UNIT_TEST_F_H

#include <AMReX_BLFort.H>

#ifdef __cplusplus
#include <AMReX.H>
extern "C"
{
#endif
    void init_unit_test_fortran(const int* name, const int* namlen);

    void init_microphysics_fortran();

#ifdef __cplusplus
}
#endif

#endif
