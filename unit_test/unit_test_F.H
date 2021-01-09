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

    void runtime_pretty_print(int* jobinfo_file_name, const int* jobinfo_file_length);

#ifdef __cplusplus
}
#endif

#endif
