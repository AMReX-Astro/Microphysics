#ifndef UNIT_TEST_H
#define UNIT_TEST_H

#include <unit_test_F.H>
#include <extern_parameters_F.H>

void write_job_info(const std::string& dir);

AMREX_INLINE
void
init_unit_test(const int* name, const int* namelen) {

    // first read in the Fortran parameters
    init_unit_test_fortran(name, namelen);

    // now sync with C++ and read in the C++ parameters
    init_extern_parameters();

    // sync any C++ changes with Fortran
    update_fortran_extern_after_cxx();

    // initialize Fortran microphysics
    init_microphysics_fortran();

}
#endif
