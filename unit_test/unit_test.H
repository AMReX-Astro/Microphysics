#ifndef UNIT_TEST_H
#define UNIT_TEST_H

#ifdef MICROPHYSICS_FORT
#include <unit_test_F.H>
#include <extern_parameters_F.H>
#endif
void write_job_info(const std::string& dir);

AMREX_INLINE
void
init_unit_test(const int* name, const int* namelen) {

#ifdef MICROPHYSICS_FORT
    // first read in the Fortran parameters
    init_unit_test_fortran(name, namelen);
#endif

    // now sync with C++ and read in the C++ parameters
    init_extern_parameters();

#ifdef MICROPHYSICS_FORT
    // sync any C++ changes with Fortran
    update_fortran_extern_after_cxx();
#endif

#ifdef MICROPHYSICS_FORT
    // initialize Fortran microphysics
    init_microphysics_fortran();
#endif

}
#endif
