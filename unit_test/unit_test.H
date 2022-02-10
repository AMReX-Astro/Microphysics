#ifndef UNIT_TEST_H
#define UNIT_TEST_H

void write_job_info(const std::string& dir);

AMREX_INLINE
void
init_unit_test() {

    // now sync with C++ and read in the C++ parameters
    init_extern_parameters();

}
#endif
