#ifndef UNIT_TEST_H
#define UNIT_TEST_H

void write_job_info(const std::string& dir);

AMREX_INLINE
extern_t
init_unit_test() {

    // now sync with C++ and read in the C++ parameters
    auto params = init_extern_parameters();
    return params;

}
#endif
