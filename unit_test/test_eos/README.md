# test_eos

This is the unit test for testing the equation of state.  This tests
both the Fortran and C++ implementations (according to the parameter
`do_cxx`).

The test sets up a cube of data (density on one axis, temperature on
another, and composition on the third) and calls the EOS in various
modes.

