# `test_sdc_vode_rhs`

This is meant to just serve as a driver for the vode SDC RHS C++
functionality.  It operates on a single zone and just sets up a
burn_t, converts it to a dvode_t, calls the vode_rhs routine, and then
prints the output.
