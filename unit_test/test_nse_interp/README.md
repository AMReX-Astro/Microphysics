# test_nse_interp

This is a simple test of the NSE cubic interpolation.

It first checks that the functions to find the closest table entry
to our input data work (e.g. `nse_get_logrho_index()`) as expected.

Then it does 1-d interpolation in each direction (rho, T, Ye)
to make sure that the interpolation there works as expected.
This is done just for Abar.

Finally, it calls the full interface that does tricubic interpolation
and prints out the interpolated state.

This is for the tabular NSE: `USE_NSE_TABLE=TRUE`
