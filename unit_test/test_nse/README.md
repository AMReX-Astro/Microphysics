# burn_cell

`burn_cell` integrates a network for a single set of initial
conditions.  The density, temperature, and composition are set in the
inputs file, as well as the maximum time to integrate.

Also calculates the nse mass fraction before and after integration.

Upon completion, the new state is printed to the screen.


## continuous integration

This is used with GitHub actions to do a test on pull requests.  The
code is built with the `subch2` network and run with `inputs_subch2`
and then diff is used to compared to the stored output in
`subch2_unit_test.out`.
