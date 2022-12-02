# test_ase

`test_ase` finds the NSE mass fractions of a given state.
The density, temperature, and composition are set in the
inputs file. Then we update the NSE mass fraction to
the current state to make sure we're in NSE. Then we
feed it into the function `in_nse` to check whether
we're in NSE or not.

Upon completion, the new state is printed to the screen.
And a statement is printed to see whether we're in NSE
or not.