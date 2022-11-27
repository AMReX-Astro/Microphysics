# burn_cell

`burn_cell` integrates a network for a single set of initial
conditions.  The density, temperature, and composition are set in the
inputs file, as well as the maximum time to integrate.

Upon completion, the new state is printed to the screen.


## continuous integration

This is used with GitHub actions to do a test on pull requests.  The
following are the tests:

* `subch_approx` network:

  ```
  make NETWORK_DIR=subch_approx
  ./main3d.gnu.ex inputs_subch_approx > test.out
  diff test.out subch_approx_unit_test.out
  ```

* `ECSN` network:

  ```
  make NETWORK_DIR=ECSN
  ./main3d.gnu.ex inputs_ecsn > test.out
  diff test.out ecsn_unit_test.out
  ```
