# jac_cell

For a single thermodynamic condition, evaluate the Jacobian using both
the analytic version from the network and the finite-difference
numerical approximation.

Note: most analytic Jacobians do not consider the derivative of the
screening factor with respect to composition.  To get a more accurate
comparison of how the analytic and numerical terms compare, considering
just the reactive part, build with `SCREEN_METHOD=null`.

