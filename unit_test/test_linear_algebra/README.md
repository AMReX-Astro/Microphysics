# `test_linear_algebra`

This test creates a simple diagonally-dominant matrix, A, and multiples by
a test x to get a righthand side vector, b.  It then uses the linear algebra
routines to solve Ax = b and compares to the original x.

This is done twice, once with the constexpr linear algebra routines in `rhs.H`
and then with the routines in `linpack.H`.
