This is a version of the Powell's hybid method.  It is a C++ port of
the Minpack hybrj algorithm.  The main change done, as compared to
that version, is that we assume that the matrix is square, n x n.

It also allows one to pass data to the function and Jacobian, of
arbitrary type.
