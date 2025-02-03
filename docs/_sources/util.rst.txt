******************************
Helper Functions and Libraries
******************************

The ``util/`` directory contains a number of external libraries and simple
utilities that are used by the different components of Microphysics.

* ``approx_math/`` : these are a set of headers that implement
  approximations to ``atan()``, ``exp()``, ``log()``, and ``pow()``.
  These can be much faster than the C++ library versions, especially
  on GPUs.

* ``autodiff/`` : this is a clone of the C++ autodiff library from
  https://github.com/autodiff/autodiff

  The header ``microphysics_autodiff.H`` provides a set of interfaces
  for working with the AMReX datatypes and interfacing with the
  autodiff library.

* ``build_scripts/`` : a set of python scripts used during the build
  process to parse the runtime parameters.

* ``cj_detonation/`` : a simple routine to compute the Chapman-Jouguet
  detonation speed for one of our networks.

* ``esum.H`` : an implementation of the exact sum algorithm based on the
  msum algorithm by Raymond Hettinger.  It is generated automatically
  by the ``esum_cxx.py`` script and creates implementations for exact
  numbers of terms (``esum3()``, ``esum4()``, ...)

* ``gcem/`` : a templated math library that provides implementations of
  the standard library math functions that can be used in ``constexpr``
  expressions.  This is from https://github.com/kthohr/gcem

  Some of the constants are redefined in 64-bit floating point in
  ``microphysics_math.H`` to avoid ``long double`` issues on some
  architectures.

* ``hybrj/`` : a C++ port of the MINPACK hybrid Powell minimization function
  to zero a set of functions.

* ``linpack.H`` : a C++ port of the LINPACK ``dgesl`` and ``dgefa`` LU
  decomposition Gaussian elimination routines.

* ``microphysics_sort.H`` : a set of sorting routines for
  ``amrex::Array1D`` data.
