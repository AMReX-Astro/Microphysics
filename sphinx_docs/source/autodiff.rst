*************************
Automatic Differentiation
*************************

Support for automatic differentiation is provided by the ``autodiff``
library :cite:p:`autodiff`, included under
``Microphysics/util/autodiff``.  We use the forward mode ``dual``
implementation, which produces the derivative of each computation along
with its output value.  This results in largely the same arithmetic
operations as manually calculating the analytical derivative of each
intermediate step, but with much less code and fewer typos.  All the
machinery needed for use in Microphysics is located in
``Microphysics/util/microphysics_autodiff.H``.

To take the derivative of some computation ``f(x)``, ``x``
must be an ``autodiff::dual``, and has to be seeded with
``autodiff::seed()`` before the function is called:

.. code-block:: c++

   autodiff::dual x = 3.14_rt;
   autodiff::seed(x);
   autodiff::dual result = f(x);

We can then use ``autodiff::val(result)`` or
``static_cast<amrex::Real>(result)`` to extract the function value,
and ``autodiff::derivative(result)`` to get the derivative with respect
to x. Which has the advantage of working on both normal and dual
numbers.

Most functions can be updated to support autodiff by adding a template
parameter for the numeric type (the current code calls it ``dual_t``).
This should be used for any values that depend on the variables we're
differentiating with respect to.  Calls to functions from ``<cmath>`` as
well as ``amrex::min`` and ``amrex::max`` can be replaced with ones in
the ``admath`` namespace.  This namespace also exports the original
functions, so they work fine on normal numeric types too.
