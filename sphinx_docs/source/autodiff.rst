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

Derivatives of single-variable functions
========================================

To take the derivative of some computation ``f(x)``, ``x``
must be an ``autodiff::dual``, and has to be seeded with
``autodiff::seed()`` before the function is called:

.. code-block:: c++

   autodiff::dual x = 3.14_rt;
   autodiff::seed(x);
   autodiff::dual result = f(x);

We can then use ``autodiff::val(result)`` or
``static_cast<amrex::Real>(result)`` to extract the function value, and
``autodiff::derivative(result)`` to get the derivative with respect to
x. The ``static_cast`` version has the advantage of working on normal
numbers as well.

Most functions can be updated to support autodiff by adding a template
parameter for the numeric type (the current code calls it ``dual_t``).
This should be used for any values that depend on the variables we're
differentiating with respect to.  Calls to functions from ``<cmath>`` as
well as ``amrex::min`` and ``amrex::max`` can be replaced with ones in
the ``admath`` namespace.  This namespace also exports the original
functions, so they work fine on normal numeric types too.


Derivatives of multi-variable functions
=======================================

.. code-block:: c++

    #include <AMReX_REAL.H>
    #include <microphysics_autodiff.H>

    using namespace amrex::literals;

    template <typename dual_t>
    dual_t f(const dual_t& x, const dual_t& y) {
        return y * admath::sin(x) + admath::exp(y);
    }

    int main() {
        using dual_t = autodiff::dual_array<1, 2>;
        amrex::Real x = 2.41;
        amrex::Real y = 0.38;
        dual_t result;
        {
            dual_t x_dual = x;
            dual_t y_dual = y;
            // seed the inputs
            x_dual.grad(1) = 1.0;
            y_dual.grad(2) = 1.0;
            // compute the function and derivatives
            result = f(x_dual, y_dual);
        }
        amrex::Real dfdx = result.grad(1);
        amrex::Real dfdy = result.grad(2);
        std::cout << "f(" << x << ", " << y << ") = " << autodiff::val(result) << "\n";
        std::cout << "df/dx = " << dfdx << "\n";
        std::cout << "df/dy = " << dfdy << "\n";
    }
