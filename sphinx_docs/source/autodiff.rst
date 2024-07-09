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

Most functions can be updated to support ``autodiff`` by adding a
template parameter for the numeric type (the current code calls it
``dual_t``).  This should be used for any values that depend on the
variables we're differentiating with respect to.  Calls to functions
from ``<cmath>`` as well as ``amrex::min`` and ``amrex::max`` can be
replaced with ones in the ``admath`` namespace.  This namespace also
exports the original functions, so they work fine on normal numeric
types too.

To manually check whether a type is a dual number or not, use
``autodiff::detail::isDual<dual_t>``.

Derivatives of single-variable functions
========================================

To take the derivative of some computation ``f(x)``, ``x`` must be an
``autodiff::dual``, and has to be seeded with ``autodiff::seed()``
before the function is called:

.. code-block:: c++

   autodiff::dual x = 3.14_rt;
   autodiff::seed(x);
   autodiff::dual result = f(x);

We can then use ``autodiff::val(result)`` or
``static_cast<amrex::Real>(result)`` to extract the function value, and
``autodiff::derivative(result)`` to get the derivative with respect to
x.  The ``static_cast`` version has the advantage of working on normal
numbers as well.

``autodiff::seed(x)`` sets the derivative term of ``x`` to 1 (it is equivalent
to ``x.grad = 1.0``), which effectively tells the code that
:math:`\frac{\partial x}{\partial x} = 1`. This propagates through any
operations that ``f(x)`` does, and we end up with :math:`\frac{\partial
f(x)}{\partial x}` in the derivative term of ``result``.


Derivatives of multi-variable functions
=======================================

It is possible to calculate derivatives with respect to several input
variables in a single pass, by using a (mathematical) vector of
derivative terms instead of a single number.  The interface is very
similar to the single-variable case, with slight modifications: for
``N`` input variables, we use ``autodiff::dual_array<1, N>`` in place of
``autodiff::dual``, and pass the variables in order to
``autodiff::seed_array()``.  After calling the function, we can extract
the derivatives using structured binding declarations, as shown in the
following example program:

.. code-block:: c++

    #include <iostream>
    #include <AMReX_REAL.H>
    #include <microphysics_autodiff.H>

    using namespace amrex::literals;

    template <typename dual_t>
    dual_t f(const dual_t& x, const dual_t& y) {
        return y * admath::sin(x) + admath::exp(y);
    }

    int main() {
        using dual_t = autodiff::dual_array<1, 2>;
        dual_t result;
        dual_t x_dual = 2.41, y_dual = 0.38;
        // seed the inputs
        autodiff::seed_array(x_dual, y_dual);
        // compute the function and both derivatives in a single pass
        dual_t result = f(x_dual, y_dual);

        auto [dfdx, dfdy] = autodiff::derivative(result);
        std::cout << "f(" << x << ", " << y << ") = " << autodiff::val(result) << "\n";
        std::cout << "df/dx = " << dfdx << "\n";
        std::cout << "df/dy = " << dfdy << "\n";
    }
