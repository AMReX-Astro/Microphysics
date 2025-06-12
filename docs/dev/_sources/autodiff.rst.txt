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
``number_t``).  This should be used for any values that depend on the
variables we're differentiating with respect to.  Calls to functions
from ``<cmath>`` as well as ``amrex::min``, ``amrex::max``, and
``amrex::Math::powi`` can be
replaced with ones in the ``admath`` namespace.  This namespace also
exports the original functions, so they work fine on normal numeric
types too.

.. tip::
   The ``autodiff`` library uses expression templates to optimize some
   expressions, which can cause compilation errors when trying to pass
   expressions using ``autodiff::dual`` to other templated functions.
   To avoid this, either use the corresponding function from ``admath``,
   or wrap the expression with ``autodiff::eval()`` before passing it to
   the templated function.

To manually check whether a type is a dual number or not, use
``autodiff::detail::isDual<number_t>``.

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
:math:`{\partial x}/{\partial x} = 1`. This propagates through any
operations that ``f(x)`` does, and we end up with :math:`{\partial
f(x)}/{\partial x}` in the derivative term of ``result``.


Derivatives of multi-variable functions
=======================================

It is possible to calculate derivatives with respect to several input
variables in a single pass, by using a (mathematical) vector of
derivative terms instead of a single number.  This is known as "vector
mode" automatic differentiation, and can save time when the function
evaluation code is much more expensive than the derivative calculations.

The interface is very
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

    template <typename number_t>
    number_t f(const number_t& x, const number_t& y) {
        return y * admath::sin(x) + admath::exp(y);
    }

    int main() {
        using number_t = autodiff::dual_array<1, 2>;
        number_t result;
        number_t x_dual = 2.41, y_dual = 0.38;
        // seed the inputs
        autodiff::seed_array(x_dual, y_dual);
        // compute the function and both derivatives in a single pass
        number_t result = f(x_dual, y_dual);

        auto [dfdx, dfdy] = autodiff::derivative(result);
        std::cout << "f(" << x << ", " << y << ") = " << autodiff::val(result) << "\n";
        std::cout << "df/dx = " << dfdx << "\n";
        std::cout << "df/dy = " << dfdy << "\n";
    }


Partial arrays
--------------

If an intermediate value only depends on some of the variables in an
``autodiff::dual_array``, we can save time and memory by not calculating
or storing the derivatives with respect to the other variables.  This is
done by assigning such values to a smaller ``dual_array`` with indices
that match the full ``dual_array``.  Any operations between these
partial arrays will skip the unused components at compile time.
Operations that combine different ranges will produce a promoted type
that holds all the components between the minimum and maximum indices
present.  For example, multiplying a ``dual_array<1, 1>`` and a
``dual_array<2, 3>`` will produce a ``dual_array<1, 3>``.  Partial
arrays can be implicitly converted into a wider array type, but
converting to a narrower array type must explicitly use the
``autodiff::narrow_array`` function (to prevent derivative components
from being accidentally discarded).  This function takes the target type
as a template parameter, and returns its argument converted to that
type.

Partial ``dual_arrays`` can be created with the
``autodiff::make_partial_arrays()`` helper function, which takes the
full array range as template parameters and the values for each variable
as arguments.  It returns a tuple of seeded partial arrays, which can be
unpacked using a structured binding declaration.  The following line
defines ``dual_array<1, 1> x = 1.0_rt`` and ``dual_array<2, 2> y =
2.0_rt`` and seeds both of their derivative terms.

.. code-block:: c++

    auto [x, y] = autodiff::make_partial_arrays<1, 2>(1.0_rt, 2.0_rt);

The recommended way to use partial arrays is to declare each
intermediate value with ``auto`` and wrap each expression with
``autodiff::eval()``, as shown in this example code (using the variables
declared above):

.. code-block:: c++

    using autodiff::eval;

    // this only computes the derivative terms with respect to x
    // (type is dual_array<1, 1>)
    auto x_squared = eval(x * x);

    // this only computes the derivative terms with respect to y
    // (type is dual_array<2, 2>)
    auto sin_2y = eval(admath::sin(2.0_rt * y));

    // partial arrays are promoted as needed by overloaded operators
    // (type is dual_array<1, 2>)
    auto z = eval(x_squared * sin_2y);

The following is equivalent to the code above, but requires more care to
keep the types and expressions in sync:

.. code-block:: c++

   using autodiff::dual_array;

    // this only computes the derivative terms with respect to x
    dual_array<1, 1> x_squared = x * x;

    // this only computes the derivative terms with respect to y
    dual_array<2, 2> sin_2y = admath::sin(2.0_rt * y);

    // partial arrays are promoted as needed by overloaded operators
    dual_array<1, 2> z = x_squared * sin_2y;
