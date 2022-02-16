*********************************
Templated Network Righthand Sides
*********************************

A network consists of evaluating the rates and screening, and then for
each rate, adding or subtracting the appropriate contribution to the
evolution equation for each nucleus it affects.  That second part is
something that can be automated, greatly simplifying the construction
of the righthand side function (and Jacobian) as well as enabling
optimizations.

The templated network righthand side functionality builds the
righthand side function at *compile time* using C++ templates.  For
each reaction, we need to provide only a small amount of metadata.
This greatly simplifies the networks.  This infrastructure works for
approximate networks as well as regular networks.

Metadata
========

Consider a reaction of the form:

.. math::

   n_A A + n_B B + n_C C \rightarrow n_D D + n_E E + n_F F

where :math:`A`, :math:`B`, :math:`C`, :math:`D`, :math:`E`, and
:math:`F` are nuclei and :math:`n_A`, :math:`n_B`, :math:`n_C`,
:math:`n_D`, :math:`n_E`, and :math:`n_F` are the number of each
nuclei that participate in the reacton (these can be zero).  This form
of a reaction covers all of the cases we use in Microphysics.

For this reaction, the evolution equation for :math:`A` would be:

.. math::

   \frac{dY(A)}{dt} = -n_A \, \rho^{n_A + n_B + n_C - 1} \, Y(A)^{n_A} \, Y(B)^{n_B} \, Y(C)^{n_C} \frac{\langle \sigma v \rangle_{ABC,DEF}}{n_A! n_B! n_C!}

where :math:`\langle \sigma v \rangle_{ABC,DEF}` is the temperature
portion of the reaction rate (cross section averaged over the velocity
distribution).

For approximate networks, it may be the case that the exponent on
:math:`Y` is not the same as the number of nuclei involved, e.g.,
:math:`n_A`, so we will carry the exponents on :math:`Y(A)`
... :math:`Y(F)` as :math:`a`, :math:`b`, :math:`c`, :math:`d`,
:math:`e`, :math:`f`, and instead write this as:

.. math::

   \frac{dY(A)}{dt} = -n_A \, \rho^{n_A + n_B + n_C - 1} \, Y(A)^a \, Y(B)^b \, Y(C)^c \frac{\langle \sigma v \rangle_{ABC,DEF}}{a! b! c!}




Loop over Rates
===============


Jacobian
========

