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

   \frac{dY(A)}{dt} = -n_A \, \rho^{n_A + n_B + n_C - 1} \, Y(A)^{n_A} \, Y(B)^{n_B} \, Y(C)^{n_C} \frac{N_A \langle \sigma v \rangle_{ABC,DEF}}{n_A! n_B! n_C!}

where :math:`N_A \langle \sigma v \rangle_{ABC,DEF}` is the temperature
portion of the reaction rate (cross section averaged over the velocity
distribution), provide by nuclear reaction rate tables or fits from the nuclear
experimental community.

For approximate networks, it may be the case that the exponent on
:math:`Y` is not the same as the number of nuclei involved, e.g.,
:math:`n_A`, so we will carry the exponents on :math:`Y(A)`,
..., :math:`Y(F)` as :math:`a`, :math:`b`, :math:`c`, :math:`d`,
:math:`e`, :math:`f`, and instead write this as:

.. math::

   \frac{dY(A)}{dt} = -n_A \, \rho^{a + b + c - 1} \, Y(A)^a \, Y(B)^b \, Y(C)^c \frac{N_A \langle \sigma v \rangle_{ABC,DEF}}{a! b! c!}


Metadata
========

With the above formulation, adding a rate to the network means
supplying the nuclei (:math:`A`, ...), coefficients (:math:`n_A`,
...), and the exponents (:math:`a`, ...).

The network provides a function that returns a ``rhs_t`` given a rate
(passed in as an integer from an ``enum`` of all the rates that makeup
the network).  The function signature is:

.. code:: c++

   AMREX_GPU_HOST_DEVICE AMREX_INLINE
   constexpr rhs_t rhs_data (int rate)


For example, consider the reaction :math:`\isotm{He}{4} + \isotm{C}{12} \rightarrow \isotm{O}{16} + \gamma`.  The metadata for this is initialized as:

.. code:: c++

   rhs_t data{};

   data.species_A = C12;
   data.species_B = He4;
   data.species_D = O16;

   data.number_A = 1;
   data.number_B = 1;
   data.number_D = 1;

   data.exponent_A = 1;
   data.exponent_B = 1;
   data.exponent_D = 1;

There are some additional fields in ``rhs_t`` that can be used in
special cases (for approximate nets).


Loop over Rates
===============

The main logic for constructing RHS is contained in ``Microphysics/networks/rhs.H``

The basic flow 

pairing terms


Jacobian
========

With the same rate infrastructure, we are able to provide an analytic
Jacobian for our reaction networks.

.. note::

   We do one approximation to the species derivatives in the Jacobian.
   Some approximate networks have compound rates where :math:`N_A
   \langle \sigma v \rangle` can depend on composition, :math:`Y`.  We
   neglect this derivative in our Jacobian.

   Testing has shown that this does not greatly affect the performance
   of the network.
