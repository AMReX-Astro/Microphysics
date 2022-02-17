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
special cases (e.g., for approximate nets):

* ``forward_branching_ratio``, ``reverse_branching_ratio`` :

  Some reactions have multiple possible outcomes (or branches).  For example,
  in ``aprox13``, we have:

  .. math::

     \isotm{C}{12} + \isotm{O}{16} \rightarrow \left \{
       \begin{array}{c} \isotm{Mg}{24} + \isotm{He}{4} \\
                        \isotm{Si}{28} + \gamma \end{array} \right .

  To capture this, we would include this as:

  .. code:: c++

      case C12_O16_to_Mg24_He4:
          data.species_A = C12;
          data.species_B = O16;
          data.species_D = Mg24;
          data.species_E = He4;

          data.number_A = 1;
          data.number_B = 1;
          data.number_D = 1;
          data.number_E = 1;

          data.exponent_A = 1;
          data.exponent_B = 1;
          data.exponent_D = 1;
          data.exponent_E = 1;

          data.forward_branching_ratio = 0.5_rt;
          break;

      case C12_O16_to_Si28:
          data.species_A = C12;
          data.species_B = O16;
          data.species_D = Si28;

          data.number_A = 1;
          data.number_B = 1;
          data.number_D = 1;

          data.exponent_A = 1;
          data.exponent_B = 1;
          data.exponent_D = 1;

          data.forward_branching_ratio = 0.5_rt;
          break;

  This indicates that each branch happens 50% of the time.


* ``apply_identical_particle_factor`` : 

   Normally for rates involving identical nuclei, we divide
   the rate by a factor (:math:`n!`, where `n` is the number of the same nuclei participating).  This
   avoids double-counting.

   For some approximate networks, we want to skip this, since although
   the net reaction appears to have identical particles, it
   participates via a chain that does not need the identical particle
   factor.

   An example of this (from ``aprox19``) is the rate ``P_P_N_N_to_He4``, which represents

   .. math::

      p + p + n + n \rightarrow \isotm{He}{4} + 3 \gamma

   In the approximation used in this network, this rate proceeds as the sequence:

   .. math::

      p(n,\gamma) d(n,\gamma) \isotm{He}{3} (p, \gamma) \isotm{He}{4}

   or

   .. math::

      p(n,\gamma) d(p,\gamma) \isotm{He}{3} (n, \gamma) \isotm{He}{4}

   and none of the reactions in the sequence involve like-nuclei fusing.

* ``rate_can_be_tabulated`` :

  To save on computation, we can create a table of reaction rates
  by evaluating over a grid of temperature and then interpolating
  in temperature as needed.  This operates only on the 
  :math:`N_A \langle \sigma v \rangle` portion of the rate.

  Some rates are more complex than fits into the rate tabulation
  scheme, and therefore we turn off the ability to tabulate by
  setting ``rate_can_be_tabulated = 0``.  For instance, this applies
  to weak rates and any rate that is actually a compound rate build
  up via a combination of intermediate rates.

* ``additional_reaction_1``, ``additional_reaction_2``, ``additional_reaction_3`` :


* ``screen_forward_reaction``, ``screen_reverse_reaction`` :

  These fields indicate whether to compute and apply the screening
  factors to the reaction rates.  Usually we will do this on all
  rates, but sometimes if a rate involves additional rates in
  a sequence, the screening is instead applied to those rates.

.. note::

   Every intermediate reaction is some other reaction's additional
   reaction, but not necessarily vice versa.


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
