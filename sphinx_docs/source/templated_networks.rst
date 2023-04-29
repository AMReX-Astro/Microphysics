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

The properties of a templated network are contained in the network's
``actual_network.H`` file.

Rate ``enum``
=============

A network using the templated RHS construction needs to provide an
``enum`` called ``Rates::NetworkRates`` that lists all of the reaction
rates in the network.  Only the forward rates need to be listed.  The
reverse rates will use the same rate index.

.. note::

   For some of the networks, we may not actually want to model the
   reverse rates, so the function that computes the reverse rate will
   simply return `0`.  But it is always carried in the reaction
   infrastructure.

Note: some of the reactions listed involve nuclei that are not present
in the actual network (but represented as ``__extra_`` nuclei.  We call
these reactions *intermediate reactions*.  These are used as steps in building
compound reaction sequences in approximate networks and are called upon via the
*additional reaction* mechanism provided by the metadata (described below).

.. note::

   Every intermediate reaction is some other reaction's additional
   reaction, but not necessarily vice versa.



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


At the core, this is just a ``switch`` statement on the rate index:

.. code:: c++

   rhs_t data{};

   switch rate {

     // fill the metadata of individual rates
     ...

   }

For example, consider the reaction :math:`\isotm{He}{4} + \isotm{C}{12} \rightarrow \isotm{O}{16} + \gamma`.  The metadata for this is initialized as:

.. code:: c++

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

  Consider burning :math:`\isotm{Mg}{24}` to :math:`\isotm{Si}{28}`.  We can imagine two
  sequences:

  .. math::

     \isotm{Mg}{24} (\alpha, \gamma) \isotm{Si}{28}

  or

  .. math::

     \isotm{Mg}{24} (\alpha, p) \isotm{Al}{27} (p, \gamma) \isotm{Si}{28}

  In the approximate networks, we combine these two reaction sequences
  into a single effective rate.  If we use the notation:

  .. math::

     \lambda_{\alpha\gamma} \rightarrow \isotm{Mg}{24} (\alpha, \gamma) \isotm{Si}{28}

  .. math::

     \lambda_{\alpha p} \rightarrow \isotm{Mg}{24} (\alpha, p) \isotm{Al}{27}

  .. math::

     \lambda_{p\gamma} \rightarrow \isotm{Al}{27} (p, \gamma) \isotm{Si}{28}

  and the reverse rate of :math:`\lambda_{\alpha p}` is noted as :math:`\lambda_{p\alpha}`.

  The effective rate combining these two channels is:

  .. math::

     (\lambda_{\alpha\gamma})_\mathrm{effective} =
         \lambda_{\alpha\gamma} + \lambda_{\alpha p} \left [ 1 - \frac{\lambda_{p\alpha}}{\lambda_{p\alpha} + \lambda_{p\gamma}} \right ]

  To capture this approximation with our rate metadata, we specify the
  additional rates that are needed, and then in the loop over rates
  that follows their evaluation, we will arrange them into this
  approximation.

  The metadata for this reaction appears as:

  .. code:: c++

     case Mg24_He4_to_Si28:
         data.species_A = Mg24;
         data.species_B = He4;
         data.species_D = Si28;

         data.number_A = 1;
         data.number_B = 1;
         data.number_D = 1;

         data.exponent_A = 1;
         data.exponent_B = 1;
         data.exponent_D = 1;

         data.additional_reaction_1 = Mg24_He4_to_Al27_P;
         data.additional_reaction_2 = Al27_P_to_Si28;
         break;

* ``screen_forward_reaction``, ``screen_reverse_reaction`` :

  These fields indicate whether to compute and apply the screening
  factors to the reaction rates.  Usually we will do this on all
  charged-particle reactions.  For neutron or gamma capture reactions,
  screening should be manually disabled.

  Additionally, if a rate involves additional rates in a sequence, we
  sometimes disable screening, as the screening is instead applied to
  those additional rates.


Loop over Rates
===============

The main logic for constructing RHS is contained in ``Microphysics/networks/rhs.H`` as ``rhs()``:

.. code:: c++

   AMREX_GPU_HOST_DEVICE AMREX_INLINE
   void rhs (burn_t& state, Array1D<Real, 1, neqs>& ydot)

The basic flow is:

#. Convert the incoming mass fractions into molar fractions

#. Compute the common temperature factors (used by rates) and
   composition factors (used by screening)

#. Compute all of the intermediate rates.

   Since these rates are used multiple times, we compute them once and cache them.
   This is done solely for performance reasons, since computing the rates is expensive.

#. Loop over rates

   a. Compute the rate if it is not an intermediate rate, otherwise,
      get the rate from the cache.

   b. Find any additional rates needed if this rate is actual a rate
      sequence.  These additional rates are drawn from the cache.

   c. Post-process the rate.  This is only done for rates that have
      additional rates.  This is where we would combine the primary
      rate and additional rates into a single effective rate (like
      shown for the :math:`(\alpha,\gamma)` and
      :math:`(\alpha,p)(p,\gamma)` sequence above.

      If the rate has no additional rates, then this step is a no-op.

   d. Loop over species:

      * If the species is used by this rate, then add it to the
        ``ydot`` for the species.

        Note: we always add the forward and reverse rates paired
        together here.  This greatly reduces roundoff error when the
        rates should drive us toward equilibrium.

#. Evaluate the neutrino cooling term

#. Compute the energy term from the ``ydot`` s.

Note that all of construction logic is done using ``constexpr`` expressions
(including the for-loops), allowing all of this logic to be evaluated
at compile time.  This effectively means that the compiler writes out
the full RHS for the network, leaving only the rate evaluation for
runtime.



Jacobian
========

With the same rate infrastructure, we are able to provide an analytic
Jacobian for our reaction networks.  The overall structure is the same
as the ``rhs()`` function described above.

.. note::

   We do one approximation to the species derivatives in the Jacobian.
   Some approximate networks have compound rates where :math:`N_A
   \langle \sigma v \rangle` can depend on composition, :math:`Y`.  We
   neglect this derivative in our Jacobian.

   Testing has shown that this does not greatly affect the performance
   of the network.


Linear Algebra
==============

The VODE integrator needs routines to do LU factorization and back
substitution.  We build off of the linpack ``dgefa`` and ``dgesl``
routines, but because we know at compile time which Jacobian terms are
non-zero, we are able to use ``constexpr`` for-loops to only do the
calculations on non-zero elements.  This greatly reduces the amount of work
in the linear algebra.

Note:

* Currently we are still storing a dense Jacobian -- we just skip computation
  on the elements that are 0.

* These routines do not perform pivoting.  This does not seem to be an
  issue for the types of matrices we solve with reactions (since they are
  all of the form :math:`I - \tau J`, where :math:`tau` is the timestep).
