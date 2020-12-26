*****************
Reaction Networks
*****************

Network Requirements and Structure
==================================

A network both defines the composition advected by the hydro code as
well as describes the burning processes between those isotopes.
Evolving the species in a network requires an integrator. The design
of Microphysics decouples the integrator from the network, allowing
for the ability to swap integrators as desired. We discuss the
integrators in a later section.

At a minimum, a network needs to provide:

* ``nspec`` : the number of species in the network

* ``nrates`` : the number of reaction rates. This is used to
  allocate space in the ``rate_t`` type

* ``num_rate_groups`` : the number of components for each reaction
  rate we want to store in the ``rate_t`` type

* ``naux`` : the number of auxiliary quantities needed by the network (these are not evolved).

* ``aion(:)`` : the atomic weight (in atomic mass units) of the species

* ``zion(:)`` : the atomic number of the species

* ``spec_names(:)`` : a descriptive name of the species (e.g. "hydrogen-1")

* ``short_spec_names(:)`` : a shorten version of the species name (e.g. "H1")

* ``short_aux_names(:)`` : the names of the auxiliary quantities

* ``network_name`` : a descriptive name for the network

Most of these quantities are Fortran parameters.

.. note::

   A convention adopted in Microphysics is that each network is
   responsible for determining the energy release from a change in
   composition. Most networks will provide an array of the species
   binding energies and a routine to compute the energy yield from the
   reaction rates.

There are three primary files within each network directory.

* ``actual_network.f90``:

   This is the Fortran module actual_network with routines:

   * ``actual_network_init()``

   * ``actual_network_finalize()``

   This supplies the number and names of species and auxiliary
   variables, as well as other initializing data, such as their mass
   numbers, proton numbers, and binding energies. It needs to define
   the ``nspec`` and ``naux`` quantities as integer
   parameters. Finally, it
   must also define nrates, the number of reaction rates linking
   the isotopes in the network.

* ``actual_rhs.f90``:

   This is the Fortran module ``actual_rhs_module``, with routines:

   * ``actual_rhs_init()``

   * ``actual_rhs(state)``

   * ``actual_jac(state)``

   This supplies an interface for computing the right-hand-side of the
   network, the time-derivative of each species (and the temperature
   and nuclear energy release), as well as the analytic Jacobian.
   Both ``actual_rhs`` and ``actual_jac`` take a single argument,
   a burn_t state. They set the time-derivatives and Jacobian
   elements in this derived type directly.

   Note: some networks do not provide an analytic Jacobian and instead
   rely on the numerical difference-approximation to the Jacobian. In
   this case, the interface ``actual_jac`` is still needed to compile.

* ``actual_burner``:

   This is the Fortran module ``actual_burner_module``, with routines:

   * ``actual_burner_init()``

   * ``actual_burner(state_in, state_out, dt, time)``

   This contains the interface for doing an actual burn. Here,
   ``state_in`` and ``state_out`` are ``burn_t`` objects. In
   general, you will want to call integrator to use one of the
   pre-defined ODE integrators, but you could also write a custom
   integration here. This is covered in more detail in § \ `5 <#ch:networks:integrators>`__.

Notice that all three of these modules have initialization routines:

* ``actual_network_init()``

* ``actual_rhs_init()``

* ``actual_burner_init()``

These must be called upon initialization. These should be not called
within OpenMP parallel regions, because in general they will modify
shared module data.

Note, depending on the network, some of these may do nothing, but
these interfaces are all required for maximum flexibility.

Available Networks
==================

iso7, aprox13, aprox19, and aprox21
-----------------------------------

These are alpha-chains (with some other nuclei) from Frank Timmes.
These networks share common rates (from ``Microphysics/rates``),
plasma neutrino loses (from ``Microphysics/neutrinos``), and
electron screening (from ``Microphysics/screening``).

Energy generation.
^^^^^^^^^^^^^^^^^^

These networks store the total binding energy of the nucleus in MeV as
``bion(:)``. They then compute the mass of each nucleus in grams as:

.. math:: M_k = (A_k - Z_k) m_n + Z_k (m_p + m_e) - B_k

where :math:`m_n`, :math:`m_p`, and :math:`m_e` are the neutron, proton, and electron
masses, :math:`A_k` and :math:`Z_k` are the atomic weight and number, and :math:`B_k`
is the binding energy of the nucleus (converted to grams). :math:`M_k`
is stored as ``mion(:)`` in the network.

The energy release per gram is converted from the rates as:

.. math:: \edot = -N_A c^2 \sum_k \frac{dY_k}{dt} M_k - \edotnu

where :math:`N_A` is Avogadro’s number (to convert this to “per gram”)
and :math:`\edotnu` is the neutrino loss term.

breakout
--------

general_null
------------

``general_null`` is a bare interface for a nuclear reaction network;
no reactions are enabled, and no auxiliary variables are accepted. The
data in the Fortran module is defined at compile type by specifying an
inputs file. For example,
``Networks/general_null/triple_alpha_plus_o.net`` would describe the
triple-:math:`\alpha` reaction converting helium into carbon, as
well as oxygen and iron.

At compile time, the network module ``actual_network.f90``
is written using the python script ``write_network.py``
and the template ``network.template``. The make rule
for this is contained in ``Make.package``. The name of the inputs file
is specified by the variable ``GENERAL_NET_INPUTS``.

A version of this network comes with MAESTRO and CASTRO, so you do
not usually need to worry about the version in Microphysics.

ignition_chamulak
-----------------

This network was introduced in our paper on convection in white dwarfs
as a model of Type Ia supernovae :cite:`wdconvect`. It models
carbon burning in a regime appropriate for a simmering white dwarf,
and captures the effects of a much larger network by setting the ash
state and energetics to the values suggested in :cite:`chamulak:2008`.


.. _energy-generation.-1:

Energy generation.
^^^^^^^^^^^^^^^^^^

The binding energy, :math:`q`, in this
network is interpolated based on the density. It is stored as the
binding energy (ergs/g) *per nucleon*, with a sign convention that
binding energies are negative. The energy generation rate is then:

.. math:: \edot = q \frac{dX(\isotm{C}{12})}{dt} = q A_{\isotm{C}{12}} \frac{dY(\isotm{C}{12})}{dt}

(this is positive since both :math:`q` and :math:`dY/dt` are negative)

ignition_reaclib
----------------

ignition_simple
---------------

This is the original network used in our white dwarf convection
studies :cite:`lowMach4`. It includes a single-step
:math:`^{12}\mathrm{C}(^{12}\mathrm{C},\gamma)^{24}\mathrm{Mg}` reaction.
The carbon mass fraction equation appears as

.. math::

   \frac{D X(^{12}\mathrm{C})}{Dt} = - \frac{1}{12} \rho X(^{12}\mathrm{C})^2
       f_\mathrm{Coul} \left [N_A \left <\sigma v \right > \right]

where :math:`N_A \left <\sigma v\right>` is evaluated using the reaction
rate from (Caughlan and Fowler 1988). The Coulomb screening factor,
:math:`f_\mathrm{Coul}`, is evaluated using the general routine from the
Kepler stellar evolution code (Weaver 1978), which implements the work
of (Graboske 1973) for weak screening and the work of (Alastuey 1978
and Itoh 1979) for strong screening.

kpp
---

powerlaw
--------

This is a simple single-step reaction rate.
We will consider only two species, fuel, :math:`f`, and ash, :math:`a`, through
the reaction: :math:`f + f \rightarrow a + \gamma`. Baryon conservation
requres that :math:`A_f = A_a/2`, and charge conservation requires that :math:`Z_f
= Z_a/2`. We take
our reaction rate to be a powerlaw in temperature. The standard way
to write this is in terms of the number densities, in which case we
have

.. math:: \frac{d n_f}{d t} = -2\frac{d n_a}{d t} = -r

with

.. math:: r = r_0 n_X^2 \left( \frac{T}{T_0} \right )^\nu

Here, :math:`r_0` sets the overall rate, with units of
:math:`[\mathrm{cm^3~s^{-1}}]`, :math:`T_0` is a reference temperature scale, and
:math:`\nu` is the temperature exponent, which will play a role in setting
the reaction zone thickness. In terms of mass fractions, :math:`n_f = \rho
X_a / (A_a m_u)`, our rate equation is

.. math::

   \begin{align}
    \frac{dX_f}{dt} &= - \frac{r_0}{m_u} \rho X_f^2 \frac{1}{A_f} \left (\frac{T}{T_0}\right)^\nu \equiv \omegadot_f  \\
    \frac{dX_a}{dt} &= \frac{1}{2}\frac{r_0}{m_u} \rho X_f^2 \frac{A_a}{A_f^2} \left (\frac{T}{T_0}\right)^\nu = \frac{r_0}{m_u} \rho X_f^2 \frac{1}{A_f} \left (\frac{T}{T_0}\right)^\nu 
   \end{align}

We define a new rate constant, :math:`\rt` with units of :math:`[\mathrm{s^{-1}}]` as

.. math::

   \rt =  \begin{cases}
     \dfrac{r_0}{m_u A_f} \rho_0 & \text{if $T \ge T_a$} \\[1em]
     0                          & \text{if $T < T_a$}
    \end{cases}

where :math:`\rho_0` is a reference density and :math:`T_a` is an activation
temperature, and then our mass fraction equation is:

.. math:: \frac{dX_f}{dt} = -\rt X_f^2 \left (\frac{\rho}{\rho_0} \right ) \left ( \frac{T}{T_0}\right )^\nu

Finally, for the
energy generation, we take our reaction to release a specific energy,
:math:`[\mathrm{erg~g^{-1}}]`, of :math:`\qburn`, and our energy source is

.. math:: \edot = -\qburn \frac{dX_f}{dt}

There are a number of parameters we use to control the constants in
this network. This is one of the few networks that was designed
to work with gamma_law_general as the EOS.

rprox
-----

This network contains 10 species, approximating hot CNO,
triple-\ :math:`\alpha`, and rp-breakout burning up through :math:`^{56}\mathrm{Ni}`,
using the ideas from :cite:`wallacewoosley:1981`, but with modern
reaction rates from ReacLib :cite:`ReacLib` where available.
This network was used for the X-ray burst studies in
:cite:`xrb:II`, :cite:`xrb:III`, and more details are contained in those papers.

triple_alpha_plus_cago
----------------------

This is a 2 reaction network for helium burning, capturing the :math:`3`-:math:`\alpha`
reaction and :math:`\isotm{C}{12}(\alpha,\gamma)\isotm{O}{16}`. Additionally,
:math:`^{56}\mathrm{Fe}` is included as an inert species.


xrb_simple
----------

This is a simple 7 isotope network approximating the burning that
takes place in X-ray bursts (6 isotopes participate in reactions, one
additional, :math:`^{56}\mathrm{Fe}`, serves as an inert
composition). The 6 reactions modeled are:

* :math:`3\alpha + 2p \rightarrow \isotm{O}{14}` (limited by the 3-\ :math:`\alpha` rate)

* :math:`\isotm{O}{14} + \alpha \rightarrow \isotm{Ne}{18}` (limited by :math:`\isotm{O}{14}(\alpha,p)\isotm{F}{17}` rate)

* :math:`\isotm{O}{15} + \alpha + 6 p \rightarrow \isotm{Si}{25}` (limited by :math:`\isotm{O}{15}(\alpha,\gamma)\isotm{Ne}{19}` rate)

* :math:`\isotm{Ne}{18} + \alpha + 3p \rightarrow \isotm{Si}{25}` (limited by :math:`\isotm{Ne}{18}(\alpha,p)\isotm{Na}{21}` rate)

* :math:`\isotm{O}{14} + p \rightarrow \isotm{O}{15}` (limited by :math:`\isotm{O}{14}(e+\nu)\isotm{N}{14}` rate)

* :math:`\isotm{O}{15} + 3p \rightarrow \isotm{O}{14} + \alpha`  (limited by :math:`\isotm{O}{15}(e+\nu)\isotm{N}{15}` rate)

All reactions conserve mass. Where charge is not conserved, fast weak
interactions are assumed. Weak rates are trivial, fits to the 4
strong rates to a power law in :math:`T_9 \in [0.3, 1]`, linear in density.

subch
-----

This is a 10 isotope network including rates from reactions suggested
by Shen and Bildsten in their 2009 paper on helium burning on a white
dwarf :cite:`ShenBildsten`.  The reactions included in
this networks are as follows:

.. math::

   \begin{aligned}
       \isotm{He}{4} &\rightarrow  \isotm{C}{12} + 2\gamma \\
       \isotm{C}{12} + \isotm{He}{4} &\rightarrow \isotm{O}{16} + \gamma \\
       \isotm{N}{14} + \isotm{He}{4} &\rightarrow \isotm{F}{18} + \gamma \label{chemeq:1.1} \\
       \isotm{F}{18} + \isotm{He}{4} &\rightarrow \isotm{Ne}{21} +  \text{p} \label{chemeq:1.2} \\
       \isotm{C}{12} + p+ &\rightarrow \isotm{N}{13} + \gamma  \label{chemeq:2.1} \\
       \isotm{N}{13} + \isotm{He}{4} &\rightarrow \isotm{O}{16} + \text{p} \label{chemeq:2.2} \\
       \isotm{O}{16} + \isotm{He}{4} &\rightarrow \isotm{Ne}{20} + \gamma \\
       \isotm{C}{14} + \isotm{He}{4} &\rightarrow \isotm{O}{18} + \gamma \label{chemeq:3.2}
   \end{aligned}

The main reactions suggested by Shen and Bildsten were the :math:`\isotm{N}{14}(\alpha,\gamma)\isotm{F}{18}`,
leading into :math:`\isotm{F}{18}(\alpha,p)\isotm{Ne}{21}`,
:math:`\isotm{C}{12}(p,\gamma)\isotm{N}{13}` leading into :math:`\isotm{N}{13}(\alpha,p)\isotm{O}{16}`,
and :math:`\isotm{C}{14}(\alpha,\gamma)\isotm{O}{18}` :cite:`ShenBildsten`.
The rates of these reactions are shown in the figure below.
Notably, the reaction :math:`\isotm{N}{13}(\alpha,p)\isotm{O}{16}`, is high and may produce :math:`\isotm{O}{16}` more quickly than reactions involving only :math:`\isotm{He}{4}` and :math:`\isotm{C}{12}`,


.. figure:: subch.png
   :alt: pynucastro plot of the reaction rates of the subch network.
   :scale: 80%
   :align: center

   pynucastro plot of the reaction rates of the subch network.
