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

* ``nspec_evolve`` : the number of species that are actually integrated
  in the network.

  Usually this is ``nspec``, but in general any ``nspec_evolve``
  :math:`\le` ``nspec`` is allowed. Those species not evolved are held
  constant in the integration.

  Note that the convention is that the first ``nspec_evolve`` out
  of the ``nspec`` species are the ones evolved.

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

**A convention adopted in Microphysics is that each network
is responsible for determining the energy release from a change
in composition**. Most networks will provide an array of the species
binding energies and a routine to compute the energy yield from
the reaction rates.

There are three primary files within each network directory.

* ``actual_network.f90``:

   This is the Fortran module actual_network with routines:

   * ``actual_network_init()``

   * ``actual_network_finalize()``

   This supplies the number and names of species and auxiliary
   variables, as well as other initializing data, such as their mass
   numbers, proton numbers, and binding energies. It needs to define
   the ``nspec`` and ``naux`` quantities as integer
   parameters. Additionally it must define ``nspec_evolve``, the
   number of species that are actually evolved during a burn; in most
   cases, this should have the same value as ``nspec``. Finally, it
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

aprox13, aprox19, and aprox21
-----------------------------

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

CONe2NSE
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
for this is contained in ``Make.package`` (for C++ AMReX) and
``GPackage.mak`` (for F90 AMReX). The name of the inputs file
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

This network has ``nspec`` = 3, but ``nspec_evolve`` = 1. Only a
single reaction is modeled, converting :math:`^{12}\mathrm{C}` into
“ash”.

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

iso7
----

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
    \frac{dX_f}{dt} &= - \frac{r_0}{m_u} \rho X_f^2 \frac{1}{A_f} \left (\frac{T}{T_0}\right)^\nu \equiv \omegadot_f \label{eq:Xf} \\
    \frac{dX_a}{dt} &= \frac{1}{2}\frac{r_0}{m_u} \rho X_f^2 \frac{A_a}{A_f^2} \left (\frac{T}{T_0}\right)^\nu = \frac{r_0}{m_u} \rho X_f^2 \frac{1}{A_f} \left (\frac{T}{T_0}\right)^\nu  \label{eq:Xa}
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

This network has ``nspec`` = 4, but ``nspec_evolve`` = 3.

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

The main reactions suggested by Shen and Bildsten were the reaction series of
chemical equation `[chemeq:1.1] <#chemeq:1.1>`__ leading into equation `[chemeq:1.2] <#chemeq:1.2>`__,
chemical equation `[chemeq:2.1] <#chemeq:2.1>`__ leading into equation `[chemeq:2.2] <#chemeq:2.2>`__,
and chemical equation `[chemeq:3.2] <#chemeq:3.2>`__ :cite:`ShenBildsten`.
The rates of these reactions are shown in the figure below.
Notably, the reaction rate of chemical equation `[chemeq:2.2] <#chemeq:2.2>`__ is high and may produce Oxygen-16 more quickly than reactions involving only Helium-4, and Carbon-12.

.. raw:: latex

   \centering

.. figure:: subch.png
   :alt: pynucastro plot of the reaction rates of the subch network.
   :scale: 80%
   :align: center

   pynucastro plot of the reaction rates of the subch network.

Reaction ODE System
===================

Note: the integration works on the state :math:`\rho`, :math:`T`, and :math:`X_k`, e.g., the
mass fractions, but the individual networks construct the rates in terms
of the molar fractions, :math:`Y_k`. The wrappers between the integrators and
network righthand side routines do the conversion of the state to mass
fractions for the integrator.

The equations we integrate to do a nuclear burn are:

.. math::
   \frac{dX_k}{dt} = \omegadot_k(\rho,X_k,T)
   :label: eq:spec_integrate

.. math::
   \frac{d\enuc}{dt} = f(\dot{X}_k)
   :label: eq:enuc_integrate

.. math::
   \frac{dT}{dt} =\frac{\edot}{c_x}.
   :label: eq:temp_integrate

Here, :math:`X_k` is the mass fraction of species :math:`k`, :math:`\enuc` is the specifc
nuclear energy created through reactions, :math:`T` is the
temperature [1]_ , and :math:`c_x` is the specific heat for the
fluid. The function :math:`f` provides the energy release based on the
instantaneous reaction terms, :math:`\dot{X}_k`. As noted in the previous
section, this is implemented in a network-specific manner.

In this system, :math:`\enuc` is not necessarily the total specific internal
energy, but rather just captures the energy release during the burn. In
this system, it acts as a diagnostic,

.. math:: \enuc = \int \edot dt

so we know how much energy was released (or required) over the burn.

While this is the most common way to construct the set of
burn equations, and is used in most of our production networks,
all of them are ultimately implemented by the network itself, which
can choose to disable the evolution of any of these equations by
setting the RHS to zero. The integration software provides some
helper routines that construct common RHS evaluations, like the RHS
of the temperature equation given :math:`\dot{e}`, but these calls
are always explicitly done by the individual networks rather than
being handled by the integration backend. This allows you to write a
new network that defines the RHS in whatever way you like.

Interfaces
----------

.. note::

   StarKiller integrates the reaction system in terms of mass fractions,
   :math:`X_k`, but most astrophysical networks use molar fractions,
   :math:`Y_k`.  As a result, we expect the networks to return the
   righthand side and Jacobians in terms of molar fractions.  The StarKiller
   routines will internally convert to mass fractions as needed for the
   integrators.

The righthand side of the network is implemented by
``actual_rhs()`` ``in actual_rhs.f90``, and appears as:

::

      subroutine actual_rhs(state)
        type (burn_t) :: state

All of the necessary integration data comes in through state, as:

* ``state % xn(:)`` : the ``nspec`` mass fractions (note: for
  the case that ``nspec_evolve`` < ``nspec``, an algebraic constraint
  may need to be enforced. See § \ `3.3 <#ch:networks:nspec_evolve>`__).

* ``state % e`` : the current value of the ODE system’s energy
  release, :math:`\enuc`—note: as discussed above, this is not
  necessarily the energy you would get by calling the EOS on the
  state. It is very rare (never?) that a RHS implementation would need
  to use this variable.

* ``state % T`` : the current temperature

* ``state % rho`` : the current density

Note that we come in with the mass fractions, but the molar fractions can
be computed as:

::

      double precision :: y(nspec)
      ...
      y(:) = state % xn(:) / aion(:)

The actual_rhs() routine’s job is to fill the righthand side vector
for the ODE system, state % ydot(:). Here, the important
fields to fill are:

* ``state % ydot(1:nspec_evolve)`` : the change in *molar
  fractions* for the ``nspec_evolve`` species that we are evolving,
  :math:`d({Y}_k)/dt`

* ``state % ydot(net_ienuc)`` : the change in the energy release
  from the net, :math:`d\enuc/dt`

* ``state % ydot(net_itemp)`` : the change in temperature, :math:`dT/dt`

The righthand side routine is assumed to return the change in *molar fractions*,
:math:`dY_k/dt`. These will be converted to the change in mass fractions, :math:`dX_k/dt`
by the wrappers that call the righthand side routine for the integrator.
If the network builds the RHS in terms of mass fractions directly, :math:`dX_k/dt`, then
these will need to be converted to molar fraction rates for storage, e.g.,
:math:`dY_k/dt = A_k^{-1} dX_k/dt`.

The Jacobian is provided by actual_jac(state), and takes the
form:

::

      subroutine actual_jac(state)
        type (burn_t) :: state

The Jacobian matrix elements are stored in ``state % jac`` as:

* ``state % jac(m, n)`` for :math:`\mathrm{m}, \mathrm{n} \in [1, \mathrm{nspec\_evolve}]` :
  :math:`d(\dot{Y}_m)/dY_n`

* ``state % jac(net_ienuc, n)`` for :math:`\mathrm{n} \in [1, \mathrm{nspec\_evolve}]` :
  :math:`d(\edot)/dY_n`

* ``state % jac(net_itemp, n)`` for :math:`\mathrm{n} \in [1, \mathrm{nspec\_evolve}]` :
  :math:`d(\dot{T})/dY_n`

* ``state % jac(m, net_itemp)`` for :math:`\mathrm{m} \in [1, \mathrm{nspec\_evolve}]` :
  :math:`d(\dot{Y}_m)/dT`

* ``state % jac(net_ienuc, net_itemp)`` :
  :math:`d(\edot)/dT`

* ``state % jac(net_itemp, net_itemp)`` :
  :math:`d(\dot{T})/dT`

* ``state % jac(p, net_ienuc)`` :math:`= 0` for :math:`\mathrm{p} \in [1, \mathrm{neqs}]`, since nothing
  should depend on the integrated energy release

The form looks like:

.. math::
   \left (
   \begin{matrix}
      \ddots  & \vdots                          &          & \vdots & \vdots \\
      \cdots  & \partial \dot{Y}_m/\partial Y_n & \cdots   & 0      & \partial \dot{Y}_m/\partial T    \\
              & \vdots                          & \ddots   & \vdots & \vdots  \\
      \cdots  & \partial \edot/\partial Y_n     & \cdots   & 0      & \partial \edot/\partial T   \\
      \cdots  & \partial \dot{T}/\partial Y_n   & \cdots   & 0      & \partial \dot{T}/\partial T   \\
   \end{matrix}
   \right )

This shows that all of the derivatives with respect to the nuclear
energy generated, :math:`e_\mathrm{nuc}` are zero. Again, this is because
this is just a diagnostic variable.

Note: a network is not required to compute a Jacobian if a numerical
Jacobian is used. This is set with the runtime parameter
``jacobian`` = 2, and implemented in
``integration/numerical_jacobian.f90`` using finite-differences.

Thermodynamics and :math:`T` Evolution
--------------------------------------

Burning Mode
^^^^^^^^^^^^

There are several different modes under which the burning can be done, set
via the burning_mode runtime parameter:

* ``burning_mode`` = 0 : hydrostatic burning

  :math:`\rho`, :math:`T` remain constant

* ``burning_mode = 1`` : self-heating burn

  :math:`T` evolves with the burning according to the temperature
  evolution equation. This is the “usual” way of thinking of the
  burning—all three equations (:eq:`eq:spec_integrate`,
  :eq:`eq:enuc_integrate`, and :eq:`eq:temp_integrate`) are solved
  simultaneously.

* ``burning_mode = 2`` : hybrid approach

  This implements an approach from :cite:`raskin:2010` in which we do
  a hydrostatic burn everywhere, but if we get a negative energy
  change, the burning is redone in self-heating mode (the logic being
  that a negative energy release corresponds to NSE conditions)

* ``burning_mode = 3`` : suppressed burning

  This does a self-heating burn, but limits all values of the RHS by a
  factor :math:`L = \text{min}(1, f_s (e / \dot{e}) / t_s)`, such that
  :math:`\dot{e} = f_s\, e / t_s`, where :math:`f_s` is a safety
  factor, set via burning_mode_factor.

When the integration is started, the burning mode is used to identify
whether temperature evolution should take place. This is used to
set the self_heat field in the burn_t type passed
into the RHS and Jacobian functions.

EOS Calls
^^^^^^^^^

The evolution of the thermodynamic quantities (like specific heats and
other partial derivatives) can be frozen during the integration to
their initial values, updated for every RHS call, or something
in-between. Just before we call the network-specific RHS routine, we
update the thermodynamics of our state (by calling
update_thermodynamics) [2]_ The thermodynamic quantity update depends on two runtime
parameters, call_eos_in_rhs and dT_crit:

* ``call_eos_in_rhs`` = T:

  We call the EOS just before every RHS evaluation, using :math:`\rho,
  T` as inputs. Therefore, the thermodynamic quantities will always be
  consistent with the input state.

* ``call_eos_in_rhs = F``

  Here we keep track of the temperature, :math:`T_\mathrm{old}`, at
  which the EOS was last called (which may be the start of
  integration).

  If

  .. math:: \frac{T - T_\mathrm{old}}{T} > \mathtt{dT\_crit}

  then we update the thermodynamics. We also compute :math:`d(c_v)/dT`
  and :math:`d(c_p)/dT` via differencing with the old thermodynamic
  state and store these in the integrator. If this inequality is not
  met, then we don’t change the thermodynamics, but simply update the
  composition terms in the EOS state, e.g., :math:`\bar{A}`.

  We interpret ``dT_crit`` as the fractional change needed in the
  temperature during a burn to trigger an EOS call that updates the
  thermodynamic variables. Note that this is fully independent of
  ``call_eos_in_rhs``.

:math:`T` Evolution
^^^^^^^^^^^^^^^^^^^

A network is free to write their own righthand side for the
temperature evolution equation in its ``actual_rhs()`` routine.
But since this equation only really needs to know the instantaneous
energy generation rate, :math:`\dot{e}`, most networks use the helper
function, ``temperature_rhs`` (in
``integration/temperature_integration.f90``):

::

      subroutine temperature_rhs(state)
        type (burn_t) :: state

This function assumes that ``state % ydot(net_ienuc)`` is already
filled and simply fills ``state % ydot(net_itemp)`` according to
the prescription below.

We need the specific heat, :math:`c_x`. Note that because we are evaluating
the temperature evolution independent of any hydrodynamics, we do not
incorporate any advective or :math:`pdV` terms in the evolution. Therefore,
for our approximation here, we need to decide which specific heat we
want—usually either the specific heat at constant pressure, :math:`c_p`,
or the specific heat at constant volume, :math:`c_v`. The EOS generally
will provide both of these specific heats; which one to use is a
choice the user needs to make based on the physics of their problem.
This is controlled by the parameter do_constant_volume_burn,
which will use :math:`c_v` if .true. and :math:`c_p` is .false.. See
:cite:`maestro:III` for a discussion of the temperature evolution
equation.

A fully accurate integration of Equation `[eq:temp_integrate] <#eq:temp_integrate>`__
requires an evaluation of the specific heat at each integration step,
which involves an EOS call given the current temperature. This is done
if ``call_eos_in_rhs`` = T, as discussed above.
This may add significantly to the expense of the calculation,
especially for smaller networks where construction of the RHS is
inexpensive

For ``call_eos_in_rhs`` = F, we can still capture some evolution
of the specific heat by periodically calling the EOS (using
``dT_crit`` as described above) and extrapolating to the current
temperature as:

.. math:: c_x = (c_x)_0 + \frac{T - T_0}{d(c_x)/dT|_0}

where the ‘:math:`_0`’ quantities are the values from when the EOS was last
called. This represents a middle ground between fully on and fully
off.

Note: if ``state % self_heat`` = F, then we do not evolve
temperature.

The runtime parameter integrate_temperature can be used to disable
temperature evolution (by zeroing out ``ydot(net_itemp)``).

Energy Integration
^^^^^^^^^^^^^^^^^^

The last equation in our system is the nuclear energy release,
:math:`\edot`. Because of the operator-split approach to this ODE system,
this is not the true specific internal energy, :math:`e` (since it only
responds only to the nuclear energy release and no pdV work).

At initialization, :math:`e` is set to the value from the EOS consistent
with the initial temperature, density, and composition:

.. math:: e_0 = e(\rho_0, T_0, {X_k}_0)

In the integration routines, this is termed the “energy offset”.

As the system is integrated, :math:`e` is updated to account for the
nuclear energy release,

.. math:: e(t) = e_0 + \int_{t_0}^t f(\dot{Y}_k) dt

Note that thermodynamic consistency will no longer be maintained
(because density doesn’t evolve and the :math:`T` evolution is approximate)
but :math:`e` will represent an approximation to the current specific
internal energy, including the nuclear energy generation release.

Upon exit, we subtract off this initial offset, so ``% e`` in
the returned ``burn_t`` type from the ``actual_integrator``
call represents the energy *release* during the burn.

.. _ch:networks:nspec_evolve:

nspec_evolve Implementation
---------------------------

For networks where ``nspec_evolve`` < ``nspec``, it may be
necessary to reset the species mass fractions each time you enter the
RHS routine. As an example, the network ignition_chamulak has
3 species, :math:`^{12}\mathrm{C}`, :math:`^{16}\mathrm{O}`, and :math:`^{24}\mathrm{Mg}`. In this
network, :math:`^{16}\mathrm{O}` is not evolved at all and any change in
:math:`^{12}\mathrm{C}` is reflected in :math:`^{24}\mathrm{Mg}`. So we can evolve only the
equation for :math:`^{12}\mathrm{C}`. The algebraic relation between the
unevolved mass fractions that must be enforced then is:

.. math:: X(\isotm{Mg}{24}) = 1 - X(\isotm{C}{12}) - X(\isotm{O}{16})

This is implemented in the routine ``update_unevolved_species``:

::

      subroutine update_unevolved_species(state)
        type (burn_t) :: state

This needs to be explicitly called in ``actual_rhs`` before
the mass fractions from the input state are accessed. It is
also called directly by the integrator at the end of integration,
to make sure the final state is consistent.



Renormalization
---------------

The ``renormalize_abundances`` parameter controls whether we
renormalize the abundances so that the mass fractions sum to one
during a burn. This has the positive benefit that in some cases it can
prevent the integrator from going off to infinity or otherwise go
crazy; a possible negative benefit is that it may slow down
convergence because it interferes with the integration
scheme. Regardless of whether you enable this, we will always ensure
that the mass fractions stay positive and larger than some floor
``small_x``.

Tolerances
==========

Tolerances dictate how accurate the ODE solver must be while solving
equations during a simulation.  Typically, the smaller the tolerance
is, the more accurate the results will be.  However, if the tolerance
is too small, the code may run for too long or the ODE solver will
never converge.  In these simulations, rtol values will set the
relative tolerances and atol values will set the absolute tolerances
for the ODE solver.  Often, one can find and set these values in an
input file for a simulation.

:numref:`fig:tolerances` shows the results of a simple simulation using the
burn_cell unit test to determine
what tolerances are ideal for simulations.
For this investigation, it was assumed that a run with a tolerance of :math:`10^{-12}`
corresponded to an exact result,
so it is used as the basis for the rest of the tests.
From the figure, one can infer that the :math:`10^{-3}` and :math:`10^{-6}` tolerances
do not yeild the most accurate results
because their relative error values are fairly large.
However, the test with a tolerance of :math:`10^{-9}` is accurate
and not so low that it takes incredible amounts of computer time,
so :math:`10^{-9}` should be used as the default tolerance in future simulations.

.. _fig:tolerances:
.. figure:: tolerances.png
   :alt: Relative error plot
   :width: 100%

   Relative error of runs with varying tolerances as compared
   to a run with an ODE tolerance of :math:`10^{-12}`.

.. _ch:networks:integrators:

Stiff ODE Solvers
=================

The integration tolerances for the burn are controlled by
``rtol_spec``, ``rtol_enuc``, and ``rtol_temp``,
which are the relative error tolerances for
:eq:`eq:spec_integrate`, :eq:`eq:enuc_integrate`, and
:eq:`eq:temp_integrate`, respectively. There are corresponding
``atol`` parameters for the absolute error tolerances. Note that
not all integrators handle error tolerances the same way—see the
sections below for integrator-specific information.

We integrate our system using a stiff ordinary differential equation
integration solver. The absolute error tolerances are set by default
to :math:`10^{-12}` for the species, and a relative tolerance of :math:`10^{-6}`
is used for the temperature and energy. The integration yields the
new values of the mass fractions, :math:`Y_k^\outp`.

There are several options for integrators. Each should be capable of
evolving any of the networks, but varying in their approach. Internally,
the integrators uses different data structures to store the integration
progress (from the old-style rpar array in VODE to derived
types), and each integrator needs to provide a routine to convert
from the integrator’s internal representation to the ``burn_t``
type required by the ``actual_rhs`` and ``actual_jac`` routine.

The name of the integrator can be selected at compile time using
the ``INTEGRATOR_DIR`` variable in the makefile. Presently,
the allowed options are BS, VBDF, VODE, and VODE90.

actual_integrator
-----------------

The entry point to the integrator is actual_integrator:

::

      subroutine actual_integrator(state_in, state_out, dt, time)

        type (burn_t), intent(in   ) :: state_in
        type (burn_t), intent(inout) :: state_out
        real(dp_t),    intent(in   ) :: dt, time

A basic flow chart of this interface is as follows (note: there are
many conversions between ``eos_t``, ``burn_t``, and any
integrator-specific type implied in these operations):

#. Call the EOS on the input state, using :math:`\rho, T` as the input
   variables.

   This involves:

   #. calling ``burn_to_eos`` to produce an ``eos_t``
      with the thermodynamic information.

   #. calling the EOS

   #. calling ``eos_to_XX``, where XX is, e.g.
      bs, the integrator type. This copies all of the relevant
      data into the internal representation used by the integrator.

   We use the EOS result to define the energy offset for :math:`e`
   integration.

#. Compute the initial :math:`d(c_x)/dT` derivatives, if necessary, by
   finite differencing on the temperature.

#. Call the main integration routine to advance the inputs state
   through the desired time interval, producing the new, output
   state.

#. If necessary (integration failure, burn_mode demands)
   do any retries of the integration

#. Subtract off the energy offset—we now store just the
   energy release as ``state_out % e``

#. Convert back to a ``burn_t`` type (by ``calling XX_to_burn``).

#. update any unevolved species, for ``nspec_evolve`` <
   ``nspec``

#. normalize the abundances so they sum to 1

Righthand side wrapper
----------------------

Each integrator does their own thing to construct the solution,
but they will all need to assess the RHS in ``actual_rhs``,
which means converting from their internal representation
to the ``burn_t`` type. This is handled in a file
called ``XX_rhs.f90``, where XX is the integrator name.
The basic outline of this routine is:

#. call ``clean_state``

   This function operates on the ODE integrator vector directly
   (accessing it from the integrator’s internal data structure). It
   makes sure that the mass fractions lie between :math:`[10^{-30}, 1]` and
   that the temperature lies between :math:`[{\tt small\_temp}, 10^{11}]`. The
   latter upper limit is arbitrary, but is safe for the types of problems
   we support with these networks.

#. renormalize the species, if ``renormalize_abundances`` = T

#. update the thermodynamic quantities if we are doing
   ``call_eos_in_rhs`` or the ``dT_crit`` requires

#. convert to a ``burn_t`` type and call the actual RHS

#. convert derives to mass-fraction-based (if
   ``integrate_molar_fraction`` = F and zero out the temperature or
   energy derivatives (if ``integrate_temperature`` = F or
   ``integrate_energy`` = F, respectively).

#. limit the rates if ``burning_mode`` = 3

#. convert back to the integrator’s internal representation

.. _sec:BS:

BS
--

This integrator is based on the stiff-ODE methods from :cite:`NR`, but
written with reaction network integration in mind (so it knows about
species), and in a modular / threadsafe fashion to work with our data
structures. This integrator appears quite robust.

bs_t data structure.
^^^^^^^^^^^^^^^^^^^^

The ``bs_t`` type is the main data structure for the BS integrator.
This holds the integration variables (as ``y(1:neqs)``) and data
associated with the timestepping. It also holds a ``burn_t`` type
as ``bs_t % burn_s``. This component is used to interface with
the networks. The conversion routines ``bs_to_burn`` and
``burn_to_bs`` simply sync up ``bs_t % y(:)`` and ``bs_t % burn_s``.

The ``upar(:)`` component contains the meta-data that is not held in
the ``burn_t`` but nevertheless is associate with the current
state. This is an array that can be indexed via the integers define
in the ``rpar_indices`` module. Note that because the ``bs_t``
contains its own ``burn_t`` type, the BS integrator does not need
as much meta-data as some other integrators. The fields of upar
are:

* ``bs_t % upar(irp_nspec : irp_nspec-1+n_not_evolved)``

  These are the mass fractions of the ``nspec`` - ``nspec_evolve``
  species that are not evolved in the ODE system.

* ``bs_t % upar(irp_y_init : irp_y_init-1+neqs)``

  This is the initial values of the ODE integration vector.

* ``bs_t % upar(irp_t_sound)``

  This is the sound-crossing time for a zone.

* ``bs_t % upar(irp_t0)``

  This is the simulation time at the start of integration. This can be
  used as an offset to convert between simulation time and integration
  time (we always start the integration at :math:`t = 0`).

Error criteria.
^^^^^^^^^^^^^^^

There is a single relative tolerance used for all ODEs, instead of a
separate one for species, temperature, and energy, it is simply the
maximum of {``rtol_spec``, ``rtol_temp``, ``rtol_enuc``}. The absolute
tolerance parameters are ignored.

A relative tolerance needs a metric against which to compare. BS
has two options, chosen by the runtime parameter scaling_method.
Considering a vector :math:`{\bf y} = (Y_k, e, T)^\intercal`, the scales
to compare against, :math:`{\bf y}_\mathrm{scal}`, are:

* ``scaling_method`` = 1 :

  .. math:: {\bf y}_\mathrm{scal} = |{\bf y}| + \Delta t  |\dot{\bf y}| + \epsilon

  This is an extrapolation of :math:`{\bf y}` in time. The quantity
  :math:`\epsilon` is a small number (hardcoded to :math:`10^{-30}`)
  to prevent any scale from being zero.

* ``scaling_method`` = 2 :

  .. math:: ({y}_\mathrm{scal})_j = \max \left (|y_j|, \mathtt{ode\_scale\_floor} \right )

  for :math:`j = 1, \ldots, {\tt neq}`. Here, ode_scale_floor is a
   runtime parameter that sets a lower-limit to the scaling for each
   variable in the vector :math:`{\bf y}_\mathrm{scal}`. The default
   value is currently :math:`10^{-6}` (although any network can
   override this using priorities). The effect of this scaling is that
   species with an abundance :math:`\ll` ``ode_scal_floor`` will not be
   used as strongly in assessing the accuracy of a step.

These correspond to the options presented in :cite:`NR`.

A final option, use_timestep_estimator enables the
timestep estimator from VODE to determine a good starting
timestep for integration.

.. _sec:VODE:

VODE
----

VODE is a classic integration package described in :cite:`vode`. We
use the implicit integration method in the package.

data structures.
^^^^^^^^^^^^^^^^

VODE does not allow for derived types for its internal representation
and instead simply uses a solution vector, ``y(neqs)``, and an array of
parameters, ``rpar(:)``. The indices into ``rpar`` are defined in the
``rpar_indices`` module.

tolerances.
^^^^^^^^^^^

Our copy of VODE is made threadsafe by use of the OpenMP
threadprivate directive on Fortran common blocks. However, due to
the use of computed go tos, we have not ported it to GPUs using
OpenACC.

.. _sec:VBDF:

VBDF
----

VBDF is a modern implementation of the methods in VODE, written by
Matt Emmett. It supports integrating a vector of states, but we do
not use that capability.

The error criteria is the same as VODE—separate relative, ``rtol``,
and absolute, ``atol``, error tolerances are specified for each
variable that is being integrated. A weight is constructed as:

.. math:: W_m = \frac{1}{{\tt rtol}_m |y_m| + {\tt atol}_m}

where needed, the error, :math:`\epsilon`, is constructed by computing an :math:`L_2`
norm:

.. math:: \epsilon = \left [ \frac{1}{N} \sum_m (y_m W_m)^2 \right ]^{1/2}

where :math:`m = 1, \ldots, N` indexes the ODE solution vector. With this
weighting, :math:`\epsilon < 1` means we’ve achieved our desired accuracy.

Retries
-------

Overriding Parameter Defaults on a Network-by-Network Basis
===========================================================

Any network can override or add to any of the existing runtime
parameters by creating a ``_parameters`` file in the network directory
(e.g., ``networks/triple_alpha_plus_cago/_parameters``). As noted in
Chapter [chapter:parameters], the fourth column in the ``_parameter``
file definition is the *priority*. When a duplicate parameter is
encountered by the scripts writing the ``extern_probin_module``, the value
of the parameter with the highest priority is used. So picking a large
integer value for the priority in a network’s ``_parameter`` file will
ensure that it takes precedence.

.. raw:: latex

   \centering

|image|

.. [1]
   Note that in previous versions of our networks in
   CASTRO and MAESTRO, there was another term in the temperature
   equation relating to the chemical potential of the gas as it came
   from the EOS. We have since decided that this term should
   analytically cancel to zero in all cases for our nuclear networks,
   and so we no longer think it is correct to include a numerical
   approximation of it in the integration scheme. So the current
   results given by our networks will in general be a little different
   than in the past.

.. [2]
   Note: each integrator provides its
   own implementation of this, since it operates on the internal
   data-structure of the integrator, but the basic procedure is the
   same.

.. |image| image:: doxygen_network.png


