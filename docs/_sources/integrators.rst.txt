*********************
Integrating a Network
*********************

Reaction ODE System
===================

The equations we integrate to do a nuclear burn are:

.. math::
   \frac{dX_k}{dt} = \omegadot_k(\rho,X_k,T)
   :label: eq:spec_integrate

.. math::
   \frac{d\enuc}{dt} = f(\dot{X}_k)
   :label: eq:enuc_integrate

.. math::
   \frac{dT}{dt} =\frac{\edot}{c_x} .
   :label: eq:temp_integrate

Here, :math:`X_k` is the mass fraction of species :math:`k`, :math:`\enuc` is the specifc
nuclear energy created through reactions, :math:`T` is the
temperature [1]_ , and :math:`c_x` is the specific heat for the
fluid. The function :math:`f` provides the energy release based on the
instantaneous reaction terms, :math:`\dot{X}_k`. As noted in the previous
section, this is implemented in a network-specific manner.

In this system, :math:`\enuc` is not necessarily the total specific internal
energy, but rather captures the energy release during the burn:

.. math:: \enuc = \int \edot dt

so we know how much energy was released (or required) over the burn.

.. note::

   This is the energy release that is used in computing the final
   energy generation rate for the burning.  By integrating a separate
   equation for :math:`\enuc`, we can account for neutrino losses as
   well as the energy release from the changing binding energy of the
   fusion products.

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

.. index:: react_boost

The standard reaction rates can all be boosted by a constant factor by
setting the ``react_boost`` runtime parameter.  This will simply
multiply the righthand sides of each species evolution equation (and
appropriate Jacobian terms) by the specified constant amount.

Interfaces
==========

There are both C++ and Fortran interfaces to all of the networks and
integrators.  For the most part, they implement similar data
structures to describe the evolution.

.. note::

   StarKiller integrates the reaction system in terms of mass fractions,
   :math:`X_k`, but most astrophysical networks use molar fractions,
   :math:`Y_k`.  As a result, we expect the networks to return the
   righthand side and Jacobians in terms of molar fractions.  The StarKiller
   routines will internally convert to mass fractions as needed for the
   integrators.


Fortran
-------


The righthand side of the network is implemented by
``actual_rhs()`` in ``actual_rhs.f90``, and appears as

.. code-block:: fortran

      subroutine actual_rhs(state)
        type (burn_t) :: state

All of the necessary integration data comes in through state, as:

* ``state % xn(:)`` : the ``nspec`` mass fractions.

* ``state % aux(:)`` : the ``naux`` auxillary data (only available if ``NAUX_NET`` > 0)

* ``state % e`` : the current value of the ODE system’s energy
  release, :math:`\enuc`—note: as discussed above, this is not
  necessarily the energy you would get by calling the EOS on the
  state. It is very rare (never?) that a RHS implementation would need
  to use this variable.

* ``state % T`` : the current temperature

* ``state % rho`` : the current density

Note that we come in with the mass fractions, but the molar fractions can
be computed as:

.. code-block:: fortran

      double precision :: y(nspec)
      ...
      y(:) = state % xn(:) * aion_inv(:)

The ``actual_rhs()`` routine’s job is to fill the righthand side vector
for the ODE system, ``state % ydot(:)``. Here, the important
fields to fill are:

* ``state % ydot(1:nspec)`` : the change in *molar
  fractions* for the ``nspec`` species that we are evolving,
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

.. code-block:: fortran

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

C++
---

.. note::

   Currently, only the VODE solver supports C++, so the interfaces
   here are specific to that integrator.

The righthand side implementation of the network has the interface:

.. code-block:: c++

   AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
   void actual_rhs(burn_t& state, Array1D<Real, 1, neqs>& ydot)

.. note::

   In the C++ implementation of the integrator, we use 1-based
   indexing for ``ydot`` to allow for easier conversion between
   Fortran and C++ networks.

The Jacobian has the form:

.. code-block:: c++

   template<class MatrixType>
   AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
   void actual_jac(burn_t& state, MatrixType& jac)

Here, ``MatrixType`` is either a ``SparseMatrix`` type or a
``RArray2D`` type (set in the ``dvode_t`` type in ``vode_type.H``.
This allows a network to use either sparse or dense linear algebra.
In the case of dense linear algebra, ``RArray2D`` is essentially a 2-d
array indexed from ``1`` to ``VODE_NEQS`` in each dimension.


Thermodynamics and :math:`T` Evolution
======================================


EOS Calls
---------

The evolution of the thermodynamic quantities (like specific heats and
other partial derivatives) can be frozen during the integration to
their initial values, updated for every RHS call, or something
in-between. Just before we call the network-specific RHS routine, we
update the thermodynamics of our state (by calling
``update_thermodynamics``) [2]_ The thermodynamic quantity update depends on two runtime
parameters, call_eos_in_rhs and dT_crit:

* ``call_eos_in_rhs = T``:

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
-------------------

A network is free to write their own righthand side for the
temperature evolution equation in its ``actual_rhs()`` routine.
But since this equation only really needs to know the instantaneous
energy generation rate, :math:`\dot{e}`, most networks use the helper
function, ``temperature_rhs``.  The Fortran implementation is in
``integration/utils/temperature_integration.f90``:

.. code-block:: fortran

      subroutine temperature_rhs(state)
        type (burn_t) :: state

This function assumes that ``state % ydot(net_ienuc)`` is already
filled and simply fills ``state % ydot(net_itemp)`` according to
the prescription below.

The C++ implementation is in ``integration/utils/temperature_integration.H``:

.. code-block:: c++

     AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
     void temperature_rhs (burn_t& state, Array1D<Real, 1, neqs>& ydot)

We need the specific heat, :math:`c_x`. Note that because we are evaluating
the temperature evolution independent of any hydrodynamics, we do not
incorporate any advective or :math:`pdV` terms in the evolution. Therefore,
for our approximation here, we need to decide which specific heat we
want—usually either the specific heat at constant pressure, :math:`c_p`,
or the specific heat at constant volume, :math:`c_v`. The EOS generally
will provide both of these specific heats; which one to use is a
choice the user needs to make based on the physics of their problem.
This is controlled by the parameter ``do_constant_volume_burn``,
which will use :math:`c_v` if ``.true.`` and :math:`c_p` is ``.false.``. See
:cite:`maestro:III` for a discussion of the temperature evolution
equation.

A fully accurate integration of Equation :eq:`eq:temp_integrate`
requires an evaluation of the specific heat at each integration step,
which involves an EOS call given the current temperature. This is done
if ``call_eos_in_rhs = T``, as discussed above.
This may add significantly to the expense of the calculation,
especially for smaller networks where construction of the RHS is
inexpensive

For ``call_eos_in_rhs = F``, we can still capture some evolution
of the specific heat by periodically calling the EOS (using
``dT_crit`` as described above) and extrapolating to the current
temperature as:

.. math:: c_x = (c_x)_0 + \frac{T - T_0}{d(c_x)/dT|_0}

where the ‘:math:`_0`’ quantities are the values from when the EOS was last
called. This represents a middle ground between fully on and fully
off.

Note: if ``state % self_heat = F`` (Fortran) or ``state.self_heat =
false`` (C++), then we do not evolve temperature.

The runtime parameter ``integrate_temperature`` can be used to disable
temperature evolution (by zeroing out ``ydot(net_itemp)``).

Energy Integration
==================

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



Renormalization
===============

The ``renormalize_abundances`` parameter controls whether we
renormalize the abundances so that the mass fractions sum to one
during a burn. This has the positive benefit that in some cases it can
prevent the integrator from going off to infinity or otherwise go
crazy; a possible negative benefit is that it may slow down
convergence because it interferes with the integration
scheme. Regardless of whether you enable this, we will always ensure
that the mass fractions stay positive and larger than some floor
``small_x``.


.. _ch:networks:integrators:

Stiff ODE Solvers
=================

We use high-order implicit ODE solvers for integrating the reaction
system.  There are several options for integrators. Each should be capable of
evolving any of the networks, but varying in their approach. Internally,
the integrators uses different data structures to store the integration
progress (from the old-style rpar array in VODE to derived
types), and each integrator needs to provide a routine to convert
from the integrator’s internal representation to the ``burn_t``
type required by the ``actual_rhs`` and ``actual_jac`` routine.

The name of the integrator can be selected at compile time using
the ``INTEGRATOR_DIR`` variable in the makefile. Presently,
the allowed options are:

* ``ForwardEuler``: an explicit first-order forward-Euler method.  This is
  meant for testing purposes only.

* ``VODE``: the VODE (:cite:`vode`) integration package.  We ported this
  integrator to C++ and removed the non-stiff integration code paths.

We recommend that you use the VODE solver, as it is the most
robust and has both Fortran and C++ implementations.

.. note::

   In the implementation details shown below, we write the flow in
   terms of the VODE solver routine names.

Tolerances
----------

Tolerances dictate how accurate the ODE solver must be while solving
equations during a simulation.  Typically, the smaller the tolerance
is, the more accurate the results will be.  However, if the tolerance
is too small, the code may run for too long or the ODE solver will
never converge.  In these simulations, ``rtol`` values will set the
relative tolerances and ``atol`` values will set the absolute tolerances
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

The integration tolerances for the burn are controlled by
``rtol_spec``, ``rtol_enuc``, and ``rtol_temp``,
which are the relative error tolerances for
:eq:`eq:spec_integrate`, :eq:`eq:enuc_integrate`, and
:eq:`eq:temp_integrate`, respectively. There are corresponding
``atol`` parameters for the absolute error tolerances. Note that
not all integrators handle error tolerances the same way—see the
sections below for integrator-specific information.

The absolute error tolerances are set by default
to :math:`10^{-12}` for the species, and a relative tolerance of :math:`10^{-6}`
is used for the temperature and energy.


Fortran interfaces
------------------

``integrator``
^^^^^^^^^^^^^^

The entry point to the integrator is ``integrator()`` in
``integration/integrator.F90``.  This does some setup and then calls
the specific integration routine, e.g., ``vode_integrator()`` in
``integration/VODE/vode_integrator.F90``.

.. code-block:: fortran

      subroutine vode_integrator(state_in, state_out, dt, time, status)

        type (burn_t), intent(in   ) :: state_in
        type (burn_t), intent(inout) :: state_out
        real(rt),    intent(in   ) :: dt, time
        type (integration_status_t), intent(inout) :: status

A basic flow chart of this interface is as follows (note: there are
many conversions between ``eos_t``, ``burn_t``, and any
integrator-specific type implied in these operations):

#. Call the EOS on the input state, using :math:`\rho, T` as the input
   variables.

   This involves:

   #. calling ``burn_to_eos`` to produce an ``eos_t``
      with the thermodynamic information.

   #. calling the EOS

   #. calling ``eos_to_vode`` to produce a ``dvode_t`` type
      containing all of the relevant
      data into the internal representation used by the integrator.
      Data that is not part of the integration state is stored in an ``rpar``
      array that is indexed using the integer keys in ``vode_rpar_indices``.

   We use the EOS result to define the energy offset for :math:`e`
   integration.

#. Compute the initial :math:`d(c_x)/dT` derivatives, if necessary, by
   finite differencing on the temperature.

#. Call the main integration routine, ``dvode()``, passing in the
   ``dvode_t`` state to advance the inputs state through the desired
   time interval, producing the new, output state.

#. If necessary (integration failure, burn_mode demands)
   do any retries of the integration

#. Subtract off the energy offset—we now store just the
   energy release as ``state_out % e``

#. Convert back to a ``burn_t`` type (by ``calling vode_to_burn``).

#. normalize the abundances so they sum to 1

Righthand side wrapper
^^^^^^^^^^^^^^^^^^^^^^

Each integrator does their own thing to construct the solution,
but they will all need to assess the RHS in ``actual_rhs``,
which means converting from their internal representation
to the ``burn_t`` type. This is handled in a file
called ``vode_rhs.F90``.
The basic outline of this routine is:

#. call ``clean_state``

   This function operates on the ODE integrator vector directly
   (accessing it from the integrator’s internal data structure). It
   makes sure that the mass fractions lie between ``SMALL_X_SAFE`` and 1  and
   that the temperature lies between :math:`[{\tt small\_temp}, {\tt MAX_TEMP}]`. The
   latter upper limit is arbitrary, but is safe for the types of problems
   we support with these networks.

   It also renormalizes the species, if ``renormalize_abundances = T``

#. update the thermodynamic quantities by calling
   ``update_thermodynamics()``

   among other things, this will handle the ``call_eos_in_rhs`` option
   or if the ``dT_crit`` requires the EOS call.

#. call ``vode_to_burn`` to convert to a ``burn_t``

#. call the actual RHS

#. convert derivatives to mass-fraction-based (since we integrate :math:`X`),
   and zero out the temperature or
   energy derivatives (if ``integrate_temperature = F`` or
   ``integrate_energy = F``, respectively).

#. apply any boosting to the rates if ``react_boost`` > 0.

#. convert back to the integrator’s internal representation by calling ``burn_to_vode``


Jacobian wrapper
^^^^^^^^^^^^^^^^

Similar to the RHS, the Jacobian wrapper is handled in the same
``vode_rhs.f90``.
The basic outline of this routine is:

.. note::

   It is assumed that the thermodynamics are already correct when
   calling the Jacobian wrapper, likely because we just called the RHS
   wrapper above which did the ``clean_state`` and
   ``update_thermodynamics`` calls.

#. call ``vode_to_burn`` to convert to a ``burn_t`` type.

#. call the actual Jacobian routine

#. convert derivatives to mass-fraction-based (since we integrate :math:`X`),
   and zero out the temperature or
   energy derivatives (if ``integrate_temperature = F`` or
   ``integrate_energy = F``, respectively).

#. apply any boosting to the rates if ``react_boost`` > 0.

#. convert back to the integrator’s internal representation by calling ``burn_to_vode``



C++ interfaces
--------------

``burner``
^^^^^^^^^^

The main entry point for C++ is ``burner()`` in
``interfaces/burner.H``.  This simply calls the ``integrator()``
routine, which at the moment is only provided by VODE.

.. code-block:: c++

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void burner (burn_t& state, Real dt)


The basic flow of the ``integrator()`` routine mirrors the Fortran one.

.. note::

   The C++ VODE integrator does not use a separate ``rpar`` array as
   part of the ``dvode_t`` type.  Instead, any auxillary information
   is kept in the original ``burn_t`` that was passed into the
   integration routines.  For this reason, we often need to pass both
   the ``dvode_t`` and ``burn_t`` objects into the network routines.

#. Call the EOS on the input ``burn_t`` state.  This involves:

   #. calling ``burn_to_eos`` to convert the ``burn_t`` to an ``eos_t``

   #. calling the EOS with :math:`\rho` and :math:`T` as input

   #. calling ``eos_to_burn`` to convert the ``eos_t`` back to a ``burn_t``

#. Fill the integrator type by calling ``burn_to_vode`` to create a
   ``dvode_t`` from the ``burn_t``

   .. note::

      unlike the Fortran interface, there is no ``vode_to_eos`` routine in C++

#. Compute the initial :math:`d(c_x)/dt` derivatives

#. call the ODE integrator, ``dvode()``, passing in the ``dvode_t`` _and_ the
   ``burn_t`` --- as noted above, the auxillary information that is
   not part of the integration state will be obtained from the
   ``burn_t``.

#. subtract off the energy offset---we now store just the energy released
   in the ``dvode_t`` integration state.

#. convert back to a ``burn_t`` by calling ``vode_to_burn``

#. normalize the abundances so they sum to 1.


Righthand side wrapper
^^^^^^^^^^^^^^^^^^^^^^

#. call ``clean_state`` on the ``dvode_t``

#. update the thermodynamics by calling ``update_thermodynamics``.  This takes both
   the ``dvode_t`` and the ``burn_t``.

#. call ``vode_to_burn`` to update the ``burn_t``

#. call ``actual_rhs``

#. convert the derivatives to mass-fraction-based (since we integrate :math:`X`)
   and zero out the temperature and energy derivatives if we are not integrating
   those quantities.

#. apply any boosting if ``react_boost`` > 0

#. convert back to the ``dvode_t`` by calling ``burn_to_vode``


Jacobian wrapper
^^^^^^^^^^^^^^^^

.. note::

   It is assumed that the thermodynamics are already correct when
   calling the Jacobian wrapper, likely because we just called the RHS
   wrapper above which did the ``clean_state`` and
   ``update_thermodynamics`` calls.

#. call ``vode_to_burn`` to update the ``burn_t``

#. call ``actual_jac()`` to have the network fill the Jacobian array

#. convert the derivative to be mass-fraction-based

#. apply any boosting to the rates if ``react_boost`` > 0

#. call ``burn_to_vode`` to update the ``dvode_t`` 




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
