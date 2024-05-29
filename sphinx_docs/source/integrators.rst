*******************
Reaction ODE System
*******************

.. important::

   This describes the integration done when doing Strang operator-splitting, which is the
   default mode of coupling burning to application codes.

The equations we integrate to do a nuclear burn are:

.. math::
   \frac{dX_k}{dt} = \dot{\omega}_k(\rho,X_k,T)
   :label: eq:spec_integrate

.. math::
   \frac{de}{dt} = f(\rho,X_k,T)
   :label: eq:enuc_integrate

Here, :math:`X_k` is the mass fraction of species :math:`k`, :math:`e` is the specific
nuclear energy created through reactions. Also needed are density :math:`\rho`,
temperature :math:`T`, and the specific heat. The function :math:`f` provides the energy release from reactions and can often be expressed in terms of the
instantaneous reaction terms, :math:`\dot{X}_k`. As noted in the previous
section, this is implemented in a network-specific manner.

In this system, :math:`e` is equal to the total specific internal
energy. This allows us to easily call the EOS during the burn to obtain the temperature.

.. note::

   The energy generation rate includes a term for neutrino losses in addition
   to the energy release from the changing binding energy of the
   fusion products.

.. index:: integrator.use_number_densities

.. note::

   By setting ``integrator.use_number_densities=1``, number densities will be
   integrated instead of mass fractions.  This makes the system:

   .. math::
      \frac{dn_k}{dt} = \dot{\omega}_k(\rho,n_k,T)
      :label: eq:spec_n_integrate

   .. math::
      \frac{de}{dt} = f(\rho,n_k,T)
      :label: eq:enuc_n_integrate

   The effect of this flag in the integrators is that we don't worry
   about converting between mass and molar fractions when calling the
   righthand side function and Jacobian, and we don't do any normalization
   requiring $\sum_k X_k = 1$.


While this is the most common way to construct the set of
burn equations, and is used in most of our production networks,
all of them are ultimately implemented by the network itself, which
can choose to disable the evolution of any of these equations by
setting the RHS to zero. The integration software provides some
helper routines that construct common RHS evaluations, like the routine
that converts a temperature update to :math:`\dot{e}`, but these calls
are always explicitly done by the individual networks rather than
being handled by the integration backend. This allows you to write a
new network that defines the RHS in whatever way you like.

.. index:: integrator.react_boost

The standard reaction rates can all be boosted by a constant factor by
setting the ``integrator.react_boost`` runtime parameter.  This will simply
multiply the righthand sides of each species evolution equation (and
appropriate Jacobian terms) by the specified constant amount.

Interfaces
==========

The interfaces to all of the networks and integrators are written in C++.

``burner``
----------

The main entry point for C++ is ``burner()`` in
``interfaces/burner.H``.  This simply calls the ``integrator()``
routine (at the moment this can be ``VODE``, ``BackwardEuler``, ``ForwardEuler``, ``QSS``, or ``RKC``).

.. code-block:: c++

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void burner (burn_t& state, Real dt)

The input is a ``burn_t``.

.. note::

   For the thermodynamic state, only the density, temperature, and
   mass fractions are used directly--we compute the internal energy
   corresponding to this input state through the equation of state
   before integrating.

When integrating the system, we often need auxiliary information to
close the system.  This is kept in the original ``burn_t`` that was
passed into the integration routines.  For this reason, we often need
to pass both the specific integrator's type (e.g. ``dvode_t``) and
``burn_t`` objects into the lower-level network routines.

The overall flow of the integrator is (using VODE as the example):

#. Call the EOS directly on the input ``burn_t`` state using :math:`\rho` and :math:`T` as inputs.

#. Fill the integrator type by calling ``burn_to_integrator()`` to create a
   ``dvode_t``.

#. call the ODE integrator, ``dvode()``, passing in the ``dvode_t`` _and_ the
   ``burn_t`` --- as noted above, the auxiliary information that is
   not part of the integration state will be obtained from the
   ``burn_t``.

#. subtract off the energy offset---we now store just the energy released
   in the ``dvode_t`` integration state.

#. convert back to a ``burn_t`` by calling ``integrator_to_burn``

#. normalize the abundances so they sum to 1.

.. index:: integrator.subtract_internal_energy

.. note::

   Upon exit, ``burn_t burn_state.e`` is the energy *released* during
   the burn, and not the actual internal energy of the state.

   Optionally, by setting ``integrator.subtract_internal_energy=0``
   the output will be the total internal energy, including that released
   burning the burn.

Network Routines
----------------

.. important::

   Microphysics integrates the reaction system in terms of mass
   fractions, :math:`X_k`, but most astrophysical networks use molar
   fractions, :math:`Y_k`.  As a result, we expect the networks to
   return the righthand side and Jacobians in terms of molar
   fractions.  The integration wrappers will internally
   convert to mass fractions as needed for the integrators.

Righthand size implementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The righthand side of the network is implemented by
``actual_rhs()`` in ``actual_rhs.H``, and appears as

.. code-block:: c++

   void actual_rhs(burn_t& state, Array1D<Real, 1, neqs>& ydot)

All of the necessary integration data comes in through state, as:

* ``state.xn[NumSpec]`` : the mass fractions.

* ``state.aux[NumAux]`` : the auxiliary data (only available if ``NAUX_NET`` > 0)

* ``state.e`` : the current internal energy. It is very rare (never?) that a RHS
  implementation would need to use this variable directly -- even though this is
  the main thermodynamic integration variable, we obtain the temperature from the
  energy through an EOS evaluation.

* ``state.T`` : the current temperature

* ``state.rho`` : the current density

Note that we come in with the mass fractions, but the molar fractions can
be computed as:

.. code-block:: c++

      Array1D<Real, 1, NumSpec> y;
      ...
      for (int i = 1; i <= NumSpec; ++i) {
          y(i) = state.xn[i-1] * aion_inv[i-1];
      }

.. warning::

   We use 1-based indexing for ``ydot`` for legacy reasons, so watch out when filling in
   this array based on 0-indexed C arrays.

The ``actual_rhs()`` routine’s job is to fill the righthand side vector
for the ODE system, ``ydot(neqs)``. Here, the important
fields to fill are:

* ``state.ydot(1:NumSpec)`` : the change in *molar
  fractions* for the ``NumSpec`` species that we are evolving,
  :math:`d({Y}_k)/dt`

* ``state.ydot(net_ienuc)`` : the change in the internal energy
  from the net, :math:`de/dt`

The righthand side routine is assumed to return the change in *molar fractions*,
:math:`dY_k/dt`. These will be converted to the change in mass fractions, :math:`dX_k/dt`
by the wrappers that call the righthand side routine for the integrator.
If the network builds the RHS in terms of mass fractions directly, :math:`dX_k/dt`, then
these will need to be converted to molar fraction rates for storage, e.g.,
:math:`dY_k/dt = A_k^{-1} dX_k/dt`.

Righthand side wrapper
^^^^^^^^^^^^^^^^^^^^^^

The integrator provides a wrapper that sits between the integration
routines and the network's implementation of the righthand side.  Its
flow is (for VODE):

#. call ``clean_state`` on the ``dvode_t``

#. update the thermodynamics by calling ``update_thermodynamics``.  This takes both
   the ``dvode_t`` and the ``burn_t`` and computes the temperature that matches the
   current state.

#. call ``actual_rhs``

#. convert the derivatives to mass-fraction-based (since we integrate :math:`X`)
   and zero out the temperature and energy derivatives if we are not integrating
   those quantities.

#. apply any boosting if ``react_boost`` > 0


Jacobian implementation
^^^^^^^^^^^^^^^^^^^^^^^

.. index:: integrator.jacobian

Either an analytic or numerical Jacobian is used for the implicit
integrators, selected via the ``integrator.jacobian`` runtime
parameter (``1`` = analytic; ``2`` = numerical).  For VODE, the
numerical Jacobian is computed internally.  For the other integrators,
a difference method is implemented in
``integration/utils/numerical_jacobian.H``.

The analytic Jacobian is specific to each network and is provided by
``actual_jac(state, jac)``.  It takes the form:

.. code-block:: c++

   void actual_jac(burn_t& state, MathArray2D<1, neqs, 1, neqs>& jac)

The Jacobian matrix elements are stored in ``jac`` as:

* ``jac(m, n)`` for :math:`\mathrm{m}, \mathrm{n} \in [1, \mathrm{NumSpec}]` :
  :math:`d(\dot{Y}_m)/dY_n`

* ``jac(net_ienuc, n)`` for :math:`\mathrm{n} \in [1, \mathrm{NumSpec}]` :
  :math:`d(\dot{e})/dY_n`

* ``jac(m, net_ienuc)`` for :math:`\mathrm{m} \in [1, \mathrm{NumSpec}]` :
  :math:`d(\dot{Y}_m)/de`

* ``jac(net_ienuc, net_ienuc)`` :
  :math:`d(\dot{e})/de`

The form looks like:

.. math::
   \left (
   \begin{matrix}
      \ddots  & \vdots                          &          & \vdots \\
      \cdots  & \partial \dot{Y}_m/\partial Y_n & \cdots   & \partial \dot{Y}_m/\partial e    \\
              & \vdots                          & \ddots   & \vdots  \\
      \cdots  & \partial \dot{e}/\partial Y_n   & \cdots   & \partial \dot{e}/\partial e   \\
   \end{matrix}
   \right )

.. note::

   A network is not required to provide a Jacobian if a numerical
   Jacobian is used.


Jacobian wrapper
^^^^^^^^^^^^^^^^

The integrator provides a wrapper that sits between the integration
routines and the network's implementation of the Jacobian.  Its
flow is (for VODE):

.. note::

   It is assumed that the thermodynamics are already correct when
   calling the Jacobian wrapper, likely because we just called the RHS
   wrapper above which did the ``clean_state`` and
   ``update_thermodynamics`` calls.

.. index:: integrator.react_boost

#. call ``integrator_to_burn()`` to update the ``burn_t``

#. call ``actual_jac()`` to have the network fill the Jacobian array

#. convert the derivative to be mass-fraction-based

#. apply any boosting to the rates if ``integrator.react_boost`` > 0





Thermodynamics and :math:`e` Evolution
======================================

The thermodynamic equation in our system is the evolution of the internal energy,
:math:`e`.

.. note::

   When the system is integrated in an operator-split approach, the
   energy equation accounts for only the nuclear energy release and
   not pdV work.

At initialization, :math:`e` is set to the value from the EOS consistent
with the initial temperature, density, and composition:

.. math::

   e_0 = e(\rho_0, T_0, {X_k}_0)

In the integration routines, this is termed the *energy offset*.

As the system is integrated, :math:`e` is updated to account for the
nuclear energy release,

.. math:: e(t) = e_0 + \int_{t_0}^t f(\dot{Y}_k) dt

As noted above, upon exit, we subtract off this initial offset, so ``state.e`` in
the returned ``burn_t`` type from the ``actual_integrator``
call represents the energy *release* during the burn.

Integration of Equation :eq:`eq:enuc_integrate`
requires an evaluation of the temperature at each integration step
(since the RHS for the species is given in terms of :math:`T`, not :math:`e`).
This involves an EOS call and is the default behavior of the integration.
Note also that for the Jacobian, we need the specific heat, :math:`c_v`, since we
usually calculate derivatives with respect to temperature (as this is the form
the rates are commonly provided in).

.. index:: integrator.call_eos_in_rhs

.. note::

   If desired, the EOS call can be skipped and the temperature and $c_v$ kept
   frozen over the entire time interval of the integration by setting ``integrator.call_eos_in_rhs=0``.
