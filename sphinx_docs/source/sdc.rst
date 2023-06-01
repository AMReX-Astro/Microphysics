*****************************
Spectral Deferred Corrections
*****************************

Introduction
============

The Simplified-SDC method provides a means to more strongly couple the
reactions to the hydrodynamics by evolving the reactions together with
an approximation of the advection over the timestep.  The full details
of the algorithm are presented in :cite:`castro_simple_sdc`.

We want to solve the coupled equations:

.. math:: \Uc_t = \Advs{\Uc} + \Rb(\Uc)

where :math:`\Uc` is the conserved state vector, :math:`\Uc = (\rho,
(\rho X_k), (\rho \Ub), (\rho E))^\intercal`.  Here, :math:`\rho` is
the density, :math:`X_k` are the mass fractions constrained via
:math:`\sum_k X_k = 1`, :math:`\Ub` is the velocity vector, and
:math:`E` is the specific total energy, related to the specific
internal energy, :math:`e`, via

.. math::

   E = e + |\Ub|^2/2 .

The quantity :math:`\Advs{\Uc}` represents the advective source term (including any
hydrodynamical sources),

.. math:: \Advs{\Uc} = - \nabla \cdot \Fb(\Uc) + \Shydro

and :math:`\Rb(\Uc)`
is the reaction source term.

Interface and Data Structures
=============================

burn_t
------

To accommodate the increased information required to evolve the
coupled system, the `burn_t` used for Strang integration is extended
to include the conserved state and the advective sources.  This is
used to pass information to/from the integration routine from the
hydrodynamics code.

ODE system
==========

The reactions don’t modify the total density, :math:`\rho`, or momentum,
:math:`\rho \Ub`, so our ODE system is just:

.. math::

   \frac{d}{dt}\left ( 
      \begin{array}{c} \rho X_k \\ \rho e \end{array} 
   \right ) = 
   \left ( \begin{array}{c}
      \Adv{\rho X_k}^{n+1/2} \\ \Adv{\rho e}^{n+1/2} \\
   \end{array} \right ) +
   \left (
      \begin{array}{c} \rho \omegadot_k \\ \rho \Sdot \end{array}
   \right )

Here the advective terms are piecewise-constant (in time)
approximations to the change in the state due to the hydrodynamics,
computed with the during the hydro step.

However, to define the temperature, we need the density at any
intermediate time, :math:`t`. We construct these as needed from the
time-advanced momenta:

.. math::

   \rho(t) = \rho^n + \Adv{\rho}^{n+1/2} (t - t^n)

Interfaces
==========

actual_integrator
-----------------

The main driver, ``actual_integrator``, is nearly identical to the Strang counterpart.  The
main difference is that it interprets the absolute tolerances in terms of :math:`(\rho X_k)`
instead of :math:`X_k`.

The flow of this main routine is:

#. Convert from the ``burn_t`` type to the integrator’s internal
   representation using ``burn_to_int()``.

   This copies the state variables into the
   integration type and stores the initial density.

#. Call the main integration routine to advance the inputs state
   through the desired time interval, producing the new, output state.

#. Convert back from the internal representation to the ``burn_t`` by
   calling ``int_to_burn()``.

Righthand side wrapper
----------------------

The manipulation of the righthand side is a little more complex
now.  Each network only provides the change in molar fractions and
internal energy, but we need to convert these to the conservative
system we are integrating, including the advective terms.

.. note::

   Presently only the ``VODE`` and ``BackwardEuler`` integrators supports SDC evolution.

#. Get the current density by calling ``update_density_in_time()``

#. Call ``clean_state`` to ensure that the mass fractions are valid

#. Convert the integrator-specific data structures to a ``burn_t`` via ``int_to_burn()``.

   a. Update the density (this may be redundant).

   b. Fill the ``burn_t``'s ``xn[]``, auxiliary data, internal energy,
      and call the EOS to get the temperature.

#. Call the ``actual_rhs()`` routine to get just the reaction sources
   to the update. In
   particular, this returns the change in molar fractions,
   :math:`\dot{Y}_k` and the nuclear energy release, :math:`\dot{S}`.

#. Convert back to the integrator’s internal representation via ``rhs_to_int``
   This converts the ``ydot``s to mass fractions and adds the advective terms
   to all ``ydots``.

Jacobian
--------

The Jacobian of this system is :math:`{\bf J} = \partial \Rb /
\partial \Uc`, since :math:`\Advs{\Uc}` is held constant during the
integration.  We follow the approach of :cite:`castro_simple_sdc` and factor
the Jacobian as:

.. math::

   {\bf J} = \frac{\partial \Rb}{\partial \Uc} = \frac{\partial \Rb}{\partial {\bf w}}
             \frac{\partial {\bf w}}{\partial \Uc}

where :math:`{\bf w} = (X_k, T)^\intercal` are the more natural variables
for a reaction network.

