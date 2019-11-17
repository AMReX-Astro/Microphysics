*****************************
Spectral Deferred Corrections
*****************************

Introduction
============

The Simplified-SDC method provides a means to more strongly couple the
reactions to the hydrodynamics by evolving the reactions together with
an approximation of the advection over the timestep.

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

sdc_t
-----

To accommodate the increased information required to evolve the
coupled system, we introduce a new data type, ``sdc_t``. This is
used to pass information to/from the integration routine from the
hydrodynamics code.

ODE system
==========

The reactions don’t modify the total density, :math:`\rho`, or momentum,
:math:`\rho \Ub`, so our ODE system is just:

.. math::

   \frac{d}{dt}\left ( 
      \begin{array}{c} \rho X_k \\ \rho E \\  \rho e \end{array} 
   \right ) = 
   \left ( \begin{array}{c}
      \Adv{\rho X_k}^{n+1/2} \\ \Adv{\rho E}^{n+1/2} \\ \Adv{\rho e}^{n+1/2} \\
   \end{array} \right ) +
   \left (
      \begin{array}{c} \rho \omegadot_k \\ \rho \Sdot \\ \rho \Sdot \end{array}
   \right )

where we include :math:`e` in addition to :math:`E` to provide
additional thermodynamic information to help find a consistent
:math:`T`. Here the advective terms are piecewise-constant (in time)
approximations to the change in the state due to the hydrodynamics,
computed with the during the hydro step.

However, to define the temperature, we need the kinetic energy (and
hence momentum and density) at any intermediate time, :math:`t`. We construct
these as needed from the time-advanced momenta:

.. math::

   (\rho \Ub)(t) = (\rho \Ub)^n + \Adv{\rho \Ub}^{n+1/2} (t - t^n)

Interfaces
==========

actual_integrator
-----------------

Simplified-SDC integration provides its own implementation of the main entry
point, ``actual_integrator``.

::

      subroutine actual_integrator(state_in, state_out, dt, time)

        type (sdc_t), intent(in   ) :: state_in
        type (sdc_t), intent(inout) :: state_out
        real(dp_t),    intent(in   ) :: dt, time

The main difference here is that the type that comes in and out of the
interface is no longer a ``burn_t``, but instead is an
``sdc_t``.

The flow of this main routine is simpler than the non-SDC version:

#. Convert from the ``sdc_t`` type to the integrator’s internal
   representation (e.g., ``sdc_to_bs`` converts from a ``bs_t``
   for the BS integrator).

   This copies the state variables and advective sources into the
   integration type. Since we only actually integrate :math:`(\rho X_k),
   (\rho e), (\rho E)`, the terms corresponding to density and momentum
   are carried in an auxillary array (indexed through the rpar
   mechanism).

#. Call the main integration routine to advance the inputs state
   through the desired time interval, producing the new, output state.

#. Convert back from the internal representation (e.g., a
   ``bs_t``) to the ``sdc_t`` type.

Righthand side wrapper
----------------------

The manipulation of the righthand side is considerably more complex
now. Each network only provides the change in molar fractions and
internal energy, but we need to convert these to the conservative
system we are integrating, including the advective terms.

.. note::

   In the example here, we'll assume we are using the VODE90
   integrator, but the VODE and BS integrators also support
   simplified-SDC integration.


#. Call ``clean_state``

#. Convert the integrator-specific data structures (``y`` and ``rpar``
   for VODE90) to a ``burn_t`` type in ``vode_to_burn``.

   a. call ``fill_unevolved_variables`` to update the density
      and momentum. Since these don’t depend on reactions, this is a
      simple algebraic update based on the piecewise-constant-in-time
      advective term:

      .. math::

         \begin{aligned}
               \rho(t) &= \rho^n + t \cdot \left [ \mathcal{A}(\rho) \right]^{n+1/2} \\
               (\rho \Ub)(t) &= (\rho \Ub)^n + t  \cdot \left [ \mathcal{A}(\rho\Ub) \right]^{n+1/2} 
             \end{aligned}

      where here we assume that we are integrating in the ODE system
      starting with :math:`t=0`.

   b. compute the specific internal energy, :math:`e`, from either the
      total minus kinetic energy or internal energy carried by the
      integrator (depending on the value of ``T_from_eden``).

   c. get the temperature from the equation of state

   d. convert to a ``burn_t`` type via ``eos_to_burn``.  This
      ``burn_t`` variable has fields for ``rho``, ``T``, and
      ``xn``—everything that is needed to evaluate the reaction rates.

#. Call the ``network_rhs`` routine to get just the reaction sources
   to the update (this simply calls the network's ``actual_rhs``). In
   particular, this returns the change in molar fractions,
   :math:`\dot{Y}_k` and the nuclear energy release, :math:`\dot{S}`.

#. Convert back to the integrator’s internal representation (e.g.,
   a ``y``, ``ydot``, and ``rpar`` for VODE via ``rhs_to_vode``)

    a. call ``fill_unevolved_variables``

    b. fill the ydot array in the integrator type (e.g.,
       ``rpar``) with the advective sources that originally came into the
       intergrator through the ``sdc_t``.

    c. Add the reacting terms. This is done as:

      .. math::

         \begin{aligned}
               \dot{y}_{\rho X_k} &= \Adv{\rho X_k}^{n+1/2} + \rho A_k \dot{Y}_k \\
               \dot{y}_{\rho e} &= \Adv{\rho e}^{n+1/2} +\rho \dot{S} \\
               \dot{y}_{\rho E} &= \Adv{\rho E}^{n+1/2} + \rho \dot{S}
             \end{aligned}

Jacobian
--------

The Jacobian of this system is :math:`{\bf J} = \partial \Rb /
\partial \Uc`, since :math:`\Advs{\Uc}` is held constant during the
integration.  We follow the approach of :cite:`castro_sdc` and factor
the Jacobian as:

.. math::

   {\bf J} = \frac{\partial \Rb}{\partial \Uc} = \frac{\partial \Rb}{\partial {\bf w}}
             \frac{\partial {\bf w}}{\partial \Uc}

where :math:`{\bf w} = (\rho, X_k, T)^\intercal` are the more natural variables
for a reaction network.

