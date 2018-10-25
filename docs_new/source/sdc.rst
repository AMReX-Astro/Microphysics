***
SDC
***

Introduction
============

SDC provides a means to more strongly couple the reactions to the
hydrodynamics by evolving the reactions together with an approximation
of the advection over the timestep.

We want to solve the coupled equations:

.. math:: {\bm{\mathcal{U}}}_t = {\mathcal{A} \left({\bm{\mathcal{U}}}\right)} + {\bf R}({\bm{\mathcal{U}}})

where :math:`{\bm{\mathcal{U}}}` is the conserved state vector, :math:`{\bm{\mathcal{U}}}= (\rho, (\rho X_k),
(\rho {\bf U}), (\rho E))^\intercal`, :math:`X_k` are the mass fractions
constrained via :math:`\sum_k X_k = 1`, :math:`{\bf U}` is the velocity vector, and
:math:`E` is the specific total energy, related to the specific internal
energy, :math:`e`, via :math:`E = e + |{\bf U}|^2/2`. Here :math:`{\mathcal{A} \left({\bm{\mathcal{U}}}\right)}` is the
advective source term (including any hydrodynamical sources),

.. math:: {\mathcal{A} \left({\bm{\mathcal{U}}}\right)} = - \nabla \cdot \mathbf{F}({\bm{\mathcal{U}}}) + \Sb_\mathrm{hydro}

and :math:`{\bf R}({\bm{\mathcal{U}}})`
is the reaction source term.

Interface and Data Structures
=============================

sdc_t
-----

To accommodate the increased information required to evolve the
coupled system, we introduce a new data type, sdc_t. This is
used to pass information to/from the integration routine from the
hydrodynamics code.

ODE system
==========

The reactions don’t modify the total density, :math:`\rho`, or momentum,
:math:`\rho {\bf U}`, so our ODE system is just:

.. math::

   \frac{d}{dt}\left ( 
      \begin{array}{c} \rho X_k \\ \rho E \\  \rho e \end{array} 
   \right ) = 
   \left ( \begin{array}{c}
      {\left [\mathcal{A} \left(\rho X_k\right)\right]}^{n+1/2} \\ {\left [\mathcal{A} \left(\rho E\right)\right]}^{n+1/2} \\ {\left [\mathcal{A} \left(\rho e\right)\right]}^{n+1/2} \\
   \end{array} \right ) +
   \left (
      \begin{array}{c} \rho \dot{\omega}_k \\ \rho \dot{S}\\ \rho \dot{S}\end{array}
   \right )

where we include :math:`e` in addition to :math:`E` to provide additional thermodynamic
information to help find a consistent :math:`T`. Here the advective courses
are piecewise-constant approximations to the change in the state due
to the hydrodynamics, computed with the during the hydro step.

However, to define the temperature, we need the kinetic energy (and
hence momentum and density) at any intermediate time, :math:`t`. We construct
these as needed from

Interfaces
==========

actual_integrator
-----------------

SDC integration provides its own implementation of the main entry
point, actual_integrator.

::

      subroutine actual_integrator(state_in, state_out, dt, time)

        type (sdc_t), intent(in   ) :: state_in
        type (sdc_t), intent(inout) :: state_out
        real(dp_t),    intent(in   ) :: dt, time

The main difference here is that the type that comes in and out of the
interface is no longer a burn_t, but instead is an
sdc_t.

The flow of this main routine is simpler than the non-SDC version:

#. Convert from the sdc_t type to the integrator’s internal
   representation (e.g., sdc_to_bs converts from a bs_t
   for the BS integrator).

   This copies the state variables and advective sources into the
   integration type. Since we only actually integrate :math:`(\rho X_k),
     (\rho e), (\rho E)`, the terms corresponding to density and momentum
   are carried in an auxillary array (indexed through the rpar
   mechanism).

#. Call the main integration routine to advance the inputs state
   through the desired time interval, producing the new, output state.

#. Convert back from the internal representation (e.g., a
   bs_t) to the sdc_t type.

Righthand side wrapper
----------------------

The manipulation of the righthand side is considerably more complex
now. Each network only provides the change in molar
fractions and internal energy, but
we need to convert these to the conservative system we are
integrating, including the advective terms.

#. Convert to a burn_t type, for instance via bs_to_burn:

   #. call fill_unevolved_variables to update the density
      and momentum. Since these don’t depend on reactions, this is a
      simple algebraic update based on the piecewise-constant-in-time
      advective term:

      .. math::

         \begin{aligned}
               \rho(t) &= \rho^n + t \cdot \left [ \mathcal{A}(\rho) \right]^{n+1/2} \\
               (\rho {\bf U})(t) &= (\rho {\bf U})^n + t \cdot \left [ \mathcal{A}(\rho{\bf U}) \right]^{n+1/2} 
             \end{aligned}

      where here we assume that we are integrating in the ODE system
      starting with :math:`t=0`.

   #. compute the specific internal energy, :math:`e`, from either the
      total minus kinetic energy or internal energy carried by the
      integrator (depending on the value of T_from_eden).

   #. get the temperature from the equation of state

   #. convert to a burn_t type, for instance via eos_to_burn:

#. Call the network’s actual_rhs routine to get just the
   reaction sources to the update. In particular, this returns
   the change in molar fractions, :math:`\dot{Y}_k` and the nuclear energy
   release, :math:`\dot{S}`.

#. Convert back to the integrator’s internal representation (e.g.,
   a bs_t, via rhs_to_bs)

   #. call fill_unevolved_variables

   #. fill the ydot array in the integrator type (e.g.,
      bs_t) with the advective sources that originally came into the
      intergrator through the sdc_t.

   #. Add the reacting terms. This is done as:

      .. math::

         \begin{aligned}
               \dot{y}_{\rho X_k} &= {\left [\mathcal{A} \left(\rho X_k\right)\right]}^{n+1/2} + \rho A_k \dot{Y}_k \\
               \dot{y}_{\rho e} &= {\left [\mathcal{A} \left(\rho e\right)\right]}^{n+1/2} +\rho \dot{S} \\
               \dot{y}_{\rho E} &= {\left [\mathcal{A} \left(\rho E\right)\right]}^{n+1/2} + \rho \dot{S}
             \end{aligned}
