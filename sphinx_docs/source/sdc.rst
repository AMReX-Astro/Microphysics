*****************************
Spectral Deferred Corrections
*****************************

Introduction
============

Spectral deferred correction (SDC) methods strongly couple reactions
and hydrodynamics, eliminating the splitting error that arises with
Strang (operator) splitting.  Microphysics supports two different
SDC formulations.

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


.. note::

   Application codes can set the make variables ``USE_TRUE_SDC`` or
   ``USE_SIMPLIFIED_SDC``.  But in Microphysics, both SDC formulations
   are supported by the same integrators, and both of these options
   will set the ``SDC`` preprocessor flag.

"True" SDC
----------

The true SDC implementation is described in :cite:`castro_sdc`.  It divides
the timestep into temporal nodes and uses low-order approximations to update
from one temporal node to the next.  Iteration is used to increase the order of accuracy.

The update from one temporal node, $m$, to the next, $m$, for iteration
$k+1$ appears as:

.. math::

   \begin{align*}
    \Uc^{m+1,(k+1)} = \Uc^{m,(k+1)}
     &+ \delta t_m\left[\Advs{\Uc^{m,(k+1)}} - \Advs{\Uc^{m,(k)}}\right] \\
     &+ \delta t_m\left[\Rbs{\Uc^{m+1,(k+1)}} - \Rbs{\Uc^{m+1,(k)}}\right]\\
     &+ \int_{t^m}^{t^{m+1}}\Advs{\Uc^{(k)}} + \Rbs{\Uc^{(k)}}d\tau.
   \end{align*}

Solving this requires a nonlinear solve of:

.. math::

   \Uc^{m+1,(k+1)} - \delta t_m \Rbs{\Uc}^{m+1,(k+1)} = \Uc^{m,(k+1)} + \delta t_m {\bf C}

where the right-hand side is constructed only from known states, and we
define ${\bf C}$ for convenience as:

.. math::

   \begin{align}
   {\bf C} &= \left [ {\Advs{\Uc}}^{m,(k+1)} - {\Advs{\Uc}}^{m,(k)} \right ]
                  -  {\Rbs{\Uc}}^{{m+1},(k)}  \nonumber \\
               &+ \frac{1}{\delta t_m} \int_{t^m}^{t^{m+1}} \left  ( {\Advs{\Uc}}^{(k)} + {\Rbs{\Uc}}^{(k)}\right ) d\tau
   \end{align}

This can be cast as an ODE system as:

.. math::

  \frac{d\Uc}{dt} \approx \frac{\Uc^{m+1} - \Uc^m}{\delta t_m} = \Rbs{\Uc} + {\bf C}

Simplified SDC
--------------

The Simplified-SDC method uses the ideas of the SDC method above, but instead
of dividing time into discrete temporary notes, it uses a piecewise-constant-in-time
approximation to the advection update over the timestep (for instance, computed by the CTU PPM method) and solves the ODE system:

.. math::

  \frac{\partial \Uc}{\partial t} = [\Advs{\Uc}]^{n+1/2} + \Rb(\Uc)

and uses iteration to improve the advective term based on the
reactions.  The full details of the algorithm are presented in
:cite:`castro_simple_sdc`.

Common ground
-------------

Both SDC formulations result in an ODE system of the form:

.. math::

   \frac{d\Uc^\prime}{dt} = \Rb(\Uc^\prime) + \mbox{advective terms}

where $\Uc^\prime$ is only

.. math::

   \Uc^\prime = \left ( \begin{array}{c} \rho X_k \\ \rho e \end{array} \right )

since those are the only components with reactive sources.
The ``SDC`` burner includes storage for the advective terms.

Interface and Data Structures
=============================

For concreteness, we use the simplified-SDC formulation in the description below,
but the true-SDC version is directly analogous.

``burn_t``
----------

To accommodate the increased information required to evolve the
coupled system, the ``burn_t`` used for Strang integration is extended
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
   This converts the ``ydot`` to mass fractions and adds the advective terms
   to ``ydot``.

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

.. note::

   In the original "true SDC" paper :cite:`castro_sdc`, the matrix
   system was more complicated, and we included density in ${\bf w}$.
   This is not needed, and we use the Jacobian defined in
   :cite:`castro_simple_sdc` instead.
