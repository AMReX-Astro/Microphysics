**************************
Dealing with Burn Failures
**************************

Sometimes the ODE integration of a reaction network will fail.  Here
we summarize some wisdom on how to avoid integrator failures.

Error codes
===========

The error code that is output when the integrator fails can be
interpreted from the ``enum`` in ``integrator_data.H``:

* ``-1`` : bad inputs -- this indicates that the integrator was setup incorrectly.
  At the moment, this is mainly a check on the tolerances.

* ``-2`` : the timestep underflowed.  This can happen if the integrator
  needs to repeatedly cut the timestep to try to make the specified
  tolerances.

* ``-3`` :  for the RKC integrator, this means that the power method used
  to estimate the spectral radius failed to converge.

* ``-4`` : too many steps were needed (the limit is controled by
  ``integrator.ode_max_steps``

* ``-5`` : too much accuracy was requested (this indicates a problem with
  the tolerances).

* ``-6`` : the nonlinear corrector used to solve the implicit update failed
  to converge.

* ``-7`` : the LU decomposition of the Jacobian failed.

* ``-8`` : the solution in the corrector violates positivity checks.

The most common errors are ``-2`` (timestep underflow) and ``-4`` (too many steps required).

Why does the integrator struggle?
=================================

There are a few common reasons why the integrator might encounter trouble:

* The integrators don't know that the mass fractions should stay in
  $[0, 1]$.

  Note: they should ensure that the mass fractions sum to 1 as long as
  $\sum_i \dot{\omega}_k = 0$ in the righthand side function.

* The state wants to enter nuclear statistical equilibrium.  As we
  approach equilibrium, the integrator will see large, oppositely
  signed flows from one step to the next, which should cancel, but
  numerically, the cancellation is not perfect.

* The Jacobian is not good enough to allow the nonlinear solver to
  converge.

  The Jacobians used in our networks are approximate (even the
  analytic one).  For example, the analytic Jacobian neglects the
  composition dependence in the screening functions.  In the analytic we
  neglect the composition influence in screening.


Making the integration robust
=============================

Some tips for helping the integrator:

* Use a tight absolute tolerance for the species
  (``integrator.atol_spec``).  This ensures that even small species
  are tracked well, and helps prevent generating negative mass
  fractions.

* Reduce ``integrator.species_failure_tolerance``.  This is used to
  determine whether a step should be rejected.  If any of the mass
  fractions are more than ``integrator.species_failure_tolerance``
  outside of $[0, 1]$, then the integrator will reject the step and
  try again (this is implemented in ``VODE`` and ``BackwardEuler``).

  If the value is too large (like ``0.01``), then this could allow the
  mass fractions to grow too much outside of $[0, 1]$, and then a
  single step rejection is not enough to recover.  Setting this value to
  be close to the typical scale of a mass fraction that we care about
  (closer to ``integrator.atol_spec``) can help.

* Use the burn retry logic.  By setting
  ``integrator.use_burn_retry=1``, the burn will immediately be
  retried if it fails, with a slightly different configuration.

  By default the type of Jacobian is swapped (if we were analytic,
  then do numerical, and vice vera).  This is controlled by
  ``integrator.retry_swap_jacobian``.  The tolerances can also be
  changed on retry.

* Use the right integrator.  In general, VODE is the best choice.
  But if the network is only mildly stiff, then RKC can work well
  (typically, it works when the temperatures are below $10^9~\mathrm{K}$.

* If you are near NSE, then use the NSE solver.

Things we no longer recommend:

* Clipping the species (``integrator.do_species_clip``) can lead to
  instabilities.  This changes the integration state directly in the
  righthand side function, outside of the control of the integrator.
  While it sometimes may work, it often leads to problems.  A better
  way to deal with keeping the species in $[0, 1]$ is through the
  absolute tolerance.

  See the SUNDIALS docs -- in particular, they say the modifying the
  state in the RHS function is bad.  This means don't use
  _species_clip.

Debugging a burn failure
========================


parse script
