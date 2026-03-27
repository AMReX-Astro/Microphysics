**************************
Dealing with Burn Failures
**************************

Sometimes the ODE integration of a reaction network will fail.  Here
we summarize some wisdom on how to avoid integrator failures.

Error codes
===========

The error code that is output when the integrator fails can be
interpreted from the ``enum`` in ``integrator_data.H``.  See
:ref:`sec:error_codes` for a list of the possible codes.

The most common errors are ``-2`` (timestep underflow) and ``-4`` (too
many steps required).  A timestep underflow usually means that
something has gone really wrong in the integration and the integrator
keeps trying to cut the timestep to be able to advance.  Too many
steps means that the integrator hit the cap imposed by
``integrator.ode_max_steps``.  This should never be made too
large---usually if you need more than a few 1000 steps, then the some
of the solutions discussed below can help.



Why does the integrator struggle?
=================================

There are a few common reasons why the integrator might encounter trouble:

* The integrators don't know that the mass fractions should stay in
  $[0, 1]$.

  Note: they should ensure that the mass fractions sum to 1 as long as
  $\sum_i \dot{\omega}_k = 0$ in the righthand side function.

  This is a place where adjusting ``integrator.atol_spec`` can help.

* The state wants to enter nuclear statistical equilibrium.  As we
  approach equilibrium, the integrator will see large, oppositely
  signed flows from one step to the next, which should cancel, but
  numerically, the cancellation is not perfect.

  For some networks, the solution here is to use the NSE solver.

* The Jacobian is not good enough to allow the nonlinear solver to
  converge.

  The Jacobians used in our networks are approximate (even the
  analytic one).  For example, the analytic Jacobian neglects the
  composition dependence in the screening functions.  In the analytic we
  neglect the composition influence in screening.


Making the integration robust
=============================

.. index:: integrator.atol_spec, integrator.species_failure_tolerance, integrator.use_jacobian_caching, integrator.do_corrector_validation, integrator.use_burn_retry

Some tips for helping the integrator:

* Use a tight absolute tolerance for the species
  (``integrator.atol_spec``).  This ensures that even small species
  are tracked well, and helps prevent generating negative mass
  fractions.

  In general, ``atol_spec`` should be picked to be the magnitude of
  the smallest mass fraction you care about.  You should never set
  it to be larger that ``rtol_spec``.  See :ref:`sec:tolerances`
  for some discussion.

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

  See :ref:`sec:retry` for the full list of options.

* Use the right integrator.  In general, VODE is the best choice.
  But if the network is only mildly stiff, then RKC can work well
  (typically, it works when the temperatures are below $10^9~\mathrm{K}$.

* If you are near NSE, then use the NSE solver.  This is described
  in :ref:`self_consistent_nse`.

  Note: not every network is compatible with the self-consistent
  NSE solver.

* Use Jacobian-caching.  If you build on GPUs, this is disabled by
  default.  You can re-enable it by building with
  ``USE_JACOBIAN_CACHING=TRUE``.  Also make sure that
  ``integrator.use_jacobian_caching=1`` (this is the default).

  By reducing the number of times the Jacobian is evaluated, we also
  reduce the possibility of trying to evaluate it with a bad state.

* Use the corrector validation (``integrator.do_corrector_validation``).

  This checks to make sure the state is valid inside of the corrector
  loop, and if not, it bails out of the corrector, forcing the
  integrator to retry the entire step.

Things we no longer recommend:

.. index:: integrator.do_species_clip, integrator.renormalize_abundances

* Clipping the species (``integrator.do_species_clip``) can lead to
  instabilities.  This changes the integration state directly in the
  righthand side function, outside of the control of the integrator.
  While it sometimes may work, it often leads to problems.  A better
  way to deal with keeping the species in $[0, 1]$ is through the
  absolute tolerance.

  The `SUNDIALS CVODE documentation <https://sundials.readthedocs.io/en/latest/cvode/Usage/index.html#advice-on-controlling-unphysical-negative-values>`_ has
  a good discussion on this numerical instability.

* Renormalizing the species during the integration / righthand side call
  (``integrator.renormalize_abundances``).  Like clipping, this can
  cause numerical instabilities.

Debugging a burn failure
========================

When a burn fails, the entire ``burn_t`` state will be output to
stdout (CPU runs only).  This state can then be used with
``burn_cell_sdc`` to reproduce the burn failure outside of a
simulation.  For SDC, see the discussion at :ref:`sec:redo_burn_fail`.
