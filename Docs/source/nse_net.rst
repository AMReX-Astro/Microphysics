.. _self_consistent_nse:

*******************
Self-consistent NSE
*******************

The self-consistent NSE approach uses only the nuclei in the main
reaction network.  It solves for the chemical potentials of the proton
and neutron and from there gets the abundances of each of the nuclei
under the assumption of NSE, following the procedure outlined in :cite:`Calder_2007`.

.. important::

   Self-consistent NSE is only supported for pynucastro-generated network.
   Templated C++ networks, such as ``aprox13``, are not supported,
   and compiling with them will produce a compilation error.

   Additionally, the pynucastro network must meet the following requirements:

   * Every rate should have an inverse rate, recomputed with detailed balancing

   * The "Q" value should be recomputed using the masses when doing detailed
     balance.

   * Approximate rates like $(\alpha, p)(p, \gamma)$ are okay, but modified
     rates that change the stoichiometry are not supported.

.. tip::

   The test in ``Microphysics/nse_solver/nse_compatibility/`` can be used to
   check if a network is supported.  It will integrate at fixed $(\rho, T)$
   and compare the resulting mass fractions to the NSE solution.

   If a network has weak rates, then there will be a $Y_e$ evolution which
   will prevent this test from giving perfect agreement, but if the $Y_e$
   evolution is slow, this test can still give useful results.

.. tip::

   Presently it is recommended to use the ``ase`` network.

The advantage of this approach is that it can be used with any reaction network,
given that the network satisifies the above conditions.

This solver is enabled by compiling with

.. prompt:: bash

   USE_NSE_NET=TRUE

All functions related to self-consistent NSE can be found in
the directory: ``nse_solver/``.

NSE Solver Modes
================

Unlike :ref:`tabulated_nse` where NSE abundances are pre-computed
and stored into a table, NSE abundances must be solved at runtime in this mode.
In :ref:`ch:nse`, we showed all three possible constraint equations that are
required to solve the system. In order to solve the system, a root-finding
algorithm is needed. There are two solvers for solving NSE mass abundance:

1. Hybrid Powell's method. This is the default method.
   This algorithm is taken from MINPACK and
   ported to templated C++.

2. Newton-Raphson (NR) method.
   This method is not as robust compared
   to the Hybrid Powell's method through our testings,
   and generally should *NOT* be used.
   But this can be enabled via runtime parameter:
   ``nse.use_hybrid_solver=0``.

To use the above methods, analytic Jacobian of the system is required.
Given the definition of :math:`\bar{A}` and :math:`\bar{Z}` as:

.. math::

   \bar{A} = \left(\sum_k \frac{X_k}{A_k} \right)^{-1}

and

.. math::

   \bar{Z} = \bar{A} Y_e

Along with the equation of state of form (true for Helmholtz EOS):

.. math::

   e(\rho, T, X_k) = e \left(\rho, T, \bar{A}(\mu_p^\mathrm{kin}, \mu_n^\mathrm{kin}, T), \bar{Z}(\mu_p^\mathrm{kin}, \mu_n^\mathrm{kin}, T) \right)

The Jacobian of different input modes are:

1. Given the input :math:`(\rho, T, Y_e)`, the :math:`2 \times 2`  Jacobian is given by:

   .. math::

      \begin{bmatrix}
      \frac{\partial f_\rho}{\partial \mu_p^{\mathrm{kin}}}  = \sum_k \frac{\partial X_k^{\mathrm{NSE}}}{\partial \mu_p^{\mathrm{kin}}} &
      \frac{\partial f_\rho}{\partial \mu_n^{\mathrm{kin}}}  = \sum_k \frac{\partial X_k^{\mathrm{NSE}}}{\partial \mu_n^{\mathrm{kin}}} \\
      \frac{\partial f_{Y_e}}{\partial \mu_p^{\mathrm{kin}}} = \sum_k \frac{Z_k}{A_k} \frac{\partial X_k^{\mathrm{NSE}}}{\partial \mu_p^{\mathrm{kin}}} &
      \frac{\partial f_{Y_e}}{\partial \mu_n^{\mathrm{kin}}} = \sum_k \frac{Z_k}{A_k} \frac{\partial X_k^{\mathrm{NSE}}}{\partial \mu_n^{\mathrm{kin}}}
      \end{bmatrix}


2. Given the input :math:`(\rho, e, Y_e)`, the :math:`3 \times 3` Jacobian is given by:

   .. math::

     \begin{bmatrix}
     \frac{\partial f_\rho}{\partial \mu_p^{\mathrm{kin}}}  = \sum_k \frac{\partial X_k^{\mathrm{NSE}}}{\partial \mu_p^{\mathrm{kin}}} &
     \frac{\partial f_\rho}{\partial \mu_n^{\mathrm{kin}}}  = \sum_k \frac{\partial X_k^{\mathrm{NSE}}}{\partial \mu_n^{\mathrm{kin}}} &
     \frac{\partial f_\rho}{\partial T_9} = \sum_k \frac{\partial X_k^{\mathrm{NSE}}}{\partial T_9} \\
     \frac{\partial f_{Y_e}}{\partial \mu_p^{\mathrm{kin}}} = \sum_k \frac{Z_k}{A_k} \frac{\partial X_k^{\mathrm{NSE}}}{\partial \mu_p^{\mathrm{kin}}} &
     \frac{\partial f_{Y_e}}{\partial \mu_n^{\mathrm{kin}}} = \sum_k \frac{Z_k}{A_k} \frac{\partial X_k^{\mathrm{NSE}}}{\partial \mu_n^{\mathrm{kin}}} &
     \frac{\partial f_{Y_e}}{\partial T_9} = \sum_k \frac{Z_k}{A_k} \frac{\partial X_k^{\mathrm{NSE}}}{\partial T_9} \\
     \frac{\partial f_e}{\partial \mu_p^{\mathrm{kin}}} = \frac{1}{e_0} \frac{\partial e}{\partial \mu_p^{\mathrm{kin}}} &
     \frac{\partial f_e}{\partial \mu_n^{\mathrm{kin}}} = \frac{1}{e_0} \frac{\partial e}{\partial \mu_n^{\mathrm{kin}}} &
     \frac{\partial f_e}{\partial T_9} = \frac{1}{e_0} \frac{\partial e}{\partial T_9} \\
     \end{bmatrix}


Partial derivatives appear in the Jacobian can be shown as the following:

1. :math:`\frac{\partial X_k^{\mathrm{NSE}}}{\partial \mu_p^{\mathrm{kin}}} = \frac{Z_k}{k_B T} X_k^{\mathrm{NSE}}`
2. :math:`\frac{\partial X_k^{\mathrm{NSE}}}{\partial \mu_n^{\mathrm{kin}}} = \frac{N_k}{k_B T} X_k^{\mathrm{NSE}}`
3. :math:`\frac{\partial X_k^{\mathrm{NSE}}}{\partial T_9} = X_k^{\mathrm{NSE}} \times \left(\frac{d \log{G_k}}{d T_9} + \frac{1}{T_9} \left[\frac{3}{2} - \frac{Z_k \mu_p^\mathrm{kin} + N_k \mu_n^\mathrm{kin} + B_k}{k_B T}\right] - \frac{\partial}{\partial T_9} \left[\frac{\mu^c_k}{k_B T}\right]\right)`
4. :math:`\frac{\partial \bar{A}}{\partial \mu_p^{\mathrm{kin}}} = - \bar{A}^2 \left(\frac{1}{A_k} \frac{\partial X_k^{\mathrm{NSE}}}{\partial \mu_p^{\mathrm{kin}}}\right)`
5. :math:`\frac{\partial \bar{A}}{\partial \mu_n^{\mathrm{kin}}} = - \bar{A}^2 \left(\frac{1}{A_k} \frac{\partial X_k^{\mathrm{NSE}}}{\partial \mu_n^{\mathrm{kin}}}\right)`
6. :math:`\frac{\partial \bar{A}}{\partial T_9} = - \bar{A}^2 \left(\frac{1}{A_k} \frac{\partial X_k^{\mathrm{NSE}}}{\partial T_9}\right)`
7. :math:`\frac{\partial e}{\partial \mu_p^{\mathrm{kin}}} = \left( \frac{\partial e}{\partial \bar{A}} + Y_e \frac{\partial e}{\partial \bar{Z}} \right) \frac{\partial \bar{A}}{\partial \mu_p^{\mathrm{kin}}}`
8. :math:`\frac{\partial e}{\partial \mu_n^{\mathrm{kin}}} = \left( \frac{\partial e}{\partial \bar{A}} + Y_e \frac{\partial e}{\partial \bar{Z}} \right) \frac{\partial \bar{A}}{\partial \mu_n^{\mathrm{kin}}}`
9. :math:`\frac{\partial e}{\partial T_9} = \left( \frac{\partial e}{\partial \bar{A}} + Y_e \frac{\partial e}{\partial \bar{Z}} \right) \frac{\partial \bar{A}}{\partial T_9} + \frac{\partial e}{\partial T_9}`

Note that without the loss of generality, we treated
:math:`\mu_p^\mathrm{kin} + \mu_p^c \rightarrow \mu_p^\mathrm{kin}`.
Note a similar procedure is done with :cite:`lippuner_skynet_2017`, although they
work with entropy instead of internal energy.


Dynamic NSE Check
=================

We have implemented a dynamic NSE check for the self-consistent NSE procedure
following :cite:`Kushnir_2020` with slight variations.
The overall usage comes down to a single function ``in_nse(state)``.
By supplying the current state, this function returns a boolean that tells us
whether we're in NSE or not.

The overview of the steps we take are the following:

* Minimum Temperature Check: require ``T > nse.T_min_nse``,
  where ``nse.T_min_nse`` is a runtime parameter with
  a default value ``nse.T_min_nse = 4.0e9``.

* Mass Abundance Check: compare the current mass abundances of the nuclei to
  the NSE mass fractions. If neutron, hydrogen, and Helium-4, i.e.
  :math:`(n, p, \alpha)` are all present in the network, then follow
  :cite:`Kushnir_2020` with a criteria of

  .. math::

     \frac{r - r_\mathrm{NSE}}{r_\mathrm{NSE}} < 0.5

  where :math:`r = Y_\alpha/(Y_p^2 Y_n^2)` and
  :math:`r_\mathrm{NSE} = \left(Y_\alpha/(Y_p^2 Y_n^2)\right)_\mathrm{NSE}`.

  If the molar check above failed or not all :math:`(n, p, \alpha)`
  are present, then we proceed with an overall molar fraction check:

  .. math::

    \epsilon_\mathrm{rel} = \frac{Y^i - Y^i_\mathrm{NSE}}{Y^i_\mathrm{NSE}}

  where :math:`\epsilon_\mathrm{rel}` is set by ``nse.nse_rel_tol``,
  which has a default value of 0.5.

.. note::

   We require Helium-4 to be present in the network for NSE_NET.

* **Removed** :cite:`Kushnir_2020` also requires a fast reaction cycle that
  exchanges 1 :math:`\alpha` particle with 2 :math:`p` and 2 :math:`n`
  particles. We used to have this check, but currently removed as
  we think it is not necessary. However, the description is the following:
  The reaction cycle should have the following reactions and
  their reverse:

  * 1 :math:`(\alpha, \gamma)`, 2 :math:`(\gamma, p)`, 2 :math:`(\gamma, n)`
  * 1 :math:`(\alpha, p)`, 1 :math:`(\gamma, p)`, 2 :math:`(\gamma, n)`
  * 1 :math:`(\alpha, n)`, 2 :math:`(\gamma, p)`, 1 :math:`(\gamma, n)`

  To consider to be a fast reaction cycle, every reaction rate in the cycle
  needs to have :math:`Y_i/\textbf{min}(b_f, b_r) < \epsilon \ t_s` for
  :math:`i = \{n, p, \alpha\}` participated in this step,
  where :math:`b_f` and :math:`b_r` are the forward and reverse rate of the
  reaction. :math:`\epsilon` is a tolerance which has a default value of
  :math:`0.1`, and :math:`t_s` is the sound crossing time of a simulation cell.

  An example of such reaction cycle would be:

  .. math::

     \isotm{S}{32} (\gamma, p)(\gamma, p)(\gamma, n)(\gamma, n) \isotm{Si}{28}
     (\alpha, \gamma) \isotm{S}{32}

* NSE Grouping Process: Initially, :math:`p`, :math:`n`, and
  :math:`\alpha` are grouped into a single group
  called the light-isotope-group, or LIG. Other isotopes belong to their
  own group, which only contains themselves. We need to start the grouping
  process with the reaction rate that has the fastest (smallest) timescale.
  In the original :cite:`Kushnir_2020` paper, they use the group molar fraction
  for evaluating the reaction timescale. This complicates things because
  now reaction timescale changes after each successful grouping. We've
  determined that the result is roughly the same even if we just use the
  molar fraction of the isotope that is involved in the actual reaction.
  Therefore, instead of using
  :math:`t_{i,k} = \tilde{Y}_i/\textbf{min}(b_f(k), b_r(k))`, to evaluate
  the reaction timescale of the reaction, :math:`k`, where
  :math:`\tilde{Y}_i` represents the sum of molar fractions of the
  group that isotope :math:`i` belongs to, we simply use the :math:`Y_i`,
  which is the molar fraction of the isotope :math:`i`, which is the
  isotope involved in the reaction that is different from
  :math:`p`, :math:`n`, and :math:`\alpha`. After we settle on calculating
  the timescale, since :math:`Y_i` doesn't change, we can calculate all
  timescale at once and sort the reaction to determine the order at
  which we want to start merging.

  There are two requirements for us to check whether this reaction
  can be used to group the nuclei involved, which are:

  * at least 1 isotope, :math:`i`, that passes:

    .. math::

       t_{i,k} < \epsilon \ t_s

  * the forward and reverse rates satisfy:

    .. math::

      2|b_f(k) - b_r(k)|/(b_f(k) + b_r(k) < \epsilon

    where :math:`\epsilon` is set by ``nse.ase_tol``,
    which has a default value of 0.1.

  Here we only consider two cases for merging:

  1. There are exactly three distinct groups: the
     light-isotope-group and *TWO* other nuclei from TWO different groups.
     In this case, merge those two nuclei into the same group.

  2. There are exactly TWO distinct groups: light-isotope-group and ONE
     or TWO other nuclei but belong to the *same* group.
     In this case, merge those two groups into the same group.

  At the end of the grouping process, the system is considered to have reached
  NSE when only a single group remains.

.. note::
   For NSE grouping, reactions are pre-filtered using the following criteria:

   1. The reaction must have both forward and reverse rates.
   2. Weak reactions and removed rates are excluded.
   3. The reaction may have at most two reactants and two products,
      except for the triple-alpha reaction.
   4. Across all reactants and products, there must be only one or two nuclei
      not in :math:`(p, n, \alpha)`.


NSE Reactive Update
===================

Here we give an outline on how the overall reactive update is done with NSE:

* Check if NSE applies via ``in_nse``

  * If we are in an NSE region:

    * Compute the advective source term for $Y_e$ via:

      .. math::

         \Advs{\rho Y_e} = \sum_k \frac{Z_k}{A_k} \Advs{\rho X_k}

    * Compute $[\Rb(\rho Y_e)]^n$ and $[\Rb(\rho e)_{\mathrm{nuc}}]^n$ using
      $[\rho]^n$, $[T]^n$, $[Y_e]^n$ and $[e]^n$:

      * Find the NSE composition with given $[\rho]^n$, $[e]^n$,
        and $[Y_e]^n$. There are two approaches that can be selected at runtime:

        1. ``nse.nse_solve_e_mode = 1``: An EOS inversion algorithm is used so
           that we determine $[T^*]^n$ such that $[e]^n$ remains
           unchanged after switching to the NSE composition
           via $([\rho]^n, [T^*]^n, [Y_e]^n)$ solver mode.
           Here we use $[T]^n$ as the initial guess and updated
           to the solution, $[T^*]^n$, in the end. This can be more robust
           but slower compared to the second approach.

        2. ``nse.nse_solve_e_mode = 2``: Find the NSE composition directly with
           $([\rho]^n, [e]^n, [Y_e]^n)$ by using the :math:`3 \times 3`
           Jacobian matrix defined above. By using this mode, we directly solve
           :math:`T` together with the chemical potentials. It is faster,
           but slightly less robust.

      * Compute the thermal neutrino losses,
        $\epsilon_{\nu,\mathrm{thermal}}$, using the NSE composition.

      * Evaluate $\dot{Y}_{\mathrm{weak}}$ and neutrino losses,
        $\epsilon_{\nu,\mathrm{react}}$,
        from weak reactions only as they are the only contributing
        reactions in NSE.

      * Evaluate $[\Rb(\rho Y_e)]^n$ as:

        .. math::
           [\Rb(\rho Y_e)]^n = [\rho]^n \sum_k Z_k [\dot{Y}_{\mathrm{k, weak}}]^n

      * Evaluate $[\Rb(\rho e)_{\mathrm{nuc}}]^n$ as:

        .. math::
           [\Rb(\rho e)_{\mathrm{nuc}}]^n = - N_A c^2 \sum_k [\dot{Y}_{\mathrm{k, weak}}]^n m_k

        where the nuclei mass, $m_k$ is defined as:

        .. math::
           m_k c^2 = (A_k - Z_k) m_n c^2 + Z_k (m_p + m_e) c^2 - B_k

      * The full reactive source term, $[\Rb(\rho e)]^n$ is then:

        .. math::
           [\Rb(\rho e)]^n = [\Rb(\rho e)_{\mathrm{nuc}}]^n - [\rho]^n \left(\epsilon_{\nu,\mathrm{thermal}} + \epsilon_{\nu,\mathrm{react}}\right)

    * Now evolve $\rho$, $\rho e$, and $\rho Y_e$ to midpoint in time:

      .. math::
         \Uc^{\prime,n+1/2} = \Uc^{\prime,n} + \frac{\Delta t}{2} \left([\Advs{\Uc^\prime}]^{n+1/2} + [\Rb(\Uc^\prime)]^{n}\right)

      Note that there is no reactive source term for $\rho$ and the advective
      source term is already time-centered from the simplified-SDC algorithm.

    * Compute $[\Rb(\rho Y_e)]^{n+1/2}$ and
      $[\Rb(\rho e)_{\mathrm{nuc}}]^{n+1/2}$ following the same
      procedure as above. This time, it uses
      $[\rho]^{n+1/2}$, $[Y_e]^{n+1/2}$ and $[e]^{n+1/2}$ as input
      and uses the updated $[T]^n$ as initial guess for the EOS inversion
      algorithm.

    * Now evolve all thermodynamic quantities to new time, $t^{n+1}$, using the
      midpoint reactive source terms constructed in the previous step.

      .. math::

         \Uc^{\prime,n+1} = \Uc^{\prime,n} + \Delta t \left([\Advs{\Uc^\prime}]^{n+1/2} + [\Rb(\Uc^\prime)]^{n+1/2}\right)

    * Lastly, the composition is updated by finding the corresponding NSE state
      using $[\rho]^{n+1}$, $[e]^{n+1}$, and $[Y_e]^{n+1}$

  * If we are not in an NSE region:

    * Integrate the network as usual.


Additional Options
==================

.. index:: nse.nse_dx_independent, nse.nse_molar_independent, nse.nse_skip_molar, nse.T_nse_net, nse.ase_tol, nse.nse_abs_tol, nse.nse_rel_tol, nse.T_min_nse

Here we have some runtime options to allow a more cruel estimation
to the self-consistent nse check:

* ``nse.nse_dx_independent = 1`` in the input file allows the nse check
  to ignore the dependency on the cell size, ``dx``, which calculates
  the sound crossing time, ``t_s``. Naturally, we require the
  timescale of the rates to be smaller than ``t_s`` to ensure the
  states have time to achieve equilibrium. However, sometimes this
  check can be difficult to achieve, so we leave this as an option
  for the user to explore.

* ``nse.nse_molar_independent = 1`` in the input file allows the
  user to use the nse mass fractions for nse check after the first
  check (the one that ensures we're close enough to the nse mass fractions
  to get reasonable results) is passed. This allows the subsequent checks
  to only rely on the thermodynamic conditions instead of mass fractions.

* ``nse.nse_skip_molar = 1`` in the input file allows the user to skip
  the molar fraction check after the integration has failed.
  This option is used to completely forgo the requirement on molar
  fractions and allow the check to only dependent on the thermodynamic
  conditions. By only applying this after option after the
  integration failure, we hope the integrator has evolved the
  system to the NSE state the best it can. By turning on this option,
  we hope to give relief to the integrator if the system is in
  NSE thermodynamically,  which is likely the case.

* ``nse.T_nse_net`` in the input file allows the user to define a simple
  temperature threshold to determine the NSE state instead of using
  the complicated procedure that looks for a balance between the
  forward and the reverse rates. Once this quantity is set to a positive
  value, then ``in_nse`` returns ``true`` if the current temperature
  is higher than ``T_nse_net``, and ``false`` if the current
  temperature is lower than ``T_nse_net``.
  Note that we still perform a simple molar fraction check to
  ensure that the current state is close enough to the NSE state.

* ``nse.ase_tol`` is the tolerance that determines the equilibrium
  condition for forward and reverse rates. This is set to 0.1 by default.

* ``nse.nse_abs_tol`` is the absolute tolerance of checking the difference
  between current molar fraction and the NSE molar fraction.
  This is set to 0.005 by default.

* ``nse.nse_rel_tol`` is the relative tolerance of checking the
  difference between current molar fraction and the NSE molar fraction.
  This is set to 0.2 by default.

* ``nse.T_min_nse`` is the minimum temperature required to consider
  the subsequent NSE checks. This is mainly to avoid unnecessary computations
  of computing the NSE mass fractions when the current temperature is too low.
  This is set to 4.0e9 by default.
