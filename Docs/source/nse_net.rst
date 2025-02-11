.. _self_consistent_nse:

*******************
Self-consistent NSE
*******************

The self-consistent NSE approach uses only the nuclei in the main
reaction network.  It solves for the chemical potentials of the proton
and neutron and from there gets the abundances of each of the nuclei
under the assumption of NSE, following the procedure outlined in :cite:`Calder_2007`.

.. important::

   Self-consistent NSE does not support the templated C++ networks
   (like ``aprox13``).  You should use a pynucastro-generated network.

The solve is done using a port of the hybrid Powell method from
MINPACK (we ported the solver to templated C++).

The advantage of this approach is that it can be used with any
reaction network, once the integration has reached NSE.

This solver is enabled by compiling with

.. prompt:: bash

   USE_NSE_NET=TRUE

The functions to find the NSE state are then found in ``nse_solver.H``.

Dynamic NSE Check
=================

We have implemented a dynamic NSE check for the self-consistent nse procedure
that tells us whether the network has reached the NSE state.
The overall procedure is outlined in :cite:`Kushnir_2020`.
The overall usage comes down to a single function ``in_nse(state)``.
By supplying the current state, this function returns a boolean that tells us
whether we're in NSE or not. The current status of this functionality only works
for pynucastro-generated network since aprox networks have slightly
different syntax.

The overall framework is constructed following :cite:`Kushnir_2020` with slight
variations. The overview of the steps we take are the following:

* Minimum Temperature Check: require ``T > nse.T_min_nse``, where ``nse.T_min_nse`` is
  a runtime parameter with a default value ``nse.T_min_nse = 4.0e9``.

* Mass Abundance Check: compare the current mass abundances of the nuclei to
  the NSE mass fractions. A detailed criteria are the following:

  We first determine whether the current molar fraction is close to NSE
  with a criteria of:

  .. math::

     \frac{r - r_{NSE}}{r_{NSE}} < 0.5

  where :math:`r = Y_\alpha/(Y_p^2 Y_n^2)` and
  :math:`r_{NSE} = \left(Y_\alpha/(Y_p^2 Y_n^2)\right)_{NSE}` if there is
  neutron in the network.

  .. math::

     \frac{r - r_{NSE}}{r_{NSE}} < 0.25

  where :math:`r = Y_\alpha/(Y_p^2)` and
  :math:`r_{NSE} = \left(Y_\alpha/(Y_p^2)\right)_{NSE}` if neutron
  is not in the network.

  If the molar check above failed, then we proceed with an overall molar
  fraction check:

  .. math::

    \epsilon_{abs} = Y^i - Y^i_{NSE} < \mbox{nse.nse_abs_tol}

  .. math::

    \epsilon_{rel} = \frac{\epsilon_{abs}}{Y^i} < \mbox{nse.nse_rel_tol}

  where ``nse.nse_rel_tol = 0.2`` and ``nse.nse_abs_tol = 0.005`` by default.


* **Removed** :cite:`Kushnir_2020` also requires a fast reaction cycle that
  exchanges 1 :math:`\alpha` particle with 2 :math:`p` and 2 :math:`n`
  particles. We used to have this check, but currently removed as
  we think it is not necessary. However, the description is as following:
  This reaction cycle should have the following reactions or
  their reverse:

  * 1 :math:`(\alpha, \gamma)`, 2 :math:`(\gamma, p)`, 2 :math:`(\gamma, n)`
  * 1 :math:`(\alpha, p)`, 1 :math:`(\gamma, p)`, 2 :math:`(\gamma, n)`
  * 1 :math:`(\alpha, n)`, 2 :math:`(\gamma, p)`, 1 :math:`(\gamma, n)`

  To consider to be fast reaction cycle, every step in the cycle to have
  :math:`Y_i/\textbf{min}(b_f, b_r) < \epsilon t_s` for :math:`i = n, p, \alpha`
  participated in this step, where :math:`b_f` and :math:`b_r`
  are the forward and reverse rate of the reaction,
  :math:`\epsilon` is a tolerance which has a default value of
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

       t_{i,k} < \epsilon t_s

  *

    .. math::

      2|b_f(k) - b_r(k)|/(b_f(k) + b_r(k) < \epsilon

  Here we only consider two cases of reactions:

  * There are exactly two isotopes involved in reaction, :math:`k`,
    that are not in the light-isotope-group. In this case,
    if the reaction passes the two criteria mentioned above,
    we merge the groups containing those two isotopes if they're
    not yet in the same group.

  * There is only one isotope involved in reaction, :math:`k`,
    that is not in the light-isotope-group, which is not
    necessarily isotope :math:`i` that passes the first criteria.
    In this case, we merge the isotope that is not in LIG into LIG.

  Here we skip over reactions of the following due to obvious reasons:

  * Reactions that have no reverse rates.

  * Reactions that involve more than 2 reactants and products

  * Reactions that have more than 2 non-light-isotope-group.

  * The nuclei that participate in the reaction is either in LIG or in
    another group. This means that the non-LIG nuclei have already merged.

  At the end of the grouping process,
  we define that the current state have reached NSE
  when there is only a single group left, or there are two groups
  left where one of them is the light-isotope-group.

  When there is no neutron in the network, it can be difficult
  for isotopes to form a single group due to the missing neutron rates.
  Therefore, there is an alternative criteria of defining a "single group"
  when neutron is not present in the network: for isotopes,
  :math:`Z >= 14`, isotopes with odd and even :math:`N` form two
  distinct groups.


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
