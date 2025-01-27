***
NSE
***

.. important::

   NSE is only supported with the simplified-SDC method for
   coupling hydrodynamics and reactions.  We do not support
   operator-splitting (Strang) coupling with NSE.

The reaction networks in Microphysics have the ability to use NSE
instead of integrating the entire network when the conditions are
appropriate.  There are 2 different implementations of NSE in
Microphysics, that have slightly different use cases.

.. index:: USE_NSE_TABLE, USE_NSE_NET

* :ref:`tabulated_nse` : this uses a table of NSE abundances given
  :math:`(\rho, T, Y_e)` generate from a large network (96 isotopes).
  The table also returns :math:`dY_e/dt` resulting from
  electron-captures, to allow for the NSE state to evolve.  This is
  meant to be used in the cores of massive stars and works only with the
  ``aprox19`` reaction network.

  Furthermore, since the table can achieve :math:`Y_e` and
  :math:`\bar{A}` that are not representable by the 19 isotopes in
  ``aprox19``, this table requires that we use the auxiliary
  composition and advect :math:`Y_e`, :math:`\bar{A}`, and
  :math:`\langle B/A\rangle`.  All of the EOS calls will work with
  these quantities.

  This algorithm was described in :cite:`sdc-nse`.

  This is enabled via ``USE_NSE_TABLE``


* :ref:`self_consistent_nse` : this adds an NSE solver to the network that
  can be called to find the equilibrium abundances of each of the
  species defined in the network.  It works with any of the
  pynucastro-generated networks.  Unlike the tabulated NSE, there is
  no need to advect the auxiliary composition, since this only deals
  with the isotopes defined in the main reaction network.

  This is enabled via ``USE_NSE_NET``

Both solvers define a number of preprocessor variables, and both will
provide a function ``in_nse()`` that can be used to determine if a
state is currently in NSE.

=================        ======================================
make option               preprocessor variables set
-----------------        --------------------------------------
``USE_NSE_TABLE``        ``NSE``, ``NSE_TABLE``, ``AUX_THERMO``
``USE_NSE_NET``          ``NSE``, ``NSE_NET``
=================        ======================================

The directive ``NSE`` should be used whether the specific
implementation of NSE does not matter.

These two NSE solvers are described below.


.. _tabulated_nse:

Tabulated NSE and ``aprox19``
=============================

The ``aprox19`` network can be run in a manner where we blends the
standard ``aprox19`` network with a table for nuclear statistic
equilibrium resulting from a much larger network at high density and
temperatures.    This option is enabled by building with:

.. code:: bash

   NETWORK_DIR=aprox19 USE_NSE_TABLE=TRUE

Composition and EOS
-------------------

The NSE table was generated using `pynucastro
<https://pynucastro.github.io/pynucastro/>` using 96 nuclei and
electron/positron capture/decay rates from :cite:`langanke:2001`.  The
table takes $Y_e$ as the primary composition variable and provides a
set of mass fractions that is mapped into those used by ``aprox19``.
Using the value allows us to attain a lower :math:`Y_e` than
``aprox19`` can represent.

For this reason, when we are using the NSE network, we always take the
composition quantities in the EOS directly from ``eos_state.aux[]``
instead of from ``eos_state.xn[]``.  The ``AUX_THERMO`` preprocessor
variable is enabled in this case, and the equations of state interpret
this to use the auxiliary data for the composition.  This is described in :ref:`aux_eos_comp`.


NSE Table Outputs
-----------------

The NSE table provides values for the auxiliary composition,
:math:`\bar{A}`, and :math:`\langle B/A \rangle`
resulting from the full 96 nuclei network.   It also provides a set of 19
:math:`X_k` that map into the isotopes carried by ``aprox19``.
These three quantities are stored as ``aux`` data in the network and
are indexed as ``iye``, ``iabar``, and ``ibea``.  Additionally, when
coupling to hydrodynamics, we need to advect these auxiliary
quantities.

The evolution equations for the auxiliary variables are:

.. math::

   \begin{align*}
   \frac{DY_e}{Dt} &= \sum_k \frac{Z_k}{A_k} \dot{\omega}_k \\
   \frac{D\bar{A}}{Dt} &= -\bar{A}^2 \sum_k \frac{1}{A_k} \dot{\omega}_k \\
   \frac{D}{Dt} \left (\frac{B}{A} \right ) &= \sum_k \frac{B_k}{A_k} \dot{\omega}_k
   \end{align*}

Therefore each of these auxiliary equations obeys an advection equation
in the hydro part of the advancement.

The table also provides $dY_e/dt$, $(d\langle
B/A\rangle/dt)_\mathrm{weak}$, and $\epsilon_{\nu,\mathrm{react}}$, the
weak rate neutrino losses.  These quantities are used to update the
thermodynamic state as we integrate.

NSE Flow
--------

.. index:: integrator.nse_deriv_dt_factor, integrator.nse_include_enu_weak

The time integration algorithm is described in detail in :cite:`sdc-nse`.  Here
we provide an outline:

* initialize the problem, including :math:`X_k`

* fill the initial aux data with :math:`Y_e`, :math:`\bar{A}`, and :math:`(B/A)`

* in hydro, we will update these quantities simply via advection (for
  Strang-split evolution)

* for the reactive update:

  * check if NSE applies (see below)

  * if we are in an NSE region:

    * Compute the initial temperature given $\rho$, $e$, and $Y_e$,
      using an EOS inversion algorithm that understands NSE (in
      particular that the composition depends on $T$ in NSE)

    * Compute the plasma neutrino losses, $\epsilon_{\nu,\mathrm{thermal}}$

    * Use :math:`\rho`, :math:`T`, and :math:`Y_e` to evaluate the NSE
      state and construct $[\Rb(\Uc^\prime)]^n$, the source term from reactions to the
      reduced conserved state $\Uc^\prime$ (this is the state used by the SDC algorithm
      and includes the internal energy density, mass fractions, and auxiliary variables).

      This is done via finite differencing in time (through a step
      $\tau \ll \Delta t$), and the reactive sources are constructed
      to exclude the advective contributions.  The size of $\tau$ is
      controlled via ``integrator.nse_deriv_dt_factor``.

      In particular, the energy source is constructed as:

      .. math::

         R(\rho e) = N_A \frac{\Delta (\rho \langle B/A\rangle)}{\tau} + N_A \Delta m_{np} c^2 \rho \frac{dY_e}{dt} - \rho (\epsilon_{\nu,\mathrm{thermal}} + \epsilon_{\nu,\mathrm{react}})

      where $\Delta m_{np}$ is the difference between the neutron and H atom mass.

      .. important::

         It only makes sense to include the weak rate neutrino losses, $\epsilon_{\nu,\mathrm{react}}$,
         if the initial model that you are using in your simulation also included those losses.
         Otherwise, the energy loss from our NSE table will likely be too great and that simulation
         will not be in equilibrium.  This is an issue, for example, when using a MESA model
         constructed with ``aprox21``, which does not have all of the weak rates we model here.

         The weak rate neutrino losses can be disabled by ``integrator.nse_include_enu_weak=0``.

    * Predict $\Uc^\prime$ to the midpoint in time, $n+1/2$ and construct
      $[\Rb(\Uc^\prime)]^{n+1/2}$.

    * Do the final update to time $n$ as:

      .. math::

         \Uc^{\prime,n+1/2} = \Uc^{\prime,n} + \frac{\Delta t}{2} [\Advs{\Uc^\prime}]^{n+1/2} + \frac{\Delta t}{2} [\Rb(\Uc^\prime)]^{n+1/2}


      where $[\Advs{\Uc^\prime}]^{n+1/2}$ are the advective updates carried by the SDC
      algorithm.

    * Compute the energy generation rate from the change in internal energy from $\Uc^{\prime,n}$ to $\Uc^{\prime,n+1}$, excluding advection.

    * Update the total energy.

    * Set the mass fractions carried on the grid from the NSE table (with the new temperature and $Y_e$).

  * if we are not in NSE:

    * integrate the ``aprox19`` network as usual

    * update the aux quantities at the end of the burn


NSE check
---------

.. index:: network.rho_nse, network.T_nse, network.T_always_nse
.. index:: network.He_Fe_nse, network.C_nse, network.O_nse, network.Si_nse

For a zone to be consider in NSE, we require $\rho$ > ``network.rho_nse`` and *either*

* $T$ > ``network.T_nse`` together with the composition check

* $T$ > ``network.T_always_nse``

where we assume that ``T_always_nse`` > ``T_nse``.

The composition check considers the following nuclei groups:

* He-group: atomic numbers 1 to 2 (H to He)

* C-group: atomic numbers 6 to 7 (C to N)

* O-group: atomic number 8 (O)

* Si-group: atomic number 14 (Si)

* Fe-group: atomic numbers 24 to 30 (Cr to Zn)

and we then say that a composition supports NSE if:

* :math:`X(\mathrm{C}_\mathrm{group})` < ``network.C_nse``

* :math:`X(\mathrm{O}_\mathrm{group})` < ``network.O_nse``

* :math:`X(\mathrm{Si}_\mathrm{group})` < ``network.Si_nse``

* :math:`X(\mathrm{Fe}_\mathrm{group}) + X(\mathrm{He}_\mathrm{group})` > ``network.He_Fe_nse``



NSE table ranges
----------------

The NSE table was created for:

* :math:`9.4 < \log_{10}(T) < 10.4`
* :math:`7 < \log_{10}(\rho) < 10`
* :math:`0.43 < Y_e < 0.5`



.. _self_consistent_nse:

Self-consistent NSE
===================

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
-----------------

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
------------------

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
