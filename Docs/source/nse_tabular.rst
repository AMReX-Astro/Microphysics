.. _tabulated_nse:

*****************************
Tabulated NSE and ``aprox19``
*****************************

The ``aprox19`` network can be run in a manner where we blends the
standard ``aprox19`` network with a table for nuclear statistic
equilibrium resulting from a much larger network at high density and
temperatures.    This option is enabled by building with:

.. code:: bash

   NETWORK_DIR=aprox19 USE_NSE_TABLE=TRUE

Interface
=========

.. index:: nse_table_t, nse_interp

The main interface for the NSE table is defined in ``nse_table.H``:

.. code:: c++

   AMREX_GPU_HOST_DEVICE AMREX_INLINE
   void nse_interp(nse_table_t& nse_state, bool skip_X_fill=false)

Here, ``nse_table_t`` is a simple struct that holds the thermodynamic and
NSE state and the optional parameter ``skip_X_fill`` can be used to skip
interpolating the mass fractions if desired (since there are a lot of them,
and the interpolation can be expensive).

Composition and EOS
===================

The NSE table was generated using `pynucastro
<https://pynucastro.github.io/pynucastro/>`_ :cite:`pynucastro,
pynucastro2`, using 96 nuclei and electron/positron capture/decay
rates from :cite:`langanke:2001`.  The table takes $Y_e$ as the
primary composition variable and provides a set of mass fractions that
is mapped into those used by ``aprox19``.  Using the value allows us
to attain a lower :math:`Y_e` than ``aprox19`` can represent.

.. note::

   The full details of the NSE table are provided in :cite:`sdc-nse`.
   The table can be regenerated using the script ``nse_tabular/make_nse_table.py``.

When we are using the NSE network, we always take the
composition quantities in the EOS directly from ``eos_state.aux[]``
instead of from ``eos_state.xn[]``.  The ``AUX_THERMO`` preprocessor
variable is enabled in this case, and the equations of state interpret
this to use the auxiliary data for the composition.  This is described in :ref:`aux_eos_comp`.


NSE Table Outputs
=================

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
========

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
=========

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
================

The NSE table was created for:

* :math:`9.4 < \log_{10}(T) < 10.4`
* :math:`7 < \log_{10}(\rho) < 10`
* :math:`0.43 < Y_e < 0.5`
