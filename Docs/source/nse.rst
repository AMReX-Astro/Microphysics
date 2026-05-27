.. _ch:nse:

***
NSE
***

.. important::

   NSE is only supported with the simplified-SDC method for
   coupling hydrodynamics and reactions.  We do not support
   operator-splitting (Strang) coupling with NSE.

Nuclear statistical equilibrium (NSE) describes a state in which all strong
nuclear reactions are in chemical equilibrium, with each forward reaction
balanced by its reverse. The NSE condition is a state where the
chemical potential of any isotope can be expressed as
the sum of the chemical potentials of the its constituent neutrons and protons.
which happens when :math:`T \gtrsim 6 \times 10^9 \ \mathrm{K}`.
Mathematically it means:

.. math:: \mu_i = Z_i \mu_p + N_i \mu_n

Applying this assumption to the Maxwell-Boltzmann gas, the mass abundance
of the system in NSE is uniquely defined as:

.. math::

   X_i^{\mathrm{NSE}} = \dfrac{(m_i)^{5/2}}{\rho} g_i \left( \dfrac{k_B T}{2\pi \hbar^2} \right)^{3/2} \exp{\left( \dfrac{Z_i (\mu_p^\mathrm{kin} + \mu_p^c) + N_i \mu_n^{\mathrm{kin}} + B_i - \mu^c_i}{k_B T} \right)}

A more detailed derivation can be found in many literature,
e.g., :cite:`pynucastro2`.

There are three possible constraint equations that may be required to solve
the NSE equation. Depending on the specified input mode, either two or all three
constraints are needed.

1. *Conservation of mass*, or equivalently a constraint on the input density,
   :math:`\rho`.

   .. math:: f_\rho = \sum_i X_i^{\mathrm{NSE}} - 1 = 0

2. *Conservation of charge* for strong reactions, or equivalently a constraint
   on the input electron fraction, :math:`Y_e`

   .. math:: f_{Y_e}\sum_i \frac{Z_i}{A_i} X_i^{\mathrm{NSE}} - Y_e = 0

3. *Conservation of energy*, or equivalently a constraint on the
   input internal energy, :math:`e`

   .. math:: f_e = \frac{e\left(\rho, T_\mathrm{NSE}, X_i^\mathrm{NSE}\right)}{e_0} - 1 = 0

Two input modes are currently supported:

1. Given the input :math:`(\rho, T, Y_e)`, the system has two unknowns,
   the chemical potential of proton and neutron, i.e.
   :math:`(\mu_p^\mathrm{kin}, \mu_n^\mathrm{kin})`.
   In this mode, the first two constraint equations are sufficient to solve
   the system.

2. Given the input :math:`(\rho, e_0, Y_e)`, the above equation has three unknowns,
   the chemical potential of proton and neutron and temperature, i.e.
   :math:`(\mu_p^\mathrm{kin}, \mu_n^\mathrm{kin}, T_9)`.
   In this mode, all three constraint equations are needed to solve the system.


NSE Evolution Modes
===================

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
  species defined in the network.  It only works pynucastro-generated networks.
  Note that there are requirements in order for networks to be intrinsically
  compatible with NSE, see more detail in :ref:`self_consistent_nse`.
  Unlike the tabulated NSE, there is no need to advect an auxiliary composition,
  since this only deals with the isotopes defined in the main reaction network.

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

These two NSE solvers are in the next sections.

.. toctree::
   :maxdepth: 1
   :hidden:

   nse_tabular.rst
   nse_net.rst
