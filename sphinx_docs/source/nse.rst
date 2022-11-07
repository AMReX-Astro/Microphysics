***
NSE
***

The reaction networks in Microphysics have the ability to use NSE
instead of integrating the entire network when the conditions are
appropriate.  There are 2 different implementations of NSE in
Microphysics, that have slightly different use cases.

* :ref:`tabulated_nse` : this uses a table of NSE abundances given
  :math:`(\rho, T, Y_e)` generate from a large network (125 isotopes).
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

.. prompt:: bash

   NETWORK_DIR=aprox19 USE_NSE=TRUE

Composition and EOS
-------------------

The NSE table was generated using a 125 nuclei reaction network
(described in :cite:`ma:2013`), and includes electron-capture rates,
so the compositional quantities it carries, :math:`\bar{A}` and
:math:`Y_e` and not representable from the 19 isotopes we carry in the
network.  In particular, it can attain a lower :math:`Y_e` than
``aprox19`` can represent.

For this reason, when we are using the NSE network, we always take the
composition quantities in the EOS directly from ``eos_state.aux[]``
instead of from ``eos_state.xn[]``.  The ``AUX_THERMO`` preprocessor
variable is enabled in this case, and the equations of state interpret
this to use the auxillary data for the composition.  This is described in :ref:`aux_eos_comp`.


NSE Table Outputs
-----------------

The NSE table provides values for the auxiliary composition,
:math:`Y_e`, :math:`\bar{A}`, and :math:`\langle B/A \rangle`
resulting from the full 125 nuclei network.   It also provides a set of 19
:math:`X_k` that map into the isotopes carried by ``aprox19``.


These three quantities are stored as ``aux`` data in the network and
are indexed as ``iye``, ``iabar``, and ``ibea``.  Additionally, when
coupling to hydrodynamics, we need to advect these auxillary
quantities.

For Strang split coupling of hydro and reactions, :math:`DX_k/Dt = 0`,
and our evolution equations are:

.. math::

   \begin{align*}
   \frac{DY_e}{Dt} &= \sum_k \frac{Z_k}{A_k} \frac{DX_k}{Dt} = 0 \\
   \frac{D}{Dt} \frac{1}{\bar{A}} &= - \frac{1}{\bar{A}^2} \frac{D\bar{A}}{Dt} = \sum_k \frac{1}{A_k} \frac{DX_k}{Dt} = 0 \rightarrow \frac{D\bar{A}}{Dt} = 0 \\
   \frac{D}{Dt} \left (\frac{B}{A} \right ) &= \sum_k \frac{B_k}{A_k} \frac{DX_k}{Dt} = 0
   \end{align*}

Therefore each of these auxillar equations obeys an advection equation
in the hydro part of the advancement.


NSE Flow
--------

The basic flow of a simulation using ``aprox19`` + the NSE table is as follows:

* initialize the problem, including :math:`X_k`

* fill the initial aux data with :math:`Y_e`, :math:`\bar{A}`, and :math:`(B/A)`

* in hydro, we will update these quantities simply via advection (for
  Strang-split evolution)

* for the reactive update:

  * check if NSE applies (see below)

  * if we are in an NSE region:

    * use :math:`\rho`, :math:`T`, and :math:`Y_e` to call the table.
      This returns: :math:`dY_e/dt`, :math:`(B/A)_{\rm out}`, and :math:`\bar{A}_{\rm out}`.

    * update :math:`Y_e` [#fY]_ :

      .. math::

         (Y_e)_{\rm out} = (Y_e)_{\rm in} + \Delta t \frac{dY_e}{dt}

    * :math:`\bar{A}_{\rm out}` is simply the value returned from the table

    * the energy generation rate, :math:`e_{\rm nuc}` is:

      .. math::

         e_{\rm nuc} = \eta \left [ \left ( \frac{B}{A} \right )_{\rm out} -
                                    \left ( \frac{B}{A} \right )_{\rm in} \right ] * \frac{1.602 \times 10^{-6}  {\rm erg}}{{\rm MeV}} N_A \frac{1}{\Delta t}


      where :math:`\eta` is an inertia term < 1 to prevent the energy changing too much in one set.

    * the new binding energy for the zone is then:

      .. math::

         \left ( \frac{B}{A} \right )_{\rm out}  = \left ( \frac{B}{A} \right )_{\rm in} + \eta \left [ \left ( \frac{B}{A} \right )_{\rm out} - \left ( \frac{B}{A} \right )_{\rm in} \right ]

    * update the mass fractions, :math:`X_k`, using the values from the table

  * if we are not in NSE:

    * integrate the ``aprox19`` network as usual

    * update the aux quantities at the end of the burn


NSE check
---------

We determine is a zone is in NSE according to:

* :math:`\rho` > ``rho_nse``

* :math:`T` > ``T_nse``

* :math:`X(\isotm{C}{12})` < ``C_nse``

* :math:`X(\isotm{He}{4}) + X(\isotm{Cr}{48}) + X(\isotm{Fe}{52}) + X(\isotm{Fe}{54}) + X(\isotm{Ni}{56})` > ``He_Fe_nse``

.. _self_consistent_nse:

Self-consistent NSE
===================

The self-consistent NSE approach uses only the nuclei in the main
reaction network.  It solves for the chemical potentials of the proton
and neutron and from there gets the abundances of each of the nuclei
under the assumption of NSE, following the procedure outlined in :cite:`Calder_2007`.

The solve is done using a port of the hybrid Powell method from
MINPACK (we ported the solver to templated C++).

The advantage of this approach is that it can be used with any
reaction network, once the integration has reached NSE.

This solver is enabled by compiling with

.. prompt:: bash

   USE_NSE_NET=TRUE

The functions to find the NSE state are then found in ``nse_solver.H``.



.. rubric:: Footnotes

.. [#fY] The table actually provides the weak rate, which is the sum
   of all electron capture and positron decay rates times the
   appropriate abundances minus a similar rate for the beta decay and
   positron capture, [wrate] = [rectot] + [rpdtot] - [redtot] - [rpctot]

   So if electron capture dominates, then [wrate] is positive and this should
   be subtracted from :math:`Y_e`.

