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

These two NSE solvers are in the next sections.

.. toctree::
   :maxdepth: 1
   :hidden:

   nse_tabular.rst
   nse_net.rst
