***
NSE
***

The ``aprox19`` can be run in a manner where we blends the ``aprox19``
network with a table for nuclear statistic equilibrium at high density
and temperatures.  This is based on the table described in
:ref:`ma:2013`.  This option is enabled by building with ``USE_NSE=TRUE``.

The NSE table provides:

.. math::

   \begin{align*}
   Y_e &= \sum_k \frac{Z_k X_k}{A_k} \\
   \bar{A} &= \left [ \sum_k \frac{X_k}{A_k} \right ]^{-1} \\
   \frac{B}{A} &= \sum_k \frac{B_k X_k}{A_k}
   \end{align*}

where :math:`B_k` is the binding energy of nucleus :math:`k`.

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

Composition and EOS
===================

The NSE table was generated using a 125 nuclei reaction network, so
the compositional quantities it carries, :math:`\bar{A}` and
:math:`Y_e` and not representable from the 19 isotopes we carry in the
network.  For this reason, when we are using the NSE network, we
always take the composition quantities in the EOS directly from
``eos_state.aux[]`` instead of from ``eos_state.xn[]``.  The
``NSE_THERMO`` preprocessor directive is enabled in this case, and the
equations of state interpret this to use the auxillary data for the
composition.

The equation of state also needs :math:`\bar{Z}` which is easily computed as

.. math::

   \bar{Z} = \bar{A} Y_e

NSE Flow
========

The basic flow of a simulation using the NSE network is as follows:

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
=========

We determine is a zone is in NSE according to:

* :math:`\rho` > ``rho_nse``

* :math:`T` > ``T_nse``

* :math:`X(\isotm{C}{12})` < ``C_nse``

* :math:`X(\isotm{He}{4}) + X(\isotm{Cr}{48}) + X(\isotm{Fe}{52}) + X(\isotm{Fe}{54}) + X(\isotm{Ni}{56})` > ``He_Fe_nse``


.. rubric:: Footnotes

.. [#fY] The table actually provides the weak rate, which is the sum
   of all electron capture and positron decay rates times the
   appropriate abundances minus a similar rate for the beta decay and
   positron capture, [wrate] = [rectot] + [rpdtot] - [redtot] - [rpctot]

   So if electron capture dominates, then [wrate] is positive and this should
   be subtracted from :math:`Y_e`.

