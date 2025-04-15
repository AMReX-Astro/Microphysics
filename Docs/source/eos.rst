******************
Equations of State
******************

The general interface to the equation of state is:

.. code:: c++

   template <typename I, typename T>
   AMREX_GPU_HOST_DEVICE AMREX_INLINE
   void eos (const I input, T& state, bool use_raw_inputs = false)

where ``input`` specifies which thermodynamic quantities are taken as
the input and ``state`` is a C++ struct that holds all of the
thermodynamic information.


Interface and modes
===================

.. index:: eos_t, eos_re_t, eos_rep_t, eos_rh_t, chem_eos_t

The EOS is called as:

.. code:: c++

   eos(mode, eos_type)

where *mode* determines which thermodynamic quantities are inputs,
and is one of:

* ``eos_input_rt`` : density and temperature are inputs

* ``eos_input_rh`` : density and specific enthalpy are inputs

* ``eos_input_tp`` : temperature and pressure are inputs

* ``eos_input_rp`` : density and pressure are inputs

* ``eos_input_re`` : density and specific internal energy are inputs

* ``eos_input_ps`` : pressure and entropy are inputs

* ``eos_input_ph`` : pressure and specific enthalpy are inputs

* ``eos_input_th`` : temperature and specific enthalpy are inputs

The *eos_type* passed in is one of

* ``eos_t`` : provides access to all available thermodynamic information,
  including derivatives.

* ``eos_re_t`` : only provides the energy-based thermodynamic information, including
  energy derivatives.

* ``eos_rep_t`` : expands on ``eos_re_t`` to include pressure information

* ``eos_rh_t`` : expands on ``eos_rep_t`` to include enthalpy information

* ``chem_eos_t`` : adds some quantities needed for the primordial chemistry EOS
  and explicitly does not include the mass fractions.

.. tip::

   The EOS implementations make heavy use of templating to
   "compile-out" the thermodynamic quantities that are not needed
   (depending on the input type).  This can greatly increase
   performance.  As such, you should pick the smallest EOS structure
   (``eos_re_t``, ``eos_rep_t``, ...) that contains the thermodynamic
   information you need.

.. tip::

   You can also pass a ``burn_t`` struct into the EOS, although this
   will give you access to a much smaller range of data.


Composition
===========

All input modes for ``eos()`` require a composition.  Usually this is
via the set of mass fractions, ``eos_t xn[]``, but if ``USE_AUX_THERMO``
is set to ``TRUE``, then we instead use the auxiliary quantities
stored in ``eos_t.aux[]``.

.. _aux_eos_comp:

Auxiliary composition
---------------------

.. index:: USE_AUX_THERMO

.. note::

   The main use-case for the auxiliary composition is when using a reaction
   network together with the tabulated NSE state.

With ``USE_AUX_THERMO=TRUE``, we interpret the composition from the auxiliary variables.
For ``eos_t eos_state``, the auxiliary variables are


* ``eos_state.aux[AuxZero::iye]`` : electron fraction, defined as

  .. math::

     Y_e = \sum_k \frac{X_k Z_k}{A_k}

* ``eos_state.aux[AuxZero::iabar]`` : the average mass of the nuclei, :math:`\bar{A}`, defined as:

  .. math::

     \frac{1}{\bar{A}} = \sum_k \frac{X_k}{A_k}

  In many stellar evolutions texts, this would be written as :math:`\mu_I`.

* ``eos_state.aux[AuxZero::ibea]`` : the binding energy per nucleon (units of
  MeV), defined as

  .. math::

     \left \langle \frac{B}{A} \right \rangle  = \sum_k \frac{X_k B_k}{A_k}

  where :math:`B_k` is the binding energy of nucleus :math:`k`

Given a composition of mass fractions, the function:

.. code:: c++

   template <class state_t>
   AMREX_GPU_HOST_DEVICE AMREX_INLINE
   void set_aux_comp_from_X(state_t& state)

will initialize the auxiliary data.

Many equations of state also need :math:`\bar{Z}` which is easily computed as

.. math::

   \bar{Z} = \bar{A} Y_e


Composition derivatives
-----------------------

.. index:: eos_extra_t, eos_re_extra_t, eos_rep_extra_t

The derivatives $\partial p/\partial A$, $\partial p/\partial Z$,
and $\partial e/\partial A$, $\partial e/\partial Z$ are available via
the ``eos_extra_t``, ``eos_re_extra_t``, ``eos_rep_extra_t``, which
extends the non-"extra" variants with these additional fields.

The composition derivatives can be used via the ``composition_derivatives()`` function
in ``eos_composition.H``
to compute :math:`\partial p/\partial X_k |_{\rho, T, X_j}`, :math:`\partial e/\partial X_k |_{\rho, T, X_j}`, and :math:`\partial h/\partial X_k |_{\rho, T, X_j}`.  This has the interface:

.. code:: c++

   template <typename T>
   AMREX_GPU_HOST_DEVICE AMREX_INLINE
   eos_xderivs_t composition_derivatives (const T& state)



Initialization and cutoff values
================================


The EOS will make sure that the inputs are within an acceptable range,
(e.g., ``small_temp`` :math:`< T <` ``maxT``). If they are not, then it
resets them silently—no error is thrown.

If you are calling the EOS with ``eos_input_re``, and if :math:`e <
10^{-200}`, then it calls the EOS with ``eos_input_rt`` with T =
max ( ``small_temp``, T ).

.. note::

   User’s are encourage to do their own validation of inputs before calling
   the EOS.

EOS structure
=============

Each EOS should have two main routines through which it interfaces to the
rest of Microphysics.

* ``actual_eos_init()`` :  At the beginning of the simulation,
  ``actual_eos_init`` will perform any initialization steps and save
  EOS variables (mainly ``smallt``, the temperature floor, and
  ``smalld``, the density floor). These variables are stored in the
  main EOS module of the code calling these routines.

  This is also where an EOS with tables would read in the tables
  and initialize the memory they are stored in.

* ``actual_eos()`` : this is the main evaluation routine.  It should
  accept an ``eos_input_t`` specifying the thermodynamic inputs and a
  struct (like ``eos_t``) that stores the thermodynamic quantities.
