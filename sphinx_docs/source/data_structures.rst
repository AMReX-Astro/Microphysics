.. _data_structures:

***************
Data Structures
***************

All of the routines in this software package are standardized so that
you interact with them using the same type of data structure, a C++ struct.

EOS
===

``eos_t``
---------

The main data structure for interacting with the EOS is ``eos_t``.
This is a collection of data specifying the microphysical state of the
fluid that we are evaluating. This has many components. For a
particular instantiation named ``eos_state``, the most important
data is the following:

* ``eos_state.rho`` : density [:math:`\mathrm{g~cm^{-3}}`]

* ``eos_state.T`` : temperature [K]

* ``eos_state.p`` : pressure [:math:`\mathrm{erg~cm^{-3}}`]

* ``eos_state.e`` : specific internal energy [:math:`\mathrm{erg~g^{-1}}`]

* ``eos_state.h`` : specific enthalpy [:math:`\mathrm{erg~g^{-1}}`]

* ``eos_state.s`` : specific entropy [:math:`\mathrm{erg~g^{-1}~K^{-1}}`]

* ``eos_state.xn[]`` : mass fractions of species (this is an array, dimensioned to be the number of species, ``NumSpec`` )

* ``eos_state.aux[]`` : any auxiliary variables carried with the fluid (this is an array, dimensioned to be the number of auxiliary quantities, ``NumAux`` )

Note that both ``NumSpec`` and ``NumAux`` are meant to be properties of the
network, and they will come in through the ``network_properties.H`` header file.

There is a lot more information that can be saved here, such as the
partial derivatives of the thermodynamic state variables with respect
to each other. To see a complete list, examine the ``eos_type.H``
file: ``Castro/Microphysics/interfaces/eos_type.H``.

Networks
========

``burn_t``
----------

The main data structure for interacting with the reaction networks is
``burn_t``. This holds the composition (mass fractions), thermodynamic
state, and a lot of internal information used by the reaction network
(e.g. the righthand side of the ODEs, the Jacobian, etc.). Typically
the user will only need to fill/use the following information:

* ``burn_state.rho``: density [:math:`\mathrm{g~cm^{-3}}`]

* ``burn_state.T``: temperature [K]

* ``burn_state.e``: the specific internal energy [:math:`\mathrm{erg~g^{-1}}`]

   Note: this has two different contexts, depending on when it is
   accessed.

   When you call the integrator and are in the process of integrating
   the reaction system, e will be an integration variable and
   will account for the nuclear energy release.  It will also be used to
   derive the temperature via the EOS.

   Upon exit of the integration, the initial internal energy (offset)
   is subtracted off, and e now represents the specific nuclear
   energy release from the reactions.

* ``burn_state.xn[]``: the mass fractions

* ``burn_state.aux[]``: any auxiliary quantities (like :math:`Y_e`)

* ``burn_state.i``, ``.j``, ``.k``: hydrodynamic zone i, j, k for bug reporting, diagnostics

* ``burn_state.time``: the time since the start of the integration [s]

   Note this is not the same as the simulation time. Each integrator
   will also store the simulation time at the start of integration
   in their local storage—this can be used as an offset to convert
   between integration and simulation time.

``rate_t``, ``rate_fr_t``
-------------------------

The ``rate_t`` and ``rate_fr_t`` structures are used internally in a network to pass the
raw reaction rate information (usually just the temperature-dependent
terms) between various subroutines. It does not come out of the
network-specific righthand side or Jacobian routines.

You can see their definitions in ``networks/rate_type.H``.

``burn_type.H``
---------------

In addition to defining the ``burn_t`` type, the header ``burn_type.H``
also defines integer indices into the solution vector that can be used
to access the different components of the state:

* ``neqs`` : the total number of variables we are integrating.

   It is assumed that the first ``nspec`` are the species.

* ``net_ienuc`` : the index of the specific internal energy in the solution vector

Integrators
===========

Each integrator also has their own internal data structure that holds
the information needed for the integration.  Meta-data that is not
part of the integration vector of ODEs, but is attached to a
particular state (:math:`X_k`, :math:`T`, :math:`e`), is stored in the
``burn_t`` and can be passed into the righthand side routine.

Converting Between Types
========================

There is significant overlap between ``eos_t`` and ``burn_t``.
The ``burn_type.H`` header two routines,
``burn_to_eos`` and ``eos_to_burn`` that convert a ``burn_t``
state to an ``eos_t`` state, and back. Only the thermodynamic
variables that are common in the two types are copied. This is
useful, for example, if you have a burn_t state and what to get
thermodynamic information by calling the EOS.

