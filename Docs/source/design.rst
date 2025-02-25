******
Design
******

Structure
=========

The high-level directory structure delineates the types of microphysics
and the generic solvers:

* ``conductivity/``: thermal conductivity routines

* ``constants/``: fundamental constants

* ``Docs/``: the sphinx source for this documentation

* ``EOS/``: the various equations of state

* ``integration/``: the ODE integration routines used for the
  reaction networks

* ``interfaces/``: the main structs / derived types that are used to
  interface with the EOS and reaction networks.

* ``networks/``: the nuclear reaction networks. This is mostly just the
  righthand side of the network, as the actual integrators are decoupled from
  the network.

* ``neutrinos/``: neutrino loss source terms for the network energy equation.

* ``opacity/``: opacity routines for radiation transport

* ``rates/``: common nuclear reaction rate modules used by some of the
  networks.

* ``screening/``: common electron screening factors used by some of the
  reaction networks.

* ``unit_test/``: self-contained unit tests for Microphysics. These don’t
  need any application code to build, but will require AMReX.

* ``util/``: linear algebra solvers and other routines.


.. note::

   All quantities are assumed to be in CGS units, unless otherwise
   specified.

Design philosophy
=================

Any application that uses Microphysics will at minimum need to
choose an EOS and a network. These two components work together. The
design philosophy is that the EOS depends on the network, but not the
other way around. The decision was made for the network to act as the
core module, and lots of code depends on it. This avoids circular
dependencies by having the main EOS datatype, ``eos_t``, and the
main reaction network datatype, ``burn_t``, be built on top of the
network.

The network is meant to store the properties of the species (typically
nuclear isotopes) including their atomic weights and numbers, and also
describes any links between the species when burning.

The equation of state relates the thermodynamic properties of the
material. It depends on the composition of the material, typically
specified via mass fractions of the species, and uses the properties
of the species defined by the network to interpret the state.

We try to maximize code reuse in the Microphysics source, so the
solvers (ODE integration for the network and Newton-Raphson root
finding for the EOS) are separated from the specific implementations of
the microphysics.



GPU considerations
==================

.. index:: GPUs

All of the Microphysics routines are written to run on GPUs.  This is
enabled in application codes by using the AMReX lambda-capturing
mechanism (see the `AMReX GPU documentation <https://amrex-codes.github.io/amrex/docs_html/GPU.html>`_
for more information).

This means leveraging the AMReX data-structures, macros, and
functions.  The unit tests (see :ref:`sec:unit_tests`) provide a good
reference for how to interface the Microphysics solvers and physics
terms with an AMReX-based code.

There are a few places where Microphysics behaves slightly differently
when running on a CPU vs. a GPU:

* In the VODE integrator, we disable Jacobian-caching to save memory.
  See :ref:`ch:networks:integrators`.

* In general we disable printing from GPU kernels, due to register
  pressure.  Some output can be enabled by compiling with
  ``USE_GPU_PRINTF=TRUE``.
