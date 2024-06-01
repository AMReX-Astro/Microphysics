******
Design
******

Structure
=========

The high-level directory structure delineates the types of microphysics
and the generic solvers:

* ``conductivity``: thermal conductivity routines

* ``constants``: fundamental constants

* ``EOS/``: the various equations of state

* ``integration/``: the ODE integration routines used for the
  reaction networks

* ``interfaces/``: the main structs / derived types that are used to
  interface with the EOS and reaction networks.

* ``networks/``: the nuclear reaction networks. This is mostly just the
  righthand side of the network, as the actual integrators are decoupled from
  the network.

* ``neutrinos/``: neutino loss source terms for the network energy equation.

* ``opacity/``: opacity routines for radiation transport

* ``rates/``: common nuclear reaction rate modules used by some of the
  networks.

* ``screening/``: common electron screening factors used by some of the
  reaction networks.

* ``sphinx_docs``: the sphinx source for this documentation

* ``unit_test/``: self-contained unit tests for Microphysics. These don’t
  need any application code to build, but will require AMReX.

* ``util/``: linear algebra solvers and other routines.

Design Philosophy
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
finding for the EOS) is separated from the specific implementations of
the microphysics.

.. note::

   All quantities are assumed to be in CGS units, unless otherwise
   specified.
