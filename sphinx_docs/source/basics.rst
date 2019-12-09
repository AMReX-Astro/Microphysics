*****************
StarKiller Basics
*****************

Getting Started
===============

Getting started with Microphysics using either CASTRO or MAESTRO is
straightforward. Because the modules here are already in a format that
the AMReX codes understand, you only need to provide to the code
calling these routines their location on your system. The code will do
the rest. To do so, define the ``MICROPHYSICS_HOME`` environment
variable, either at a command line or (if you use the bash shell)
through your ``~/.bashrc``, e.g.::

 export MICROPHYSICS_HOME=/path/to/Microphysics

For CASTRO  the name of the EOS and network are set via the make
variables ``EOS_DIR`` and ``NETWORK_DIR``. The macros in CASTRO’s
``Make.Castro`` will know to look in Microphysics using the
``MICROPHYSICS_HOME`` variable to find the codes.

For MAESTRO, the name of the EOS and network are set via the make
variables ``EOS_DIR`` and ``NETWORK_DIR``, and the macros in MAESTRO’s
``GMaestro.mak`` file will find the code, again using the
``MICROPHYSICS_HOME`` variable.

For other codes, one can use the interfaces in
``Microphysics/interfaces/`` and sample routines in
``Microphysics/unit_test/`` to incorporate these modules into your
code. Note that there are a few AMReX files required at the moment
(mainly error handling and constants).

Structure
=========

The high-level directory structure delineates the types of microphysics
and the generic solvers:

* ``docs/``: this User’s Guide

* ``EOS/``: the various equations of state

* ``integration/``: the ODE integration routines used for the
  reaction networks

* ``interfaces/``: copies of the main derived types that are used to
  interface with the EOS and reaction networks. Note that most application
  codes will have their own local copies. These are provided for unit testing
  in Microphysics.

* ``networks/``: the nuclear reaction networks. This is mostly just the
  righthand side of the network, as the actual integrators are decoupled from
  the network.

* ``neutrinos/``: neutino loss source terms for the network energy equation.

* ``rates/``: common nuclear reaction rate modules used by some of the
  networks.

* ``screening/``: common electron screening factors used by some of the
  reaction networks.

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
finding for the EOS) is separated from the specific implmentations of
the microphysics.

**All quantities are assumed to be in CGS units**, unless otherwise
specified.
