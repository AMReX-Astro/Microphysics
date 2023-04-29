*****************
StarKiller Basics
*****************

Getting Started (Standalone)
============================

Microphysics can be used in a "standalone" fashion to run the unit
tests and explore the behavior of the reaction networks.  The main
requirement is a copy of AMReX:

.. prompt:: bash

   git clone https://github.com/AMReX-Codes/amrex.git

We use this for some data structures and the build system.  You need
to set the ``AMREX_HOME`` environment variable to point to the
``amrex/`` directory:

.. prompt:: bash

   export AMREX_HOME=/path/to/amrex

(where you change ``/path/to/amrex`` to your actual path).

A good unit test to start with is ``burn_cell`` -- this is simply a
one-zone burn.  In ``Microphysics/`` do:

.. prompt:: bash

   cd unit_test/burn_cell
   make

This will create an executable called ``main3d.gnu.ex``.  Then you can run it as:

.. prompt:: bash

   ./main3d.gnu.ex inputs_aprox21

By default, the test is built with the 21-isotope ``aprox21`` network.
Here ``inputs_aprox21`` is the inputs file that sets options. 



Getting Started (Running with MAESTROeX or CASTRO)
==================================================

Getting started with Microphysics using either CASTRO or MAESTROeX is
straightforward. Because the modules here are already in a format that
the AMReX codes understand, you only need to provide to the code
calling these routines their location on your system. The code will do
the rest. To do so, define the ``MICROPHYSICS_HOME`` environment
variable, either at a command line or (if you use the bash shell)
through your ``~/.bashrc``, e.g.:

.. code:: bash

   export MICROPHYSICS_HOME=/path/to/Microphysics

For CASTRO and MAESTROeX the name of the EOS and network are set via
the make variables ``EOS_DIR`` and ``NETWORK_DIR``. These codes then
rely on the Microphysics ``Make.Microphysics_extern`` makefile stub
(found via the ``MICROPHYSICS_HOME`` variable) to add the necessary
source to the build.  All of the interfaces that these codes use
are found in ``Microphysics/interfaces/``.

Other codes can use Microphysics in the same fashion.  Unit tests in
``Microphysics/unit_test/`` provide some examples of using the
interfaces.

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

**All quantities are assumed to be in CGS units**, unless otherwise
specified.
