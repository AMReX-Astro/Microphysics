***************
Getting Started
***************

Requirements
============

Microphysics requires

* A C++17 or later compilers
* AMReX (https://github.com/amrex-codes/amrex)
* python (>= 3.10)
* GNU make

optional dependencies are:

* CUDA (>= 11)
* ROCm (>= 6.3.1 --- earlier versions have register allocation bugs)

Microphysics is meant to be compiled into an application code as part
of its build process, with the network, EOS, integrator, and more
picked at compile time.  As such, there is not a single library that
can be built and linked against.

Below we describe how to use Microphysics in a "standalone" fashion,
using the provided unit tests, and as part of an application code.

Standalone
==========

.. index:: AMREX_HOME

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

.. index:: burn_cell

A good unit test to start with is ``burn_cell`` --- this is simply a
one-zone burn.  In ``Microphysics/`` do:

.. prompt:: bash

   cd unit_test/burn_cell
   make

This will create an executable called ``main3d.gnu.ex``.
By default, the test is built with the 13-isotope ``aprox13`` network,
``helmholtz`` EOS, and VODE integrator.


Then you can run it as:

.. prompt:: bash

   ./main3d.gnu.ex inputs_aprox13

Here ``inputs_aprox13`` is the inputs file that sets options.

This will output information about the starting and final state to the
terminal and produce a file ``state_over_time.txt`` that contains the
thermodynamic history at different points in time.

.. note::

   See the :ref:`sec:burn_cell` documentation for more details on this
   unit test and how to visualize the output.

Running with AMReX Application Code
===================================

.. index:: MICROPHYSICS_HOME

Getting started with Microphysics using either `CASTRO
<https://amrex-astro.github.io/Castro/docs/index.html>`_ or `MAESTROeX
<https://amrex-astro.github.io/MAESTROeX/docs/index.html>`_ is
straightforward. Because the modules here are already in a format that
the AMReX codes understand, you only need to provide to the code
calling these routines their location on your system. The code will do
the rest.

First we need to define the ``MICROPHYSICS_HOME`` environment
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

Other AMReX-based codes can use Microphysics in the same fashion.
