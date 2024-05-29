***************
Getting Started
***************

Standalone
==========

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



Running with AMReX Application Code
===================================

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
