************
Build System
************

Microphysics leverages the AMReX build system.  GNU make is the
primary build system, but CMake is partially supported.  Here we focus
on the GNU make build system.

.. tip::

   All of the build variables supported by the AMReX build system can
   be used with Microphysics, for example:

   * ``COMP`` : the compiler system to use
   * ``USE_CUDA`` : build for NVIDIA GPUs with CUDA
   * ``USE_HIP`` : build for AMD GPUs with HIP/ROCm

   See the `AMReX GNU Make Docs
   <https://amrex-codes.github.io/amrex/docs_html/BuildingAMReX.html>`_

There are two main makefile stubs:

.. index:: Make.Microphysics, Make.Microphysics_extern

* ``Make.Microphysics`` : this is used by the unit tests within Microphysics and it
  defines the variables needed by the AMReX build system, including specifying the
  location of the files.  It also defines the rules for some files that are created
  at build time.

* ``Make.Microphysics_extern`` : this is the core makefile stub for Microphysics
  that interprets the build variables that enable / disable different functionality.
  This is meant to be included into any application code's build system.

additionally directories have their own ``Make.package`` files that specify
the files needed for the build and some rules for making intermediary files.


Environment variables
=====================

The build system relies on some environment variables to find the source:

.. index:: AMREX_HOME, MICROPHYSICS_HOME

* ``AMREX_HOME`` : this should point to the top-level ``amrex/`` directory

* ``MICROPHYSICS_HOME`` : this is needed by application codes, and
  should point to the top level ``Microphysics/`` directory.  For
  building unit tests within the ``Microphysics/`` directory itself,
  this does not need to be explicitly set.


Automatically generated files
=============================

There are a few source files that are created automatically at
compile-time.  These are placed in the
``tmp_build_dir/microphysics_sources/`` sub-directory under the
directory you run ``make`` (if building through an application code,
the sub-directory may have a different name).

.. index:: network_properties.H, extern_parameters.H, AMReX_buildInfo.cpp

The main files are:

* ``network_properties.H`` : this defines the properties of the composition that
  make up the network (and therefore, used by the equation of state and other
  physics).

* ``extern_parameters.H`` : this defines all of the runtime parameters that are
  part of the build.  At the moment, they are treated as global variables
  (using managed memory on GPUs), but a ``struct`` that carries their values
  is also available through ``extern_type.H``.

* ``AMReX_buildInfo.cpp`` : this defines functions that return the git hashes,
  compilers, compiler flags, and more meta-data about the build.  This file
  is automatically deleted once it is built to insure it is always up-to-date.
  The functions it defines are used when writing the ``job_info`` file
  in the plotfiles that some unit tests produce.

Controlling choice of physics
=============================

The choice of physics to include in an application is done at build time, and is
controlled by a number of make variables.

.. tip::

   You can query the value of any variable in the Microphysics build system by doing
   `make print-<NAME>` where `<NAME>` is the name of the variable.

   For example,

   .. code:: bash

      make print-EOS_DIR

   will tell you what EOS is being used.


The following control whether certain physics modules are included in
the build process.  Note: an EOS and network are always required.
These can be set to ``TRUE`` to enable and ``FALSE`` to disable specific features.

* ``USE_CONDUCTIVITY`` : determines whether a conductivity routine
  should be included in the list of build packages.  If enabled, this
  also defines the ``CONDUCTIVITY`` preprocessor variable.  Default:
  ``FALSE``.

* ``USE_NEUTRINOS`` : determines whether a neutrino cooling term
  should be applied in the reaction network energy generation
  equation.  See :ref:`neutrino_loss`.  The default is set by each
  individual network.

* ``USE_NET_NET`` : determines whether the self-consistent NSE
  infrastructure is included in the build.  See
  :ref:`self_consistent_nse`.  No default is set.

* ``USE_NSE_TABLE`` : determines whether the tabular NSE
  infrastructure is included in the build.  See :ref:`tabulated_nse`.
  No default is set.

* ``USE_RATES`` : for templated reaction networks (see
  :ref:`sec:templated_rhs`) determines whether we include the
  ``rates/`` set of reaction rates in the build system.  Also defines
  the ``RATES`` preprocessor variable.  The default is set by each of
  the templated networks separately.

* ``USE_REACT`` : determines whether we need to include any of the
  source related to reaction networks or integrators and sets the
  ``REACTIONS`` preprocessor variable.  Note: even if this is set to
  ``TRUE``, the ``network_properties.H`` file is still generated.  No
  default is set.

* ``USE_SCREENING`` : determines whether the screening routines are
  included in the list of build packages.  If enabled, this also
  defines the ``SCREENING`` preprocessor variable which is used in
  some networks to disable screening completely.  Note: it is also
  possible to set the screening routine to ``null`` which would have
  the same effect (see :ref:`sec:screening`).  The default is set by
  each individual network.


The following control the choice of implementation for the different physics modules:


* ``CONDUCTIVITY_DIR`` : the name of the conductivity implementation to use,
  relative to ``Microphysics/conductivity/``.

* ``EOS_DIR`` : the name of the EOS to use, relative to ``Microphysics/EOS/``.

* ``INTEGRATOR_DIR`` : the name of the integrator to use, relative to
  ``Microphysics/integration/``.

* ``NETWORK_DIR`` : the name of the network to use, relative to ``Microphysics/networks/``.
  If ``general_null`` is chosen, then the inputs file is determined by
  either ``GENERAL_NET_INPUTS`` or ``NETWORK_INPUTS`` (see :ref:`sec:networks:general_null`).

* ``OPACITY_DIR`` : the name of the opacity implementation to use, relative
  to ``Microphysics/opacity/``.

* ``SCREEN_METHOD`` : the name of the screening implementation to use.  The choices
  are listed in :ref:`sec:screening`.


The following control the time-integration method used by the reaction
network integration:

* ``USE_SIMPLIFIED_SDC`` : enable the simplified-SDC coupling of hydro and reactions.
  See :ref:`sec:simplified_sdc`.

* ``USE_TRUE_SDC`` : enable the true-SDC coupling of hydro and reactions.
  See :ref:`sec:true_sdc`.

.. note::

   If neither of these are set to ``TRUE``, then Strang-splitting coupling
   will be used.


Targets
=======

For the unit tests, simply doing

.. code:: bash

   make

in the test directory will build the test.  There are a few other targets defined.  The most important
is ``clean``.  Doing:

.. code:: bash

   make clean

will remove all the build temporary files (including the ``tmp_build_dir/``).

.. important::

   If you want to use a different EOS or reaction network (or any other physics), then you
   should always do ``make clean`` first in the build directory.

Some other targets include:

* ``nettables`` : create the symlinks for any weak reaction rate tables that are part of the
  network.

* ``table`` : create a symlink for the ``helm_table.dat`` EOS table if the ``helmholtz`` EOS is used.

* ``nsetable`` : create a symlink for the NSE table if ``USE_NSE_TABLE=TRUE`` is set.

* ``build_status`` : report the current git versions of Microphysics and AMReX

* ``test_extern_params`` : this will simply parse the runtime parameters and execute the
  ``write_probin.py`` script that generates the headers and C++ files necessary to use
  the parameters.  These will be generated under ``tmp_build_dir/microphysics_sources/``.

* ``net_prop_debug`` : this will simply create the ``network_properties.H`` file for the
  current network and output it into ``tmp_build_dir/microphysics_sources/``.
