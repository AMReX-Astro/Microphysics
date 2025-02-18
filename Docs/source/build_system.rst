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
  This is meant to be including into any application code's build system.

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

Controling choice of physics
============================

The choice of physics to include in an application is done at build time, and is
controlled by a number of make variables.

.. tip::

   You can query the value of any variable in the Microphysics build system by doing
   `make print-<NAME>` where `<NAME>` is the name of the variable.

   For example,

   .. code:: bash

      make print-EOS_DIR

   will tell you what EOS is being used.




USE_CONDUCTIVITY
USE_NEUTRINOS
USE_SCREENING
USE_RATES
USE_REACT
USE_NONAKA_PLOT

USE_SIMPLIFIED_SDC
USE_TRUE_SDC

SCREEN_METHOD

USE_AUX_THERMO

EOS_DIR
OPACITY_DIR
NETWORK_DIR
CONDUCTIVITY_DIR

GENERAL_NET_INPUTS
NETWORK_INPUTS

USE_NSE_TABLE
USE_NET_NET

USE_RAD

Targets
=======

nettables
table
nsetable
build_status
test_extern_params
net_prop_debug

clean


runtime parameters

buildInfo


Querying the build variables
============================
