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

* ``Make.Microphysics``

* ``Make.Microphysics_extern``

additionally directories have their own ``Make.package`` files that specify
the files needed for the build.


AMREX_HOME
MICROPHYSICS_HOME


Controling choice of physics
============================

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
