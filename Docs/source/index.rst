.. Microphysics documentation main file.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

************************
AMReX-Astro Microphysics
************************

AMReX-Astro Microphysics is a collection of microphysics routines (equations of state,
reaction networks, ...) and utilities (ODE integrators, NSE solvers)
for astrophysical simulation codes.

The original design was to support codes based on the `AMReX
<https://github.com/amrex-codes/amrex>`_ adaptive mesh refinement library :cite:`amrex_joss`,
including `CASTRO
<https://github.com/amrex-astro/Castro>`_ :cite:`castro_I`, `MAESTROeX
<https://github.com/amrex-astro/MAESTROeX>`_ :cite:`maestroex`, and Quokka :cite:`quokka`. These all have a
consistent interface and the separate Microphysics repository allows
them to share the same equation of state, reaction networks, and more.
Later, Microphysics was adopted by the `Quokka <https://github.com/quokka-astro/quokka>`_
simulation code.

While there are a number of unit tests that exercise the functionality,
Microphysics is primarily intended to be used along with another simulation
code.   At the moment, the interfaces and
build stubs are compatible with the AMReX codes and use the AMReX build
system.

A number of the routines contained here we authored by other people.
We bundle them here with permission, usually changing the interfaces
to be compatible with our standardized interface. We in particular
thank Frank Timmes for numerous reaction networks and his equation
of state routines.


.. toctree::
   :maxdepth: 1
   :caption: Microphysics overview
   :hidden:

   getting_started
   design
   data_structures
   autodiff
   rp_intro

.. toctree::
   :maxdepth: 1
   :caption: EOS and transport
   :hidden:

   eos
   eos_implementations
   transport

.. toctree::
   :maxdepth: 1
   :caption: Reaction networks
   :hidden:

   networks-overview
   networks
   templated_networks
   screening
   neutrinos

.. toctree::
   :maxdepth: 1
   :caption: ODE integrators
   :hidden:

   integrators
   ode_integrators
   nse
   sdc

.. toctree::
   :maxdepth: 1
   :caption: Util / external libraries
   :hidden:

   util

.. toctree::
   :maxdepth: 1
   :caption: Unit tests
   :hidden:

   unit_tests
   comprehensive_tests
   one_zone_tests

.. toctree::
   :maxdepth: 1
   :caption: API

   filelist
   classlist
.. namespacelist

.. toctree::
   :maxdepth: 1
   :caption: References
   :hidden:

   zreferences

.. toctree::
   :maxdepth: 1
   :caption: Index
   :hidden:

   genindex
