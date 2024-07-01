**********************
Overview of Unit Tests
**********************

There are a few unit tests in Microphysics that operate on a generic
EOS, reaction network, conductivity, or some smaller component to
Microphysics.  Many of these tests exercise the main interfaces in
``Microphysics/interfaces/`` and the code that those call.

These tests compile using the AMReX build system, which assumes that
main is in C++, so each have a ``main.cpp`` driver.  The files
``Microphysics/Make.Microphysics`` and
``Microphysics/Make.Microphysics_extern`` provide the macros necessary
to build the executable. Runtime parameters are parsed in the same
fashion as in the application codes, using the ``write_probin.py``
script.

.. note::

   Most of these tests work with MPI+OpenMP, MPI+CUDA, and MPI+HIP

Tests are divided into three categories:

* *comprehensive tests* work on a cube of data (usually
  $\rho$, $T$, and composition varying along the three dimensions) and
  are meant to exercise a wide range of input conditions.

  These are mainly used for regression testing.

* *one-zone tests* allow you to evaluate the conditions for a
  particular thermodynamic state.

  These are often used for interactive explorations and within the CI.

* *infrastructure tests* test small bits of the solver, function
  interfaces, or runtime infrastructure.  These are not really meant for
  exploring the actual thermodynamic state.



Comprehensive tests
===================

* ``test_aprox_rates``

* ``test_conductivity``

* ``test_eos``

* ``test_jac``

* ``test_neutrino_cooling``

* ``test_react``

* ``test_rhs``

* ``test_screening``

* ``test_screening_templated``

* ``test_sdc``



One-zone tests
==============

* ``burn_cell``

* ``burn_cell_primordial_chem``

* ``burn_cell_sdc``

* ``eos_cell``

* ``nse_table_cell``

* ``test_ase``

* ``test_nse``

* ``test_part_func``


Infrastructure tests
====================

* ``test_linear_algebra``

* ``test_nse_interp``

* ``test_parameters``

* ``test_sdc_vode_rhs``
