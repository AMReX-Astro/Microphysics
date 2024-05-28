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

The current list of tests are:

* ``burn_cell``

* ``burn_cell_primordial_chem``

* ``burn_cell_sdc``

* ``eos_cell``

* ``nse_table_cell``

* ``test_aprox_rates``

* ``test_ase``

* ``test_conductivity``

* ``test_eos``

* ``test_jac``

* ``test_linear_algebra``

* ``test_neutrino_cooling``

* ``test_nse``

* ``test_nse_interp``

* ``test_parameters``

* ``test_part_func``

* ``test_react``

* ``test_rhs``

* ``test_screening``

* ``test_screening_templated``

* ``test_sdc``

* ``test_sdc_vode_rhs``
