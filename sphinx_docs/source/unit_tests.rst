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
  interfaces, or runtime infrastructure.

  These are not really meant for exploring the actual thermodynamic
  state.



Comprehensive tests
===================

Each of these tests sets up a cube of data, $(\rho, T, X_k)$, with the
range of $T$ and $\rho$, and the species to focus on for $X_k$ controlled
by options in the input file.

* ``test_aprox_rates`` :

  call each of the hardcoded rate functions in ``Microphysics/rates/``
  on each cell in the data cube and store the output in a plotfile.

* ``test_conductivity`` :

  call one of the conductivity routines (set via ``CONDUCTIVITY_DIR``)
  on each cell in the data cube and store the output in a plotfile.

* ``test_eos`` :

  call one of the equations of state (set via ``EOS_DIR``) on each
  cell in the data cube. We first call it with $(\rho, T, X_k)$ as
  input (``eos_input_rt``), and then test each of the other EOS modes
  (``eos_input_rh``, ``eos_input_tp``, ``eos_input_rp``,
  ``eos_input_re``, ``eos_input_ps``, ``eos_input_ph``,
  ``eos_input_th``) and for each of these modes, we compute the error
  in the recovered $T$ and/or $\rho$ (as appropriate).  The full
  thermodynamic state and errors are stored in a plotfile.

* ``test_jac`` :

  for each cell in the data cube, and for a specific network (set via
  ``NETWORK_DIR``) call the analytic Jacobian provided by the network
  and compute the numerical Jacobian (via finite differencing) and
  store the relative difference between each approximation for each
  Jacobian element in a plotfile.

* ``test_neutrino_cooling`` :

  for each cell in the data cube, call the neutrino cooling function
  and store the output in a plotfile.

* ``test_react`` :

  for each cell in the data cube, call the reaction network and
  integrate to a specified time.  Statistics about the number of RHS
  calls are reported at the end.  A lot of options can be set via the
  inputs file to control the integration.

* ``test_rhs`` :

  for each cell in the data cube, call the network's righthand side and
  Jacobian functions and store the output in a plotfile.  The network
  is controlled by the ``NETWORK_DIR`` variable.

* ``test_screening`` :

* ``test_screening_templated``

* ``test_sdc`` :

  similar to ``test_react``, except now we exercise the SDC
  integration infrastructure.  The conserved state that is input to
  the burner is chosen to have a Mach number of $0.1$ (and otherwise
  use the thermodynamics from the data cube).  No advective terms are
  modeled.


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
