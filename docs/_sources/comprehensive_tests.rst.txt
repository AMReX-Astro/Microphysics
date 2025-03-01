.. _sec:comprehensive_tests:

************************
Comprehensive Unit Tests
************************

Generally, for each test, you simply type ``make`` in the test
directory.  There are a number of runtime parameters that can
control the behavior.  These are specified (along with defaults)
in ``_parameters`` files in each test directory and can be
overridden in an inputs file or on the commandline.

Some additional details on a few of the comprehensive unit tests
are given below.

EOS test (``test_eos``)
=======================

``Microphysics/unit_test/test_eos/`` is a unit test for the equations
of state in Microphysics. It sets up a cube of data, with
:math:`\rho`, :math:`T`, and :math:`X_k` varying along the three
dimensions, and then calls the EOS in each zone. Calls are done to
exercise all modes of calling the EOS, in order:

- ``eos_input_rt``: We call the EOS using :math:`\rho, T`. this is the
  reference call, and we save the :math:`h`, :math:`e`, :math:`p`, and
  :math:`s` from here to use in subsequent calls.

- ``eos_input_rh``: We call the EOS using :math:`\rho, h`, to recover
  the original :math:`T`. To give the root finder some work to do, we
  perturb the initial temperature.

  We store the relative error in :math:`T` in the output file.

- ``eos_input_tp``: We call the EOS using :math:`T, p`, to recover the
  original :math:`\rho`. To give the root finder some work to do, we
  perturb the initial density.

  We store the relative error in :math:`\rho` in the output file.

- ``eos_input_rp``: We call the EOS using :math:`\rho, p`, to recover
  the original :math:`T`. To give the root finder some work to do, we
  perturb the initial temperature.

  We store the relative error in :math:`T` in the output file.

- ``eos_input_re``: We call the EOS using :math:`\rho, e`, to recover
  the original :math:`T`. To give the root finder some work to do, we
  perturb the initial temperature.

  We store the relative error in :math:`T` in the output file.

- ``eos_input_ps``: We call the EOS using :math:`p, s`, to recover the
  original :math:`\rho, T`. To give the root finder some work to do,
  we perturb the initial density and temperature.

  Note: entropy is not well-defined for some EOSs, so we only attempt
  the root find if :math:`s > 0`.

  We store the relative error in :math:`\rho, T` in the output file.

- ``eos_input_ph``: We call the EOS using :math:`p, h`, to recover the
  original :math:`\rho, T`. To give the root finder some work to do,
  we perturb the initial density and temperature.

  We store the relative error in :math:`\rho, T` in the output file.

- ``eos_input_th``: We call the EOS using :math:`T, h`, to recover the
  original :math:`\rho`. To give the root finder some work to do, we
  perturb the initial density.

  Note: for some EOSs, :math:`h = h(\rho)` (e.g., an ideal gas), so there
  is no temperature dependence, and we do not do this test.

  We store the relative error in :math:`\rho` in the output file.

This unit test is marked up with OpenMP directives and therefore also
tests whether the EOS is threadsafe.

To compile for a specific EOS, e.g., helmholtz, do::

    make EOS_DIR=helmholtz -j 4

Examining the output (an AMReX plotfile) will show you how big the
errors are. You can use the ``amrex/Tools/Plotfile/`` tool
``fextrema`` to display the maximum error for each variable.


Network test (``test_react``)
=============================

``Microphysics/unit_test/test_react/`` is a unit test for the nuclear
reaction networks in Microphysics. It sets up a cube of data, with
:math:`\rho`, :math:`T`, and :math:`X_k` varying along the three
dimensions (as a :math:`16^3` domain), and then calls the EOS in each
zone.  This test does the entire ODE integration of the network for
each zone.

The state in each zone of our data cube is determined by the runtime
parameters ``dens_min``, ``dens_max``, ``temp_min``, and ``temp_max``
for :math:`(\rho, T)`. Because each network carries different
compositions, we specify the composition through runtime parameters in
the ``&extern`` namelist: ``primary_species_1``,
``primary_species_2``, ``primary_species_3``. These primary species
will vary from X = 0.2 to X = 0.7 to 0.9 (depending on the number).
Only one primary species varies at a time. The non-primary species
will be set equally to share whatever fraction of 1 is not accounted
for by the primary species mass fractions.

This test calls the network on each zone, running for a time
``tmax``. The full state, including new mass fractions and energy
release is output to a AMReX plotfile.

You can compile for a specific integrator (e.g., ``VODE``) or
network (e.g., ``aprox13``) as::

    make NETWORK_DIR=aprox13 INTEGRATOR_DIR=VODE -j 4

The loop over the burner is marked up for OpenMP and CUDA and
therefore this test can be used to assess threadsafety of the burners
as well as to optimize the GPU performance of the burners.
