**********
Unit Tests
**********


Comprehensive Unit Tests
========================

There are a few unit tests in Microphysics that operate on a generic EOS
or reaction network. To allow these to compile independent of any
application code (e.g., Maestro or Castro), copies of the EOS
driver ``eos.f90`` and network interface ``network.f90`` are
provided in ``Microphysics/unit_test/``. 

These tests compile using the AMReX build system, which assumes that
main is in C++, so each have a ``main.cpp`` driver.  The file
``Microphysics/Make.Microphysics`` provides the macros necessary to build
the executable. Runtime parameters are parsed in the same fashion as
in the application codes, using the ``write_probin.py`` script.

EOS test
--------

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
errors are. You can use the ``AMReX/Tools/Postprocessing/`` tool
``fextrema`` to display the maximum error for each variable.


Network test
------------

``Microphysics/unit_test/test_react/`` is a unit test for the nuclear
reaction networks in Microphysics. It sets up a cube of data, with
:math:`\rho`, :math:`T`, and :math:`X_k` varying along the three
dimensions (as a :math:`16^3` domain), and then calls the EOS in each
zone.

The state in each zone of our data cube is determined by the runtime
parameters ``dens_min``, ``dens_max``, ``temp_min``, and
``temp_max`` for :math:`(\rho, T)`. Because each network carries different
compositions, we explicitly specify the mass fraction of each species
for every step in the composition dimension in a file ``xin.*`` for
each network. Note: it is up to the user to ensure that they species
are in the proper order in that file and sum to 1. The name of the
file to use is specified by the runtime parameter ``xin_file``.

This test calls the network on each zone, running for a time
tmax. The full state, including new mass fractions and energy
release is output to a AMReX plotfile.

You can compile for a specific integrator (e.g., ``BS``) or
network (e.g., ``aprox13``) as::

    make NETWORK_DIR=aprox13 INTEGRATOR_DIR=BS -j 4

The loop over the burner is marked up for both OpenMP and OpenACC and
therefore this test can be used to assess threadsafety of the burners
as well as to optimize the GPU performance of the burners.

Individual Network Tests
========================

Many of the networks have a subdirectory test/ (e.g.
``Microphysics/networks/triple_alpha_plus_cago/test/``). There are
usually 3 different drivers there that can be used to evaluate the
network or Jacobian on a single state:

- ``eval.f90``

- ``testburn.f90``

- ``testjacobian.f90``

These all use the F90 AMReX build system and the macros in
``GMicrophysics.mak`` to build the executable. To make
individual tests you can use the programs variable, e.g.,::

    make programs=eval



``burn_cell``
=============

``burn_cell`` is a simple one-zone burn that will evolve a state with
a network for a specified amount of time.  This can be used to
understand the timescales involved in a reaction sequence or to
determine the needed ODE tolerances.


Getting Started
---------------

The ``burn_cell`` code are located in
``Microphysics/unit_test/burn_cell``. To run a simulation, ensure that
both an input file and an initial conditions file have been created
and are in the same directory as the executable.

Input File
----------

These files are typically named as ``inputs_burn_network`` where network
is the network you wish to use for your testing.

The structure of this file is is fairly self-explanatory.  The run
prefix defined should be unique to the tests that will be run as they
will be used to identify all of the output files. Typically, the run
prefix involves the name of the network being tested.  The ``atol``
variables define absolute tolerances of the ordinary differential
equations and the ``rtol`` variables define the relative tolerances.  The
second section of the input file collects the inputs that ``main.f90``
asks for so that the user does not have to input all 5+
parameters that are required everytime the test is run.  Each input
required is defined and initialized on the lines following
``&cellparams``.  The use of the parameters is show below:

.. table:: The definition of parameters used in the burn_cell unit tests and specified in the second half of each inputs file.

   +-----------------------+----------------------------------------+
   | ``tmax``              | Maximum Time (s)                       |
   +-----------------------+----------------------------------------+
   | ``numsteps``          | Number of time subdivisions            |
   +-----------------------+----------------------------------------+
   | ``density``           | State Density (:math:`\frac{g}{cm^3}`) |
   +-----------------------+----------------------------------------+
   | ``temperature``       | State Temperature (K)                  |
   +-----------------------+----------------------------------------+
   | ``massfractions(i)``  | Mass Fraction for element i            |
   +-----------------------+----------------------------------------+

Running the Code
----------------

To run the code, enter the burn_cell directory and run::

   ./main.Linux.gfortran.exe with inputs

where ``inputs`` is the name of your inputs file.

For each of the ``numsteps`` steps defined in the inputs
file, the code will output a files into a new directory titled
``run_prefix_output`` where ``run_prefix`` is the run prefix defined in the
inputs file.  Each output file will be named using the run prefix
defined in the inputs file and the corresponding timestep.

Next, run ``burn_cell.py`` using python 3.x, giving the defined run prefix as an argument.
For example::

    python3 burn_cell.py react_aprox13

The ``burn_cell.py`` code will gather information from all of the
output files and compile them into three graphs explained below.

Graphs Output by ``burn_cell.py``
---------------------------------

The file ``run-prefix_logX.png`` and ``run-prefix_logX.eps`` will display a
graph of the chemical abundances as a function of the time, both on
logarithmic scales, for all species involved in the simulation.  An
example of this graph is shown below.

.. figure:: react_aprox13_logX.png
   :alt: An example of a plot output by the burn_cell unit test. This is the logX output cooresponding to the network aprox13.
   :width: 4.5in

   An example of a plot output by the burn_cell unit test. This is the
   logX output cooresponding to the network aprox13.



The file ``run-prefix_ydot.png`` and ``run-prefix_ydot.eps`` will display the
molar fraction (mass fraction / atomic weight) as a function of time,
both on logarithmic scales, for all species involved in the code.

The file ``run-prefix_T-edot.png`` and ``run-prefix_T-edot.eps`` will display
the temperature and the energy generation rate as a function of time.
