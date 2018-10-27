*****************
Burn Cell Utility
*****************

burn_cell
=========

Getting Started
---------------

The burn_cell codes are located in Microphysics/unit_test/burn_cell. To run a simulation, ensure that both an input file and an initial conditions file have been created and are in the same directory as the executable.

Input File
----------

These files are typically named as inputs_burn_network where network is the network you wish to use for your testing.

The structure of this file is is fairly self-explanitory.
The run prefix defined should be unique to the tests that will be run as they will be used to identify all of the output files. Typically, the run prefix involves the name of the network being tested.
The atol variables define absolute tolerances of the ordinary differential equations and the rtol variables define the relative tolerances.
The second section of the input file collects the inputs that main.f90 asks for so that the user does not have to input all 5\ :math:`+` parameters that are required everytime the test is run.
Each input required is defined and initialized on the lines following &cellparams.
The use of the parameters is defined in Table \ `[tab:init-structure] <#tab:init-structure>`__.

.. raw:: latex

   \centering

.. table:: The definition of parameters used in the burn_cell unit tests and specified in the second half of each inputs file.

   +-------------+----------------------------------------+
   |             | Maximum Time (s)                       |
   +-------------+----------------------------------------+
   | numsteps    | Number of time subdivisions            |
   +-------------+----------------------------------------+
   |             | State Density (:math:`\frac{g}{cm^3}`) |
   +-------------+----------------------------------------+
   | temperature | State Temperature (K)                  |
   +-------------+----------------------------------------+
   |             | Mass Fraction for element i            |
   +-------------+----------------------------------------+

Running the Code
----------------

To run the code, enter the burn_cell directory and run ./main.Linux.gfortran.exe with the inputs file as an argument.
For example: ./main.Linux.gfortran.exe inputs_burn_aprox13

For however many events are run, defined as numsteps in the inputs file, the code will output that many files into a new directory titled run_prefix_output where run_prefix is the run prefix defined in the inputs file.
Each output file will be named using the run prefix defined in the inputs file and the corresponding timestep.

Next, run burn_cell.py using Python3 and listing the defined run prefix as an argument.
For example: python3 burn_cell.py react_aprox13.
The burn_cell code will gather information from all of the output files and compile them into three graphs explained below.

Graphs Output by burn_cell.py
-----------------------------

The file run-prefix_logX.png and run-prefix_logX.eps will display a graph of the chemical abundances as a function of the time, both on logarithmic scales, for all species involved in the simulation.
An example of this graph is shown in Figure \ `[fig:aprox13_logX] <#fig:aprox13_logX>`__.

.. raw:: latex

   \centering

.. figure:: react_aprox13_logX.png
   :alt: An example of a plot output by the burn_cell unit test. This is the logX output cooresponding to the network aprox13.
   :width: 4.5in

   An example of a plot output by the burn_cell unit test. This is the logX output cooresponding to the network aprox13.

[fig:aprox13_logX]

The file run-prefix_ydot.png and run-prefix_ydot.eps will display the Moller fraction (mass fraction / atomic weight) as a function of time, both on logarithmic scales, for all species involved in the code.

The file run-prefix_T-edot.png and run-prefix_T-edot.eps will display the Temperature and the Energy Generation Rate as a function of time.
