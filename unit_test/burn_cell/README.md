# `burn_cell`

This is a single zone test burn that will output a time-history of the
nucleosynthesis.

# Using

  * make sure you have AMReX and that `AMREX_HOME` is set to point to
    the `amrex/` directory.

  * Compile the unit test:

    `make -j NETWORK_DIR=ignition_reaclib/URCA-simple INTEGRATOR_DIR=VODE90`

    This builds an executable that uses the `VODE90` ODE integrator to
    solve the differential equations corresponding to the reaction
    network `ignition_reaclib/URCA-simple`.  You can experiment with
    other networks and integrators too.

    The executable will be named something like `main.Linux.gfortran.exe`

  * Now we set the inputs.  The executable takes an input file and
    also prompts for the end time for the integration and the initial
    conditions of density, temperature, and isotope abundances.  We
    can redirect STDIN from a file, as shown below:

    `./main.Linux.gfortran.exe inputs_burn_urca.VODE < init_urca`

    That init and inputs file are specific to the network compiled
    above, so if you want to use `aprox13` as your network, you'll have
    to make inputs and init files for it. 

  * A lot of output will be produced, of the form `react_urca_ABCDEF`
    that correspond to integration time steps.

    The first part of that filename, `react_urca` is the "run prefix"
    that the python plotting script needs. The plotting script doesn't
    need the init file.

  * To make plots of the output, we can use the python scripts like

    `python burn_cell.py react_urca`


