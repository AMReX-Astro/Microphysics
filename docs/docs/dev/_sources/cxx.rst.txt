**************
C++ interfaces
**************

There are C++ interfaces for the equation of state.  The following conventions are used:

  * C++ support requires AMReX

  * We assume that Fortran initialization is done first (runtime
    paramters, EOS, ...) and we get what information we can on
    parameters from Fortran.

  * The EOS logic is completely duplicated in C++ so we don't need to
    call any Fortran.

  * The nuclei information is provided in a .net file and used to
    generate both the Fortran ``actual_network.F90`` module and the C++
    ``actual_network.H`` header file.

  * We don't attempt to pass the C++ struct directly into Fortran.
