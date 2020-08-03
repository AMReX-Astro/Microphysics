#include <new>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <iomanip>

#ifndef WIN32
#include <unistd.h>
#endif

#include <AMReX_CArena.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>

#include <time.h>

#include "extern_parameters.H"
#include "eos.H"
#include "network.H"

extern "C"
{
  void do_nse_F();
  void do_initialization();
}

void do_nse_cxx();


std::string inputs_name = "";

int
main (int   argc,
      char* argv[])
{

    //
    // Make sure to catch new failures.
    //
    amrex::Initialize(argc,argv);

    // save the inputs file name for later
    if (argc > 1) {
      if (!strchr(argv[1], '=')) {
        inputs_name = argv[1];
      }
    }

    // initialize the runtime parameters in Fortran and the Fortran microphysics
    do_initialization();

    // Copy extern parameters from Fortran to C++
    init_extern_parameters();

    // C++ EOS initialization (must be done after Fortran eos_init and init_extern_parameters)
    eos_init();

#ifdef CXX_REACTIONS
    // C++ Network, RHS, screening, rates initialization
    network_init();
#endif

    std::cout << "Fortran: " << std::endl;
    do_nse_F();

#ifdef CXX_REACTIONS
    std::cout << "C++: " << std::endl;
    do_nse_cxx();
#endif

    amrex::Finalize();

    return 0;
}
