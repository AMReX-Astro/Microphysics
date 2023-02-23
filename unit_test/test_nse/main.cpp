#include <iostream>
#include <cstring>
#include <vector>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
using namespace amrex;

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <nse_example.H>
#include <unit_test.H>

int main(int argc, char *argv[]) {

  amrex::Initialize(argc, argv);

  std::cout << "starting the single zone burn..." << std::endl;

  ParmParse ppa("amr");

  init_unit_test();

  // C++ EOS initialization (must be done after Fortran eos_init and init_extern_parameters)
  eos_init(small_temp, small_dens);

  // C++ Network, RHS, screening, rates initialization
  network_init();

  nse_example_c();

  amrex::Finalize();
}
