#include <iostream>
#include <cstring>
#include <vector>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
using namespace amrex;

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <burn_cell.H>
#include <unit_test.H>

using namespace unit_test_rp;

int main(int argc, char *argv[]) {

  amrex::Initialize(argc, argv);

  std::cout << "starting the single zone burn..." << std::endl;

  ParmParse ppa("amr");

  init_unit_test();

  // C++ EOS initialization (must be done after Fortran eos_init and init_extern_parameters)
  Real st{small_temp};
  Real sd{small_dens};
  eos_init(st, sd);

  // C++ Network, RHS, screening, rates initialization
  network_init();

  burn_cell_c();

  amrex::Finalize();
}
