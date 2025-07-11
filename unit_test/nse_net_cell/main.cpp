#include <iostream>
#include <cstring>
#include <vector>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <burn_cell.H>
#include <unit_test.H>

int main(int argc, char *argv[]) {

  amrex::Initialize(argc, argv);

  std::cout << "starting the single zone burn..." << std::endl;

  amrex::ParmParse ppa("amr");

  init_unit_test();

  // C++ EOS initialization (must be done after init_extern_parameters)
  eos_init(unit_test_rp::small_temp, unit_test_rp::small_dens);

  // C++ Network, RHS, screening, rates initialization
  network_init();

  burn_cell_c();

  amrex::Finalize();
}
