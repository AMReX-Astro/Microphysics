#include <iostream>

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <jac_cell.H>
#include <unit_test.H>

using namespace unit_test_rp;

int main(int argc, char *argv[]) {

  amrex::Initialize(argc, argv);

  std::cout << "starting the single zone burn..." << std::endl;

  init_unit_test();

  // C++ EOS initialization (must be done after Fortran eos_init and init_extern_parameters)
  eos_init(small_temp, small_dens);

  // C++ Network, RHS, screening, rates initialization
  network_init();

  jac_cell_c();

  amrex::Finalize();
}
