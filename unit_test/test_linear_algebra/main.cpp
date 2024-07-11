#include <iostream>

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <unit_test.H>
#include <test_linear_algebra.H>

using namespace unit_test_rp;

int main(int argc, char *argv[]) {

  amrex::Initialize(argc, argv);

  init_unit_test();

  // C++ EOS initialization (must be done after init_extern_parameters)
  eos_init(small_temp, small_dens);

  // C++ Network, RHS, screening, rates initialization
  network_init();

  linear_algebra();

  amrex::Finalize();
}
