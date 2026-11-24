#include <cstring>
#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
using namespace amrex;

#include <burn_cell.H>
#include <eos.H>
#include <extern_parameters.H>
#include <network.H>
#include <unit_test.H>

int main(int argc, char *argv[]) {

  amrex::Initialize(argc, argv);
  int success = 0;

  {
    // check that correct input file is provided
    ParmParse const pp("unit_test");
    std::string const run_prefix = "burn_cell_metal_chem_";
    std::string input_run_prefix;
    pp.query("run_prefix", input_run_prefix);
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(run_prefix == input_run_prefix,
                                     "input file is missing or incorrect!");

    std::cout << "starting the single zone burn..." << std::endl;

    ParmParse ppa("amr");

    init_unit_test();

    // C++ EOS initialization (must be done after
    // init_extern_parameters)
    eos_init(unit_test_rp::small_temp, unit_test_rp::small_dens);

    // C++ Network, RHS, screening, rates initialization
    network_init();
    actual_network_init();

    success = burn_cell_c();
  }

  amrex::Finalize();

  return (!success);
}
