#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include <AMReX_Config.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
using namespace amrex;

#include <burn_cell.H>
#include <eos.H>
#include <extern_parameters.H>
#include <network.H>
#include <unit_test.H>

int main(int argc, char *argv[]) {

  std::string compare_final_state_file;
  std::vector<char *> amrex_argv;
  amrex_argv.reserve(static_cast<std::size_t>(argc));
  amrex_argv.push_back(argv[0]);

  for (int i = 1; i < argc; ++i) {
    std::string arg(argv[i]);
    if (arg == "--compare-final-state") {
      if (i + 1 >= argc) {
        std::cerr << "ERROR: --compare-final-state requires a file path"
                  << std::endl;
        return 1;
      }
      compare_final_state_file = argv[++i];
      continue;
    }

    constexpr const char *prefix = "--compare-final-state=";
    if (arg.rfind(prefix, 0) == 0) {
      compare_final_state_file = arg.substr(std::strlen(prefix));
      if (compare_final_state_file.empty()) {
        std::cerr << "ERROR: --compare-final-state requires a non-empty file path"
                  << std::endl;
        return 1;
      }
      continue;
    }

    amrex_argv.push_back(argv[i]);
  }

  int amrex_argc = static_cast<int>(amrex_argv.size());
  char **amrex_argv_ptr = amrex_argv.data();
  amrex::Initialize(amrex_argc, amrex_argv_ptr);
  int success = 0;

  {
    // check that correct input file is provided
    ParmParse const pp("unit_test");
    std::string const run_prefix = "burn_cell_primordial_chem_";
    std::string input_run_prefix;
    pp.query("run_prefix", input_run_prefix);
    int n_zones = 128;
    pp.query("n_zones", n_zones);
    amrex::Print() << "n_zones = " << n_zones << "\n";
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(run_prefix == input_run_prefix,
                                     "input file is missing or incorrect!");

    ParmParse ppa("amr");

    init_unit_test();

    // C++ EOS initialization (must be done after
    // init_extern_parameters)
    eos_init(unit_test_rp::small_temp, unit_test_rp::small_dens);

    // C++ Network, RHS, screening, rates initialization
    network_init();

    bool enable_gpu = false;
#ifdef AMREX_USE_GPU
    enable_gpu = true;
#endif

    std::cout << "\nstarting the multi-zone burn..." << std::endl;
    if (!compare_final_state_file.empty()) {
      std::cout << "saved-state comparison file: " << compare_final_state_file
                << std::endl;
    }
    int multi_success =
        burn_cell_multi_c(n_zones, enable_gpu, compare_final_state_file);

    success = multi_success;
  }

  amrex::Finalize();

  return (!success);
}
