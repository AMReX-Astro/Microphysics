#include <iostream>
#include <cstring>
#include <vector>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
using namespace amrex;

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <eos_cell.H>
#include <eos_cell_F.H>

int main(int argc, char *argv[]) {

  amrex::Initialize(argc, argv);

  std::cout << "calling the EOS on a single zone state..." << std::endl;

  ParmParse ppa("amr");

  std::string probin_file = "probin";

  ppa.query("probin_file", probin_file);

  std::cout << "probin = " << probin_file << std::endl;

  const int probin_file_length = probin_file.length();
  Vector<int> probin_file_name(probin_file_length);

  for (int i = 0; i < probin_file_length; i++)
    probin_file_name[i] = probin_file[i];

  init_unit_test(probin_file_name.dataPtr(), &probin_file_length);

  // Copy extern parameters from Fortran to C++
  init_extern_parameters();

  // C++ EOS initialization (must be done after Fortran eos_init and init_extern_parameters)
  eos_init();

  // C++ Network, RHS, screening, rates initialization
  network_init();

  eos_cell_c();

  amrex::Finalize();
}
