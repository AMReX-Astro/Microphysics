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

int main(int argc, char *argv[]) {

  amrex::Initialize(argc, argv);

  std::cout << "starting the single zone burn..." << std::endl;

  ParmParse ppa("amr");

  init_unit_test();

  // C++ EOS initialization (must be done after Fortran eos_init and init_extern_parameters)
  eos_init(small_temp, small_dens);

  // C++ Network, RHS, screening, rates initialization
  network_init();

  // amrex::Real mu_p;
  // amrex::Real mu_n;
  
  // std::cout << "enter initial guess for mu_p" << std::endl;
  // std::cin >> mu_p;
  // std::cout << "enter initial guess for mu_n" << std::endl;
  // std::cin >> mu_n;
  
  burn_cell_c();

  amrex::Finalize();
}
