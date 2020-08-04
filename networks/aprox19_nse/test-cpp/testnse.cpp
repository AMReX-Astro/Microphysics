#include <iostream>
#include <actual_network.H>

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Array.H>
#include <AMReX_REAL.H>

#include "extern_parameters.H"
#include "nse.H"

using namespace amrex;


void do_nse_cxx() {

  // this routine simply calls the table interpolation with a fixed
  // input, so we can compare the results with those we get from the
  // entire Microphysics machinery.


  constexpr Real temp = 1.e9_rt;
  constexpr Real rho = 1.e9_rt;
  constexpr Real ye = 0.46_rt;

  Real abar, dq, dyedt, X[NumSpec];

  nse_interp(temp, rho, ye, abar, dq, dyedt, X);

  std::cout << "temp:    " << temp << std::endl;
  std::cout << "rho:     " << rho << std::endl;
  std::cout << "ye:      " << ye << std::endl;
  std::cout << "abar:    " << abar << std::endl;
  std::cout << "be/a:    " << dq << std::endl;
  std::cout << "dyedt:   " << dyedt << std::endl;
  std::cout << "X:       ";
  for (int n = 0; n < NumSpec; ++n) {
    std::cout << X[n] << " ";
  }
  std::cout << std::endl;
}

