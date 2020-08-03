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


  constexpr Real t9 = 1.0_rt;
  constexpr Real rho = 1.e9_rt;
  constexpr Real ye = 0.46;

  Real abar, dq, dyedt, X[NumSpec];

  nse_interp(t9, rho, ye, abar, dq, dyedt, X);

  std::cout << "      t9    " << "        rho    " << "       ye     " <<
    "     abar     " << "      be/a    " << "      dyedt  " << std::endl;

  std::cout << t9 << " " << rho << " " << ye << " " <<  abar << " " << dq << " " << dyedt;
  for (int n = 0; n < NumSpec; ++n) {
    std::cout << " " << X[n];
  }
  std::cout << std::endl;
}

