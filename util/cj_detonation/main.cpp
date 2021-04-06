#include <iostream>
#include <cstring>
#include <vector>

#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
using namespace amrex;

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <cj_det.H>
#include <unit_test.H>
#include <actual_rhs.H>

int main(int argc, char *argv[]) {

  amrex::Initialize(argc, argv);

  std::cout << "starting the CJ Det solve..." << std::endl;

  ParmParse ppa("amr");

  std::string probin_file = "probin";

  ppa.query("probin_file", probin_file);

  std::cout << "probin = " << probin_file << std::endl;

  const int probin_file_length = probin_file.length();
  Vector<int> probin_file_name(probin_file_length);

  for (int i = 0; i < probin_file_length; i++)
    probin_file_name[i] = probin_file[i];

  init_unit_test(probin_file_name.dataPtr(), &probin_file_length);

  // C++ EOS initialization (must be done after Fortran eos_init and
  // init_extern_parameters)
  eos_init(small_temp, small_dens);

  // C++ Network, RHS, screening, rates initialization

  network_init();

  const Real rho_min_fac = 0.9_rt;
  const Real rho_max_fac = 10.0_rt;
  const int npts_ad = 150;
  const int npts = 100;


  // set the unburned (fuel) state

  eos_t eos_state_fuel;

  eos_state_fuel.rho = 1.e7_rt;
  eos_state_fuel.T = 1.e8_rt;
  for (int n = 0; n < NumSpec; n++) {
      eos_state_fuel.xn[n] = smallx;
  }
  eos_state_fuel.xn[0] = 1.0_rt - (NumSpec - 1) * smallx;

  eos(eos_input_rt, eos_state_fuel);


  // set the ash composition

  eos_t eos_state_ash = eos_state_fuel;
  for (int n = 0; n < NumSpec; n++) {
      eos_state_ash.xn[n] = smallx;
  }
  eos_state_ash.xn[NumSpec-1] = 1.0_rt - (NumSpec - 1) * smallx;

  // get the q value -- we need the change in molar fractions
  Array1D<Real, 1, NumSpec> dymol;

  for (int n = 1; n <= NumSpec; n++) {
      dymol(n) = eos_state_ash.xn[n-1] * aion_inv[n-1] -
                 eos_state_fuel.xn[n-1] * aion_inv[n-1];
  }

  Real q_burn;
  ener_gener_rate(dymol, q_burn);

  // store the shock adiabat and the detonation adiabat

  Real rho_min = rho_min_fac * eos_state_fuel.rho;
  Real rho_max = rho_max_fac * eos_state_fuel.rho;
  Real dlogrho = (std::log10(rho_max) - std::log10(rho_min)) / static_cast<Real>(npts-1);

  // initial guess

  eos_state_ash.T = eos_state_fuel.T;


  // now let's get the CJ velocity

  cj_cond(eos_state_fuel, eos_state_ash, q_burn);

  // we get this from the mass flux: rho_1 v_1 = j

  Real D_cj = (1.0_rt / eos_state_fuel.rho) *
      std::sqrt((eos_state_ash.p - eos_state_fuel.p) /
                (1.0_rt / eos_state_fuel.rho - 1.0_rt / eos_state_ash.rho));

  Real rho_cj = eos_state_ash.rho;
  Real p_cj = eos_state_ash.p;

  // output info along with points on the Hugoniot

  std::ofstream hout("hugoniot.txt");

  hout << "# initial rho = " << eos_state_fuel.rho << std::endl;
  hout << "# initial p = " << eos_state_fuel.p << std::endl;
  hout << "# ash rho = " << eos_state_ash.rho << std::endl;
  hout << "# ash p = " << eos_state_ash.p << std::endl;
  hout << "# CJ speed = " << D_cj << std::endl;

  // test

  Real cs_det = std::sqrt(eos_state_ash.gam1 * eos_state_ash.p / eos_state_ash.rho);

  // this checks that the slope of the Rayleigh line is rho_2 * cs_2
  // std::cout << eos_state_ash.rho * cs_det << " " << eos_state_fuel.rho * D_cj << std::endl;

  for (int n = 0; n < npts_ad; n++) {

      eos_state_ash.rho = std::pow(10.0_rt, std::log10(rho_min) + static_cast<Real>(n) * dlogrho);

      int istatus;
      adiabat(eos_state_fuel, eos_state_ash, 0.0_rt, istatus);
      Real p2_shock = eos_state_ash.p;

     if (istatus == -1) {
        amrex::Error("error in adiabat solve");
     }

     adiabat(eos_state_fuel, eos_state_ash, q_burn, istatus);
     Real p2_det = eos_state_ash.p;

     if (istatus == -1) {
        amrex::Error("error in adiabat solve");
     }

     hout << eos_state_ash.rho << " " << p2_shock << " " << p2_det << std::endl;
  }


  hout.close();

  return 0;

}
