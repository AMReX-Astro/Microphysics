#ifndef BURN_CELL_H
#define BURN_CELL_H

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <burner.H>
#include <fstream>
#include <iostream>
#include <nse_solver.H>

using namespace unit_test_rp;

AMREX_INLINE
void burn_cell_c()
{

    burn_t state;


    // // Echo initial conditions at burn and fill burn state input

    state.T = temperature;
    state.rho = density;
    state.y_e = ye;

    // set initial chemical potential of proton and neutron

    state.mu_p = mu_p;
    state.mu_n = mu_n;

    // set a reference cell size.

    state.dx = 1.0e6_rt;

    std::cout << "chemical potential of proton is " << mu_p << std::endl;
    std::cout << "chemical potential of neutron is " << mu_n << std::endl;

    // find the  nse state
    bool ye_is_valid = true;
    amrex::Real eps = 1.0e-10_rt;
    auto NSE_STATE = get_actual_nse_state(state, eps, ye_is_valid);

    std::cout << "State Density (g/cm^3): " << state.rho << std::endl;
    std::cout << "State Temperature (K): " << state.T << std::endl;
    std::cout << "electron fraction is " << state.y_e << std::endl;

    std::cout << "NSE state: " << std::endl;
    for (int n = 0; n < NumSpec; ++n) {
      std::cout << short_spec_names_cxx[n] << " : " << NSE_STATE.xn[n] << std::endl;
    }

    // Let state.y equal to nse_state.y to make sure its in nse state

    for (int n = 0; n < NumSpec; ++n){
      state.y[SFS+n] = NSE_STATE.y[SFS+n];
    }

    // get eos
    //eos(eos_input_rt, state);

    if (in_nse(state)){
      std::cout << "We're in NSE. " << std::endl;
    }
    else{
      std::cout << "We're not in NSE. " << std::endl;
    }

    // now test the EOS+NSE function.  The basic idea is the following:
    //
    // 1. find the e corresponding to the current NSE state (eos_input_rt)
    // 2. perturb e and then call the NSE EOS to get the updated T and abar
    // 3. call the NSE function explicitly with the new T to get the X's
    //    (since the NSE EOS doesn't return these).
    // 4. call the EOS with this new T (eos_input_rt) and updated NSE
    //    composition and see if we recover the same perturbed e

    eos_t eos_state;
    eos_state.T = state.T;
    eos_state.rho = state.rho;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = state.y[SFS+n] / state.rho;
    }

    // get the initial e and abar corresponding to this T

    eos(eos_input_rt, eos_state);

    amrex::Real abar_start = eos_state.abar;

    // now perturb e and find the T/abar that are consistent with it
    // and NSE

    amrex::Real e_new = eos_state.e * 1.05;

    eos_state.e = e_new;

    amrex::Real T_new{eos_state.T};
    amrex::Real abar_new{};

    nse_T_abar_from_e(eos_state.rho, eos_state.e, ye,
                      T_new, abar_new,
                      state.mu_p, state.mu_n);

    std::cout << "change in T: " << eos_state.T << " " << T_new <<  std::endl;
    std::cout << "change in abar: " << abar_start << " " << abar_new <<  std::endl;

    // now try calling the EOS with this T and composition and see if
    // we get back the same energy.  we need to recompute the NSE
    // state here since we didn't get the mass fractions from the EOS call

    state.T = T_new;
    auto nse_state = get_actual_nse_state(state, 1.e-10, true);

    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = nse_state.xn[n];
    }

    eos_state.T = T_new;

    eos(eos_input_rt, eos_state);

    std::cout << "recovered energy: " << eos_state.e << " " << e_new << std::endl;

}
#endif
