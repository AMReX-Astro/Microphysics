
#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <burner.H>
#include <fstream>
#include <iostream>
#include <nse_solver.H>

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

    // Let state.xn equal to nse_state.xn to make sure its in nse state
    for (int n = 0; n < NumSpec; ++n){
      state.xn[n] = NSE_STATE.xn[n];
    }


    // get eos
    //eos(eos_input_rt, state);
    
    if (in_nse(state)){
      std::cout << "We're in NSE. " << std::endl;
    }
    else{
      std::cout << "We're not in NSE. " << std::endl;
    }
    
}
