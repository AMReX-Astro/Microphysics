
#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <burn_type.H>
#include <fstream>
#include <iostream>
#include <nse_solver.H>
#include <react_util.H>

void nse_example_c()
{

    burn_t state;

    Real massfractions[NumSpec] = {-1.0};
    std::cout << massfractions << std::endl;
    // Make sure user set all the mass fractions to values in the interval [0, 1]
    for (int n = 1; n <= NumSpec; ++n) {

        massfractions[n-1] = get_xn(n);

        if (massfractions[n-1] < 0 || massfractions[n-1] > 1) {
            amrex::Error("mass fraction for " + short_spec_names_cxx[n-1] + " not initialized in the interval [0,1]!");
        }

    }

    // Echo initial conditions at burn and fill burn state input

    std::cout << "State Density (g/cm^3): " << density << std::endl;
    std::cout << "State Temperature (K): " << temperature << std::endl;
    for (int n = 0; n < NumSpec; ++n) {
        std::cout << "Mass Fraction (" << short_spec_names_cxx[n] << "): " << massfractions[n] << std::endl;
    }

    state.T = temperature;
    state.rho = density;
    for (int n = 0; n < NumSpec; ++n) {
        state.xn[n] = massfractions[n];
    }

    
    // normalize -- just in case
    
    normalize_abundances_burn(state);

    // compute the initial Ye
    
    state.y_e = y_e;
    
    // composition(state);
    
    std::cout << "electron fraction is " << state.y_e << std::endl;

    // set initial chemical potential of proton and neutron

    state.mu_p = mu_p;
    state.mu_n = mu_n;
    
    const bool assume_ye_valid = true;
    amrex::Real eps = 1.0e-10_rt;
    
    // find the  nse state
    use_hybrid_solver = 1;
    
    auto NSE_STATE = get_actual_nse_state(state, eps, assume_ye_valid);

    std::cout << "After solving: " << std::endl;
    std::cout << "chemical potential of proton is " << state.mu_p << std::endl;
    std::cout << "chemical potential of neutron is " << state.mu_n << std::endl;
    
    std::cout << "NSE state: " << std::endl;
    for (int n = 0; n < NumSpec; ++n) {
      std::cout << short_spec_names_cxx[n] << " : " << NSE_STATE.xn[n] << std::endl;
    }
}
