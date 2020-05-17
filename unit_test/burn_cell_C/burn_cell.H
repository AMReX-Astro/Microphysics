
#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <burner.H>

void burn_cell_c()
{

    burn_t state;

    // Set mass fractions to sanitize inputs for them
    Real massfractions[NumSpec];
    for (int n = 0; n < NumSpec; ++n) {
        massfractions[n] = -1.0e0_rt;
    }

    // Make sure user set all the mass fractions to values in the interval [0, 1]
    for (int n = 1; n <= NumSpec; ++n) {
        switch (n) {

        case 1:
            massfractions[n-1] = X1;
            break;
        case 2:
            massfractions[n-1] = X2;
            break;
        case 3:
            massfractions[n-1] = X3;
            break;
        case 4:
            massfractions[n-1] = X4;
            break;
        case 5:
            massfractions[n-1] = X5;
            break;
        case 6:
            massfractions[n-1] = X6;
            break;
        case 7:
            massfractions[n-1] = X7;
            break;
        case 8:
            massfractions[n-1] = X8;
            break;
        case 9:
            massfractions[n-1] = X9;
            break;
        case 10:
            massfractions[n-1] = X10;
            break;
        case 11:
            massfractions[n-1] = X11;
            break;
        case 12:
            massfractions[n-1] = X12;
            break;
        case 13:
            massfractions[n-1] = X13;
            break;
        case 14:
            massfractions[n-1] = X14;
            break;
        case 15:
            massfractions[n-1] = X15;
            break;
        case 16:
            massfractions[n-1] = X16;
            break;
        case 17:
            massfractions[n-1] = X17;
            break;
        case 18:
            massfractions[n-1] = X18;
            break;
        case 19:
            massfractions[n-1] = X19;
            break;
        case 20:
            massfractions[n-1] = X20;
            break;
        case 21:
            massfractions[n-1] = X21;
            break;

        }

        if (massfractions[n-1] < 0 || massfractions[n-1] > 1) {
            amrex::Error("mass fraction for " + short_spec_names_cxx[n-1] + " not initialized in the interval [0,1]!");
        }

    }

    // Echo initial conditions at burn and fill burn state input
    std::cout << "Maximum Time (s): " << tmax << std::endl;
    std::cout << "State Density (g/cm^3): " << density << std::endl;
    std::cout << "State Temperature (K): " << temperature << std::endl;
    for (int n = 0; n < NumSpec; ++n) {
        std::cout << "Mass Fraction (" << short_spec_names_cxx[n] << "): " << massfractions[n] << std::endl;
    }

    state.T   = temperature;
    state.rho = density;
    for (int n = 0; n < NumSpec; ++n) {
        state.xn[n] = massfractions[n];
    }

    // normalize -- just in case
    //normalize_abundances_burn(state);

    // Initialize initial energy to zero
    state.e = 0.0_rt;
    Real energy = 0.0_rt;

    // output initial burn type data
    Real time = 0.0_rt;

    // call the EOS to set initial e
    eos_t eos_state_in;
    burn_to_eos(state, eos_state_in);
    eos(eos_input_rt, eos_state_in);
    eos_to_burn(eos_state_in, state);

    Real dt = tmax;
    burn_t state_in = state;
    burner(state, dt);
    energy += state.e;

    // call the EOS to check consistency of integrated e
    eos_t eos_state_out;
    burn_to_eos(state, eos_state_out);
    eos(eos_input_rt, eos_state_out);

    std::cout << "------------------------------------" << std::endl;
    std::cout << "successful? " << state.success << std::endl;
    std::cout << " - Hnuc = " << state.e / dt << std::endl;
    std::cout << " - integrated e = " << eos_state_in.e + energy << std::endl;
    std::cout << " - EOS e(rho, T) = " << eos_state_out.e << std::endl;
    std::cout << " - integrated/EOS percent diff. = " << 100.0e0_rt * (eos_state_in.e + energy - eos_state_out.e) / eos_state_out.e << std::endl;

    // output burn type data
    //write_burn_t(state)

    std::cout << "------------------------------------" << std::endl;
    std::cout << "EOS e(rho, T) initial = " << eos_state_in.e << std::endl;
    std::cout << "EOS e(rho, T) final =   " << eos_state_out.e << std::endl;
    Real eos_energy_generated = eos_state_out.e - eos_state_in.e;
    std::cout << "EOS e(rho, T) generated = " << eos_energy_generated << std::endl;
    Real eos_energy_rate = (eos_state_out.e - eos_state_in.e) / tmax;
    std::cout << "EOS e(rho, T) generation rate = " << eos_energy_rate << std::endl;
    std::cout << "Integrator total generated energy: " << energy << std::endl;
    std::cout << "Integrator average energy generation rate: " << energy/tmax << std::endl;
    std::cout << "(integrator - EOS)/EOS percent diff for generated energy: " << 100.0e0_rt * (energy - eos_energy_generated) / eos_energy_generated << std::endl;
    std::cout << "(integrator - EOS)/EOS percent diff for energy gen. rate: " << 100.0e0_rt * (energy / tmax - eos_energy_rate) / eos_energy_rate << std::endl;

    for (int n = 0; n < NumSpec; ++n) {
        std::cout << "omegadot(" << short_spec_names_cxx[n] << "): "
                  << (state.xn[n] - state_in.xn[n]) / dt << std::endl;
    }

    for (int n = 0; n < NumSpec; ++n) {
        std::cout << "delta(" << short_spec_names_cxx[n] << "): "
                  << (state.xn[n] - state_in.xn[n]) << std::endl;
    }

    for (int n = 0; n < NumSpec; ++n) {
        std::cout << "percent change(" << short_spec_names_cxx[n] << "): "
                  << 100.e0_rt * (state.xn[n] - state_in.xn[n]) / state_in.xn[n] << std::endl;
    }

}
