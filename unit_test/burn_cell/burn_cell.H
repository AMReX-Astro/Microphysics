#ifndef BURN_CELL_H
#define BURN_CELL_H

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <burner.H>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <react_util.H>

void burn_cell_c()
{

    Real massfractions[NumSpec] = {-1.0};

    // Make sure user set all the mass fractions to values in the interval [0, 1]
    for (int n = 1; n <= NumSpec; ++n) {

        massfractions[n-1] = get_xn(n);

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

    burn_t burn_state;

    eos_t eos_state;
    eos_state.rho = density;
    eos_state.T = temperature;
    for (int n = 0; n < NumSpec; n++) {
        eos_state.xn[n] = massfractions[n];
    }
#ifdef AUX_THERMO
    set_aux_comp_from_X(eos_state);
#endif

    eos(eos_input_rt, eos_state);

    burn_state.rho = eos_state.rho;
    burn_state.T = eos_state.T;
    for (int n = 0; n < NumSpec; ++n) {
        burn_state.xn[n] = massfractions[n];
    }
#if NAUX_NET > 0
    for (int n = 0; n < NumAux; ++n) {
        burn_state.aux[n] = eos_state.aux[n];
    }
#endif

    burn_state.i = 0;
    burn_state.j = 0;
    burn_state.k = 0;
    burn_state.T_fixed = -1.0_rt;

    // normalize -- just in case

    if (! skip_initial_normalization) {
        normalize_abundances_burn(burn_state);
    }

    // call the EOS to set initial e -- it actually doesn't matter to
    // the burn but we need to keep track of e to get a valid
    // temperature for the burn if we substep

    eos(eos_input_rt, burn_state);

    // output just the instantaneous RHS

    Array1D<Real, 1, neqs> ydot;
    actual_rhs(burn_state, ydot);

    std::cout << "RHS at t = 0" << std::endl;
    for(int n = 0; n < NumSpec; ++n){
        std::string element = short_spec_names_cxx[n];
        std::cout << std::setw(6) << element << " " << ydot(n+1) << std::endl;
    }

    // output initial burn type data
    Real time = 0.0_rt;

    std::ofstream state_over_time("state_over_time.txt");

    // we will divide the total integration time into nsteps that are
    // logarithmically spaced

    if (tfirst == 0.0_rt) {
        if (nsteps == 1) {
            tfirst = tmax;
        } else {
            tfirst = tmax / nsteps;
        }
    }
    Real dlogt = 0.0_rt;
    if (nsteps == 1) {
        dlogt = (std::log10(tmax) - std::log10(tfirst));
    } else {
        dlogt = (std::log10(tmax) - std::log10(tfirst)) / (nsteps - 1);
    }

    // save the initial state -- we'll use this to determine
    // how much things changed over the entire burn

    burn_t burn_state_in = burn_state;

    // output the data in columns, one line per timestep


    state_over_time << std::setw(25) << "# Time";
    state_over_time << std::setw(25) << "Temperature";
    for(int x = 0; x < NumSpec; ++x){
        std::string element = short_spec_names_cxx[x];
	state_over_time << std::setw(25) << element;
    }
    state_over_time << std::endl;
    state_over_time << std::setprecision(15);

    Real t = 0.0;

    state_over_time << std::setw(25) << t;
    state_over_time << std::setw(25) << burn_state.T;
    for (int x = 0; x < NumSpec; ++x){
        state_over_time << std::setw(25) << burn_state.xn[x];
    }
    state_over_time << std::endl;


    // store the initial internal energy -- we'll update this after
    // each substep

    Real energy_initial = burn_state.e;

    // loop over steps, burn, and output the current state

    int nstep_int = 0;

    for (int n = 0; n < nsteps; n++){

        // compute the time we wish to integrate to

        Real tend = std::pow(10.0_rt, std::log10(tfirst) + dlogt * n);
        Real dt = tend - t;

    	burner(burn_state, dt);

        nstep_int += burn_state.n_step;

        // state.e represents the change in energy over the burn (for
        // just this sybcycle), so turn it back into a physical energy

        burn_state.e += energy_initial;

        // reset the initial energy for the next subcycle

        energy_initial = burn_state.e;

        // get the updated T

        if (call_eos_in_rhs) {
            eos(eos_input_re, burn_state);
        }

        t += dt;

	state_over_time << std::setw(25) << t;
	state_over_time << std::setw(25) << burn_state.T;
	for (int x = 0; x < NumSpec; ++x){
	     state_over_time << std::setw(25) << burn_state.xn[x];
	}
	state_over_time << std::endl;
    }
    state_over_time.close();

    // output diagnostics to the terminal

    std::cout << "------------------------------------" << std::endl;
    std::cout << "successful? " << burn_state.success << std::endl;
    std::cout << " - Hnuc = " << (burn_state.e - burn_state_in.e) / tmax << std::endl;
    std::cout << " - added e = " << burn_state.e - burn_state_in.e << std::endl;
    std::cout << " - final T = " << burn_state.T << std::endl;


    std::cout << "------------------------------------" << std::endl;
    std::cout << "e initial = " << burn_state_in.e << std::endl;
    std::cout << "e final =   " << burn_state.e << std::endl;


    std::cout << "------------------------------------" << std::endl;
    std::cout << "new mass fractions: " << std::endl;
    for (int n = 0; n < NumSpec; ++n) {
        std::string element = short_spec_names_cxx[n];
        std::cout << element << " " << burn_state.xn[n] << std::endl;
    }

    std::cout << "------------------------------------" << std::endl;
    std::cout << "species creation rates: " << std::endl;
    for (int n = 0; n < NumSpec; ++n) {
        std::cout << "omegadot(" << short_spec_names_cxx[n] << "): "
                  << (burn_state.xn[n] - burn_state_in.xn[n]) / tmax << std::endl;
    }

    std::cout << "number of steps taken: " << nstep_int << std::endl;

}
#endif
