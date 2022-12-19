#ifndef BURN_CELL_H
#define BURN_CELL_H

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <burner.H>
#include <fstream>
#include <iostream>
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


#ifdef AUX_THERMO
    Real auxdata[NumAux] = {-1.0};

    // Make sure user set all the mass fractions to values in the interval [0, 1]
    for (int n = 1; n <= NumAux; ++n) {
        switch (n) {

        case 1:
            auxdata[n-1] = Aux1;
            break;
        case 2:
            auxdata[n-1] = Aux2;
            break;
        case 3:
            auxdata[n-1] = Aux3;
            break;

        }

    }
#endif


    Real adv_species[NumSpec] = {0.0};

    // Make sure user set all the mass fractions to values in the interval [0, 1]
    for (int n = 1; n <= NumSpec; ++n) {
        switch (n) {

        case 1:
            adv_species[n-1] = Adv_X1;
            break;
        case 2:
            adv_species[n-1] = Adv_X2;
            break;
        case 3:
            adv_species[n-1] = Adv_X3;
            break;
        case 4:
            adv_species[n-1] = Adv_X4;
            break;
        case 5:
            adv_species[n-1] = Adv_X5;
            break;
        case 6:
            adv_species[n-1] = Adv_X6;
            break;
        case 7:
            adv_species[n-1] = Adv_X7;
            break;
        case 8:
            adv_species[n-1] = Adv_X8;
            break;
        case 9:
            adv_species[n-1] = Adv_X9;
            break;
        case 10:
            adv_species[n-1] = Adv_X10;
            break;
        case 11:
            adv_species[n-1] = Adv_X11;
            break;
        case 12:
            adv_species[n-1] = Adv_X12;
            break;
        case 13:
            adv_species[n-1] = Adv_X13;
            break;
        case 14:
            adv_species[n-1] = Adv_X14;
            break;
        case 15:
            adv_species[n-1] = Adv_X15;
            break;
        case 16:
            adv_species[n-1] = Adv_X16;
            break;
        case 17:
            adv_species[n-1] = Adv_X17;
            break;
        case 18:
            adv_species[n-1] = Adv_X18;
            break;
        case 19:
            adv_species[n-1] = Adv_X19;
            break;
        case 20:
            adv_species[n-1] = Adv_X20;
            break;
        case 21:
            adv_species[n-1] = Adv_X21;
            break;
        case 22:
            adv_species[n-1] = Adv_X22;
            break;
        case 23:
            adv_species[n-1] = Adv_X23;
            break;
        case 24:
            adv_species[n-1] = Adv_X24;
            break;
        case 25:
            adv_species[n-1] = Adv_X25;
            break;
        case 26:
            adv_species[n-1] = Adv_X26;
            break;
        case 27:
            adv_species[n-1] = Adv_X27;
            break;
        case 28:
            adv_species[n-1] = Adv_X28;
            break;
        case 29:
            adv_species[n-1] = Adv_X29;
            break;
        case 30:
            adv_species[n-1] = Adv_X30;
            break;
        case 31:
            adv_species[n-1] = Adv_X31;
            break;
        case 32:
            adv_species[n-1] = Adv_X32;
            break;
        case 33:
            adv_species[n-1] = Adv_X33;
            break;
        case 34:
            adv_species[n-1] = Adv_X34;
            break;
        case 35:
            adv_species[n-1] = Adv_X35;
            break;

        }

    }


#ifdef AUX_THERMO
    Real adv_aux[NumAux] = {0.0};

    // Make sure user set all the mass fractions to values in the interval [0, 1]
    for (int n = 1; n <= NumAux; ++n) {
        switch (n) {

        case 1:
            adv_aux[n-1] = Adv_Aux1;
            break;
        case 2:
            adv_aux[n-1] = Adv_Aux2;
            break;
        case 3:
            adv_aux[n-1] = Adv_Aux3;
            break;

        }

    }
#endif


    // Echo initial conditions at burn and fill burn state input

    std::cout << "Maximum Time (s): " << tmax << std::endl;
    std::cout << "State Density (g/cm^3): " << density << std::endl;
    std::cout << "State Temperature (K): " << temperature << std::endl;
    std::cout << "State rhoe (erg/cm^3): " << rhoe << std::endl;

    for (int n = 0; n < NumSpec; ++n) {
        std::cout << "Mass Fraction (" << short_spec_names_cxx[n] << "): " << massfractions[n] << std::endl;
    }

    // load the state -- we need to create the conserved quantities.
    // Since momentum and total energy don't matter, we'll assume that
    // the velocity is zero and just zero out their advective terms

    burn_t burn_state;

    eos_t eos_state;
    eos_state.rho = density;
    eos_state.T = temperature;

    for (int n = 0; n < NumSpec; n++) {
      eos_state.xn[n] = massfractions[n];
    }
#ifdef AUX_THERMO
    if (recompute_aux) {
        set_aux_comp_from_X(eos_state);
    } else {
        for (int n = 0; n < NumAux; ++n) {
            eos_state.aux[n] = auxdata[n];
        }
    }
#endif

    if (rhoe < 0) {
        eos(eos_input_rt, eos_state);
        rhoe = density * eos_state.e;
    }

    burn_state.rho = eos_state.rho;
    burn_state.T = eos_state.T;

    for (int n = 0; n < NumSpec; n++) {
      burn_state.y[SFS+n] = burn_state.rho * eos_state.xn[n];
    }

#if NAUX_NET > 0
    for (int n = 0; n < NumAux; n++) {
        burn_state.y[SFX+n] = burn_state.rho * eos_state.aux[n];
    }
#endif

#if defined(SDC_EVOLVE_ENERGY)
    burn_state.y[SRHO] = density;

    burn_state.y[SMX] = 0.0;
    burn_state.y[SMY] = 0.0;
    burn_state.y[SMZ] = 0.0;

    burn_state.y[SEINT] = rhoe;
    burn_state.y[SEDEN] = burn_state.y[SEINT];

#elif defined(SDC_EVOLVE_ENTHALPY)
    burn_state.y[SENTH] = burn_state.rho * eos_state.h;

    burn_state.p0 = eos_state.p;
#endif

#if defined(SDC_EVOLVE_ENERGY)
    Real rhoe_old = burn_state.y[SEINT];
#elif defined (SDC_EVOLVE_ENTHALPY)
    Real rhoh_old = burn_state.y[SENTH];
#endif

    // now initialize the advective terms -- the only ones that are
    // actually used during the integration are for rho, (rho X), and
    // (rho e)

    burn_state.ydot_a[SRHO] = Adv_rho;
    burn_state.ydot_a[SEINT] = Adv_rhoe;

    for (int n = 0; n < NumSpec; n++) {
        burn_state.ydot_a[SFS+n] = adv_species[n];
    }

#ifdef AUX_THERMO
    for (int n = 0; n < NumAux; n++) {
        burn_state.ydot_a[SFX+n] = adv_aux[n];
    }
#endif

    // these need to be initialized 

    burn_state.sdc_iter = 1;
    burn_state.num_sdc_iters = 1;

    burn_state.i = 0;
    burn_state.j = 0;
    burn_state.k = 0;

    for (int n = 0; n < NumSpec; ++n) {
        burn_state.xn[n] = eos_state.xn[n];
    }

#ifdef AUX_THERMO
    for (int n = 0; n < NumAux; ++n) {
        burn_state.aux[n] = eos_state.aux[n];
    }
#endif

    burn_state.e = rhoe / burn_state.rho;


    burn_state.T_fixed = -1.0;


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

    burn_t state_in = burn_state;

    // output the data in columns, one line per timestep

    state_over_time << std::setw(10) << "# Time";
    state_over_time << std::setw(12) << "Density";
    state_over_time << std::setw(12) << "Temperature";
    for(int x = 0; x < NumSpec; ++x){
	std::string element = short_spec_names_cxx[x];
	state_over_time << std::setw(12) << element; 
    }
    state_over_time << std::endl;

    Real t = 0.0;

    state_over_time << std::setw(10) << t;
    state_over_time << std::setw(12) << burn_state.rho; 
    state_over_time << std::setw(12) << burn_state.T; 
    for (int x = 0; x < NumSpec; ++x){
        state_over_time << std::setw(12) << burn_state.xn[x]; 
    }
    state_over_time << std::endl;


    // store the initial internal energy -- we'll update this after
    // each substep

    Real energy_initial = burn_state.e;

    // loop over steps, burn, and output the current state

    int nstep_int = 0;

    std::cout << burn_state << std::endl;

    for (int n = 0; n < nsteps; n++){

        // compute the time we wish to integrate to

        Real tend = std::pow(10.0_rt, std::log10(tfirst) + dlogt * n);
        Real dt = tend - t;

        std::cout << "burning for dt = " << dt << std::endl;

    	burner(burn_state, dt);

        nstep_int += burn_state.n_step;

        t += dt;

	state_over_time << std::setw(10) << t;
	state_over_time << std::setw(12) << burn_state.rho;
	state_over_time << std::setw(12) << burn_state.T;
	for (int x = 0; x < NumSpec; ++x){
	     state_over_time << std::setw(12) << burn_state.y[SFS+x] / burn_state.rho; 
	}
	state_over_time << std::endl;
    }
    state_over_time.close();

    std::cout << "successful? " << burn_state.success << std::endl;
    std::cout << "number of steps taken: " << nstep_int << std::endl;

}

#endif
