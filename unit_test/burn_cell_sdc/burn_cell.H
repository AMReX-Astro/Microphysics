
#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <burner.H>
#include <fstream>
#include <iostream>

void burn_cell_c()
{

    Real massfractions[NumSpec] = {-1.0};

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


#ifdef NSE_THERMO
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

        }

    }


#ifdef NSE_THERMO
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
#ifdef NSE_THERMO
    set_nse_aux_from_X(eos_state);
#endif

    eos(eos_input_rt, eos_state);

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
    std::cout << "density set" << std::endl;

    burn_state.y[SRHO] = density;

    // we will pick velocities to be 10% of the sound speed
    burn_state.y[SMX] = 0.0;
    burn_state.y[SMY] = 0.0;
    burn_state.y[SMZ] = 0.0;

    burn_state.y[SEINT] = burn_state.y[SRHO] * eos_state.e;
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

    for (int n = 0; n < NumSpec; n++) {
        burn_state.ydot_a[SFS+n] = adv_species[n];
    }

    // these need to be initialized 

    burn_state.sdc_iter = 1;
    burn_state.num_sdc_iters = 1;

    for (int n = 0; n < NumSpec; ++n) {
        burn_state.xn[n] = massfractions[n];
    }

    burn_state.e = eos_state.e;


    // output initial burn type data
    Real time = 0.0_rt;


    std::ofstream state_over_time("state_over_time.txt");

    // we will divide the total integration time into nsteps that are
    // logarithmically spaced

    Real dlogt = 0.0_rt;

    if (tfirst == 0.0_rt) {
        dlogt = std::log10(tmax) / nsteps;
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
    state_over_time << std::setw(12) << burn_state.T; 
    for (int x = 0; x < NumSpec; ++x){
        state_over_time << std::setw(12) << burn_state.xn[x]; 
    }
    state_over_time << std::endl;


    // store the initial internal energy -- we'll update this after
    // each substep

    Real energy_initial = burn_state.e;

    // loop over steps, burn, and output the current state

    for (int n = 0; n < nsteps; n++){

        // compute the time we wish to integrate to

        Real tend;
        if (tfirst == 0.0_rt) {
            tend = std::pow(10.0_rt, dlogt * (n + 1));
        } else {
            tend = std::pow(10.0_rt, std::log10(tfirst) + dlogt * n);
        }

        Real dt = tend - t;

    	burner(burn_state, dt);

        // burn_state.e represents the change in energy over the burn (for
        // just this sybcycle), so turn it back into a physical energy

        burn_state.e += energy_initial;

        // reset the initial energy for the next subcycle

        energy_initial = burn_state.e;

        // get the updated T

        eos(eos_input_re, burn_state);

        t += dt;

	state_over_time << std::setw(10) << t;
	state_over_time << std::setw(12) << burn_state.T;
	state_over_time << std::setw(12) << burn_state.rho;
	for (int x = 0; x < NumSpec; ++x){
	     state_over_time << std::setw(12) << burn_state.y[SFS+x] / burn_state.rho; 
	}
	state_over_time << std::endl;
    }
    state_over_time.close();

    std::cout << "successful? " << burn_state.success << std::endl;

}
