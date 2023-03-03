#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <burner.H>
#include <fstream>
#include <iostream>

Real grav_constant = 6.674e-8;

void burn_cell_c()
{

    burn_t state;

    Real numdens[NumSpec] = {-1.0};

    for (int n = 1; n <= NumSpec; ++n) {
        switch (n) {

        case 1:
            numdens[n-1] = primary_species_1;
            break;
        case 2:
            numdens[n-1] = primary_species_2;
            break;
        case 3:
            numdens[n-1] = primary_species_3;
            break;
        case 4:
            numdens[n-1] = primary_species_4;
            break;
        case 5:
            numdens[n-1] = primary_species_5;
            break;
        case 6:
            numdens[n-1] = primary_species_6;
            break;
        case 7:
            numdens[n-1] = primary_species_7;
            break;
        case 8:
            numdens[n-1] = primary_species_8;
            break;
        case 9:
            numdens[n-1] = primary_species_9;
            break;
        case 10:
            numdens[n-1] = primary_species_10;
            break;
        case 11:
            numdens[n-1] = primary_species_11;
            break;
        case 12:
            numdens[n-1] = primary_species_12;
            break;
        case 13:
            numdens[n-1] = primary_species_13;
            break;
        case 14:
            numdens[n-1] = primary_species_14;
            break;

        }

    }

    // Echo initial conditions at burn and fill burn state input

    std::cout << "Maximum Time (s): " << tmax << std::endl;
    std::cout << "State Temperature (K): " << temperature << std::endl;
    for (int n = 0; n < NumSpec; ++n) {
        std::cout << "Number Density input (" << short_spec_names_cxx[n] << "): " << numdens[n] << std::endl;
    }

    state.T = temperature;

    //find the density in g/cm^3
    Real rhotot = 0.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
        state.xn[n] = numdens[n];
        rhotot += state.xn[n]*spmasses[n];  //spmasses contains the masses of all species, defined in EOS
    }

    state.rho = rhotot;

    // normalize -- just in case

    Real mfracs[NumSpec] = {-1.0};
    Real msum = 0.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
        mfracs[n] = state.xn[n]*spmasses[n]/rhotot;
        msum += mfracs[n];
    }

    for (int n = 0; n < NumSpec; ++n) {
        mfracs[n] /= msum;
        //use the normalized mass fractions to obtain the corresponding number densities
        state.xn[n] = mfracs[n]*rhotot/spmasses[n];
    }

#ifdef DEBUG
    for (int n = 0; n < NumSpec; ++n) {
           std::cout << "Mass fractions (" << short_spec_names_cxx[n] << "): " << mfracs[n] << std::endl;
           std::cout << "Number Density conserved (" << short_spec_names_cxx[n] << "): " << state.xn[n] << std::endl;
    }
#endif


    // call the EOS to set initial internal energy e
    eos(eos_input_rt, state);

    // name of output file
    std::ofstream state_over_time("state_over_time.txt");

    // save the initial state -- we'll use this to determine
    // how much things changed over the entire burn
    burn_t state_in = state;

    // output the data in columns, one line per timestep
    state_over_time << std::setw(10) << "# Time";
    state_over_time << std::setw(15) << "Density";
    state_over_time << std::setw(15) << "Temperature";
    for(int x = 0; x < NumSpec; ++x){
    std::string element = short_spec_names_cxx[x];
    state_over_time << std::setw(15) << element;
    }
    state_over_time << std::endl;

    Real t = 0.0;

    state_over_time << std::setw(10) << t;
    state_over_time << std::setw(15) << state.rho;
    state_over_time << std::setw(15) << state.T;
    for (int x = 0; x < NumSpec; ++x){
        state_over_time << std::setw(15) << state.xn[x];
    }
    state_over_time << std::endl;


    // loop over steps, burn, and output the current state
    // the loop below is similar to that used in krome and pynucastro
    Real dd = rhotot;
    Real dd1 = 0.0_rt;

    for (int n = 0; n < nsteps; n++){

        dd1 = dd;

        Real rhotmp = 0.0_rt;

        for (int nn = 0; nn < NumSpec; ++nn) {
            rhotmp += state.xn[nn]*spmasses[nn];
        }

        // find the freefall time
        Real tff = std::sqrt(M_PI*3.0 / (32.0*rhotmp*grav_constant));
        Real dt = tff_reduc * tff;
        // scale the density
        dd += dt*(dd/tff);

#ifdef DEBUG
        std::cout<<"step params "<<dd<<", "<<tff<<", "<<rhotmp<<std::endl;
#endif

        // stop the test if dt is very small
        if (dt < 10) {
        break; }

        // stop the test if we have reached very high densities
        if (dd > 3e-6) {
        break; }

        std::cout<<"step "<<n<<" done"<<std::endl;

        // scale the number densities
        for (int nn = 0; nn < NumSpec; ++nn) {
            state.xn[nn] *= dd/dd1;
        }

        // input the scaled density in burn state
        state.rho *= dd/dd1;

        // do the actual integration
        burner(state, dt);

       // ensure positivity and normalize
        Real inmfracs[NumSpec] = {-1.0};
        Real insum = 0.0_rt;
        for (int nn = 0; nn < NumSpec; ++nn) {
            state.xn[nn] = amrex::max(state.xn[nn], small_x);
            inmfracs[nn] = spmasses[nn]*state.xn[nn]/state.rho;
            insum += inmfracs[nn];
        }

        for (int nn = 0; nn < NumSpec; ++nn) {
            inmfracs[nn] /= insum;
            //update the number densities with conserved mass fractions
            state.xn[nn] = inmfracs[nn]*state.rho/spmasses[nn];
        }

       //update the number density of electrons due to charge conservation
       state.xn[0] = -state.xn[3] - state.xn[7] + state.xn[1] + state.xn[12] + state.xn[6] + state.xn[4] + state.xn[9] + 2.0*state.xn[11];

       //reconserve mass fractions post charge conservation
       insum = 0;
        for (int nn = 0; nn < NumSpec; ++nn) {
            state.xn[nn] = amrex::max(state.xn[nn], small_x);
            inmfracs[nn] = spmasses[nn]*state.xn[nn]/state.rho;
            insum += inmfracs[nn];
        }

        for (int nn = 0; nn < NumSpec; ++nn) {
            inmfracs[nn] /= insum;
            //update the number densities with conserved mass fractions
            state.xn[nn] = inmfracs[nn]*state.rho/spmasses[nn];
        }

        // get the updated T
        eos(eos_input_re, state);

        t += dt;

    state_over_time << std::setw(10) << t;
        state_over_time << std::setw(15) << state.rho;
    state_over_time << std::setw(12) << state.T; 
    for (int x = 0; x < NumSpec; ++x){
         state_over_time << std::setw(15) << state.xn[x]; 
    }
    state_over_time << std::endl;


    }
    state_over_time.close();

    // output diagnostics to the terminal
    std::cout << "------------------------------------" << std::endl;
    std::cout << "successful? " << state.success << std::endl;

    std::cout << "------------------------------------" << std::endl;
    std::cout << "T initial = " << state_in.T << std::endl;
    std::cout << "T final =   " << state.T << std::endl;
    std::cout << "Eint initial = " << state_in.e << std::endl;
    std::cout << "Eint final =   " << state.e << std::endl;
    std::cout << "rho initial = " << state_in.rho << std::endl;
    std::cout << "rho final =   " << state.rho << std::endl;

    std::cout << "------------------------------------" << std::endl;
    std::cout << "New number densities: " << std::endl;
    for (int n = 0; n < NumSpec; ++n) {
        std::cout << state.xn[n] << std::endl;
    }

}
