#ifndef VODE_RHS_TEST_H
#define VODE_RHS_TEST_H

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <burner.H>
#include <fstream>
#include <iostream>
#include <vode_rhs_simplified_sdc.H>

void do_vode_rhs()
{

    eos_t eos_state;

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
    std::cout << "State Density (g/cm^3): " << density << std::endl;
    std::cout << "State Temperature (K): " << temperature << std::endl;
    for (int n = 0; n < NumSpec; ++n) {
        std::cout << "Mass Fraction (" << short_spec_names_cxx[n] << "): " << massfractions[n] << std::endl;
    }

    eos_state.T   = temperature;
    eos_state.rho = density;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = massfractions[n];
    }

#ifdef NSE_THERMO
    set_nse_aux_from_X(eos_state);
#endif

    eos(eos_input_rt, eos_state);

    burn_t burn_state;

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
    burn_state.y[SRHO] = eos_state.rho;

    // we will pick velocities to be 10% of the sound speed
    burn_state.y[SMX] = burn_state.y[SRHO] * 0.1_rt * eos_state.cs;
    burn_state.y[SMY] = burn_state.y[SRHO] * 0.1_rt * eos_state.cs;
    burn_state.y[SMZ] = burn_state.y[SRHO] * 0.1_rt * eos_state.cs;

    burn_state.y[SEINT] = burn_state.y[SRHO] * eos_state.e;
    burn_state.y[SEDEN] = burn_state.y[SEINT] +
      0.5_rt * (burn_state.y[SMX] * burn_state.y[SMX] +
                burn_state.y[SMY] * burn_state.y[SMY] +
                burn_state.y[SMZ] * burn_state.y[SMZ]) / burn_state.y[SRHO];

    // get T from total energy
    burn_state.T_from_eden = true;

#elif defined(SDC_EVOLVE_ENTHALPY)
    burn_state.y[SENTH] = burn_state.rho * eos_state.h;

    burn_state.p0 = eos_state.p;
#endif

    // normalize
    Real sum_spec = 0.0_rt;
    for (int n = 0; n < NumSpec; n++) {
      sum_spec += burn_state.y[SFS+n] / burn_state.rho;
    }
    for (int n = 0; n < NumSpec; n++) {
      burn_state.y[SFS+n] /= sum_spec;
    }

#if defined(SDC_EVOLVE_ENERGY)
    Real rhoe_old = burn_state.y[SEINT];
#elif defined (SDC_EVOLVE_ENTHALPY)
    Real rhoh_old = burn_state.y[SENTH];
#endif

    // zero out the advective terms

    for (int n = 0; n < SVAR; n++) {
        burn_state.ydot_a[n] = 0.0;
    }

    // these need to be initialized 

    burn_state.sdc_iter = 1;
    burn_state.num_sdc_iters = 1;


    // now make the dvode_t

    dvode_t vode_state;

    vode_state.t = 0.0_rt;

    // initialize the EOS state using the conserved state that was
    // passed in

    sdc_burn_to_eos(burn_state, eos_state);

    eos(eos_input_rt, eos_state);

    eos_to_burn(eos_state, burn_state);

    // Fill in the initial integration state.

    burn_to_vode(burn_state, vode_state);

    // call the RHS routine

    std::cout << "just evaluate the RHS" << std::endl;

    RArray1D ydot;
    rhs(0.0_rt, burn_state, vode_state, ydot);

    std::cout << "ydot = " << std::endl;
    for (int n = 1; n <= NumSpec+2; n++) {
        std::cout << n << " " << ydot(n) << std::endl;
    }
 
    std::cout << "now do a burn" << std::endl;

    integrator(burn_state, tmax);

    for (int n = 0; n < SVAR; n++) {
        if (n == SFS) {
            std::cout << "species: " << std::endl;
#if NAUX_NET > 0
        } else if (n == SFX) {
            std::cout << "aux: " << std::endl;
#endif
        }
        std::cout << n << " " << burn_state.y[n] << std::endl;
    }

    std::cout << "number of RHS calls = " << burn_state.n_rhs << std::endl;

}
#endif
