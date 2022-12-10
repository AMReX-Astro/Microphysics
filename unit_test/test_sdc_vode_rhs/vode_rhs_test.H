#ifndef VODE_RHS_TEST_H
#define VODE_RHS_TEST_H

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#ifdef NEW_NETWORK_IMPLEMENTATION
#include <rhs.H>
#endif
#include <burner.H>
#include <fstream>
#include <iostream>
#include <integrator_rhs_simplified_sdc.H>
#include <react_util.H>

void do_vode_rhs()
{

    Real massfractions[NumSpec] = {-1.0};

    for (int n = 1; n <= NumSpec; ++n) {

        massfractions[n-1] = get_xn(n);

        if (massfractions[n-1] < 0 || massfractions[n-1] > 1) {
            amrex::Error("mass fraction for " + short_spec_names_cxx[n-1] + " not initialized in the interval [0,1]!");
        }

    }


#ifdef AUX_THERMO
    Real auxdata[NumAux] = {-1.0};

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

        }

    }


#ifdef AUX_THERMO
    Real adv_aux[NumAux] = {0.0};

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
#ifdef AUX_THERMO
    set_aux_comp_from_X(eos_state);
#endif

    eos(eos_input_rt, eos_state);

    burn_state.rho = eos_state.rho;
    burn_state.T = eos_state.T;

    for (int n = 0; n < NumSpec; n++) {
      burn_state.y[SFS+n] = burn_state.rho * eos_state.xn[n];
    }

    burn_state.T_fixed = -1.e200_rt;

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
        burn_state.xn[n] = massfractions[n];
    }

#ifdef AUX_THERMO
    set_aux_comp_from_X(burn_state);
#endif

    // now make the dvode_t

    dvode_t vode_state;

    vode_state.t = 0.0_rt;

    // initialize the EOS state using the conserved state that was
    // passed in

    eos(eos_input_rt, burn_state);

    // Fill in the initial integration state.

    vode_state.t = 0.0_rt;
    vode_state.tn = 0.0_rt;   // this is needed for directly accessing the jac() function
    vode_state.H = tmax/1000.0_rt;  // needed for jacobian = 3 scaling
    vode_state.tout = tmax;
    vode_state.NFE = 0;
    vode_state.NJE = 0;

    burn_to_int(burn_state, vode_state);

    // call the RHS routine

    std::cout << "just evaluate the RHS" << std::endl;

    RArray1D ydot;
    rhs(0.0_rt, burn_state, vode_state, ydot);

    std::cout << "ydot = " << std::endl;
    for (int n = 1; n <= VODE_NEQS; n++) {
        std::cout << n << " " << ydot(n) << std::endl;
    }

    // call the Jacobian

    RArray2D pd;
    const Real time = 0.0_rt;
    jac(time, burn_state, vode_state, pd);

    for (int row = 1; row <= VODE_NEQS; ++row) {
        for (int col = 1; col <= VODE_NEQS; ++col) {
            std::cout << "(" << std::setw(2) << row << "," << std::setw(2) << col << ") = " << pd(row,col) << std::endl;
        }
    }

    // do a full burn to make sure this zone can be integrated

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
