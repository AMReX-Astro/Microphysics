
#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <fstream>
#include <iostream>

void eos_cell_c()
{

    eos_t state;

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

    state.T   = temperature;
    state.rho = density;
    for (int n = 0; n < NumSpec; ++n) {
        state.xn[n] = massfractions[n];
    }

    eos(eos_input_rt, state);

    std::cout << state;

}
