#ifndef EOS_CELL_H
#define EOS_CELL_H

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <iostream>
#include <react_util.H>

AMREX_INLINE
void eos_cell_c()
{

    eos_t state;

    // Set mass fractions to sanitize inputs for them
    amrex::Real massfractions[NumSpec];
    for (double &X : massfractions) {
        X = -1.0e0_rt;
    }

    // Make sure user set all the mass fractions to values in the interval [0, 1]
    for (int n = 1; n <= NumSpec; ++n) {

        massfractions[n-1] = get_xn(n, unit_test_rp::uniform_xn);

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
#endif
