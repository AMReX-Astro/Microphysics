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
#include <numerical_jacobian.H>

using namespace unit_test_rp;

AMREX_INLINE
void jac_cell_c()
{

    amrex::Real massfractions[NumSpec] = {-1.0};

    // Make sure user set all the mass fractions to values in the interval [0, 1]
    for (int n = 1; n <= NumSpec; ++n) {

        massfractions[n-1] = get_xn(n, unit_test_rp::uniform_xn);

        if (massfractions[n-1] < 0 || massfractions[n-1] > 1) {
            amrex::Error("mass fraction for " + short_spec_names_cxx[n-1] + " not initialized in the interval [0,1]!");
        }

    }

    burn_t burn_state;

    burn_state.rho = density;
    burn_state.T = temperature;
    for (int n = 0; n < NumSpec; ++n) {
        burn_state.xn[n] = massfractions[n];
    }

    if (! unit_test_rp::skip_initial_normalization) {
        normalize_abundances_burn(burn_state);
    }

    // get the energy

    eos(eos_input_rt, burn_state);

    JacNetArray2D jac_analytic;
    JacNetArray2D jac_numerical;

    actual_jac(burn_state, jac_analytic);

    // we need to convert this to be in terms of X instead of Y

    for (int j = 1; j <= NumSpec; ++j) {
        for (int i = 1; i <= neqs; ++i) {
            jac_analytic(j, i) *= aion[j-1];
            jac_analytic(i, j) *= aion_inv[j-1];
        }
    }

    jac_info_t jac_info;
    jac_info.h = 1.e-5_rt;  // timestep really doesn't make sense here

    numerical_jac(burn_state, jac_info, jac_numerical);

    // output

    for (int ii = 1; ii <= neqs; ++ii) {
        std::string ilabel = (ii < neqs) ? short_spec_names_cxx[ii-1] : "e";
        for (int jj = 1; jj <= neqs; ++jj) {
            std::string jlabel = (jj < neqs) ? short_spec_names_cxx[jj-1] : "e";
            std::cout << "J(" << std::setw(4) << ilabel << ", " << std::setw(4) << jlabel << ") = "
                      << std::setw(20) << jac_numerical(ii,jj) << " "
                      << std::setw(20) << jac_analytic(ii,jj) << std::endl;
        }
        std::cout << std::endl;
    }


}
#endif
