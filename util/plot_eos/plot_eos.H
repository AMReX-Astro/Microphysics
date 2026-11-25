#ifndef EOS_CELL_H
#define EOS_CELL_H

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <iostream>
#include <react_util.H>

AMREX_INLINE
void output_plot_info()
{

    eos_t state;

    // set the mass fractions to be equal

    amrex::Real massfractions[NumSpec];
    for (double &X : state.xn) {
        X = 1.0_rt / static_cast<amrex::Real>(NumSpec);
    }

    // compute spacing between points

    amrex::Real dlogT = (std::log10(plot_rp::T_max) -
                         std::log10(plot_rp::T_min)) /
        static_cast<amrex::Real>(plot_rp::nT - 1);

    amrex::Real dlogrho = (std::log10(plot_rp::rho_max) -
                           std::log10(plot_rp::rho_min)) /
        static_cast<amrex::Real>(plot_rp::nrho - 1);

    // loop

    for (int itemp = 0; itemp <= plot_rp::nT; ++itemp) {
        state.T = std::pow(10.0_rt, std::log10(plot_rp::T_min) + itemp * dlogT);

        for (int jrho = 0; jrho <= plot_rp::nrho; ++jrho) {
            state.rho = std::pow(10.0_rt, std::log10(plot_rp::rho_min) + jrho * dlogrho);

            eos(eos_input_rt, state);

            std::cout << state;

        }
    }

}
#endif
