#ifndef NSE_CELL_H
#define NSE_CELL_H

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <burner.H>
#include <fstream>
#include <iostream>

#include <nse_table_type.H>
#include <nse_table.H>

AMREX_INLINE
void nse_cell_c()
{


    std::cout << "rho, T, Ye = "
              << unit_test_rp::density << " "
              << unit_test_rp::temperature << " "
              << unit_test_rp::ye << std::endl;

    nse_table_t nse_state;
    nse_state.T = unit_test_rp::temperature;
    nse_state.rho = unit_test_rp::density;
    nse_state.Ye = unit_test_rp::ye;

    nse_interp(nse_state);

    std::cout << "Abar = " << nse_state.abar << std::endl;
    std::cout << "<B/A> = " << nse_state.bea << std::endl;
    std::cout << "dYe/dt = " << nse_state.dyedt << std::endl;
    std::cout << "d<B/A>/dt = " << nse_state.dbeadt << std::endl;
    for (int n = 0; n < NumSpec; ++n) {
        std::cout << "X(" << short_spec_names_cxx[n] << ") = " << nse_state.X[n] << std::endl;
    }

}
#endif
