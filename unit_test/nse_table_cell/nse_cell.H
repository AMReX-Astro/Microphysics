
#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <burner.H>
#include <fstream>
#include <iostream>

#include <nse_table_type.H>
#include <nse_table.H>

using namespace unit_test_rp;

void nse_cell_c()
{


    std::cout << "rho, T, Ye = " << density << " " << temperature << " " << ye << std::endl;

    nse_table_t nse_state;
    nse_state.T = temperature;
    nse_state.rho = density;
    nse_state.Ye = ye;

    nse_interp(nse_state);

    std::cout << "Abar = " << nse_state.abar << std::endl;
    std::cout << "<B/A> = " << nse_state.bea << std::endl;
    std::cout << "dYe/dt = " << nse_state.dyedt << std::endl;
    std::cout << "d<B/A>/dt = " << nse_state.dbeadt << std::endl;
    for (int n = 0; n < NumSpec; ++n) {
        std::cout << "X(" << short_spec_names_cxx[n] << ") = " << nse_state.X[n] << std::endl;
    }

}
