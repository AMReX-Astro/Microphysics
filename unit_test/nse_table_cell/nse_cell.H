
#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <burner.H>
#include <fstream>
#include <iostream>

#include <nse_table.H>

using namespace unit_test_rp;

void nse_cell_c()
{


    std::cout << "rho, T, Ye = " << density << " " << temperature << " " << ye << std::endl;

    Real abar{};
    Real dq{};
    Real dyedt{};
    Real X[NumSpec] = {};

    nse_interp(density, temperature, ye,
               abar, dq, dyedt, X);

    std::cout << "abar = " << abar << std::endl;
    std::cout << "B/A = " << dq << std::endl;
    std::cout << "dYe/dt = " << dyedt << std::endl;
    for (int n = 0; n < NumSpec; ++n) {
        std::cout << "X(" << short_spec_names_cxx[n] << ") = " << X[n] << std::endl;
    }

}
