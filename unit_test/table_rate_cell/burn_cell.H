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

#include <network_properties.H>
#include <table_rates.H>

AMREX_INLINE
void burn_cell_c()
{

    amrex::Real massfractions[NumSpec] = {1.0_rt / NumSpec};


    std::cout << std::endl;
    std::cout << "State Density (g/cm^3): " << unit_test_rp::density << std::endl;
    std::cout << "State Temperature (K): " << unit_test_rp::temperature << std::endl;
    std::cout << std::endl;

    amrex::Real ye{};
    for (int i = 0; i < NumSpec; ++i) {
        ye += zion[i] * massfractions[i] / aion[i];
    }

    auto rhoy = unit_test_rp::density * ye;
    auto T = unit_test_rp::temperature;

    amrex::Real log_temp = std::log10(T);
    amrex::Real log_rhoy = std::log10(rhoy);

    [[maybe_unused]] amrex::Real rate, drate_dt, edot_nu, edot_gamma;

    tabular_evaluate(j_Co56_Fe56_meta, j_Co56_Fe56_rhoy, j_Co56_Fe56_temp, j_Co56_Fe56_data,
                     log_rhoy, log_temp, T, rate, drate_dt, edot_nu, edot_gamma);

    std::cout << "Co56 -> Fe56 rate = " << rate << std::endl;

    tabular_evaluate(j_Co56_Ni56_meta, j_Co56_Ni56_rhoy, j_Co56_Ni56_temp, j_Co56_Ni56_data,
                     log_rhoy, log_temp, T, rate, drate_dt, edot_nu, edot_gamma);

    std::cout << "Co56 -> Ni56 rate = " << rate << std::endl;

}
#endif
