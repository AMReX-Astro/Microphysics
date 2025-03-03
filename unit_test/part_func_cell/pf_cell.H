#ifndef PF_CELL_H
#define PF_CELL_H

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <burner.H>
#include <fstream>
#include <iostream>

#include <reaclib_rates.H>
#include <partition_functions.H>

AMREX_INLINE
void pf_cell_c()
{

    amrex::Real pf;
    amrex::Real dpf_dT;

    std::cout << "temperature = " << unit_test_rp::temperature << std::endl;

    tf_t tfactors = evaluate_tfactors(unit_test_rp::temperature);

    get_partition_function(Ni56, tfactors, pf, dpf_dT);

    std::cout << "Ni56: " << pf << " " << dpf_dT << std::endl;

    get_partition_function(Fe52, tfactors, pf, dpf_dT);

    std::cout << "Fe52: " << pf << " " << dpf_dT << std::endl;

    std::cout << "spins: " << get_spin_state(H1) << " " << get_spin_state(He4) << " " << get_spin_state(Fe52) << " " << get_spin_state(Co55) << " " << get_spin_state(Ni56) << std::endl;

}
#endif
