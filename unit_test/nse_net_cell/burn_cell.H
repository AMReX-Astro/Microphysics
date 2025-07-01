#ifndef BURN_CELL_H
#define BURN_CELL_H

#include <extern_parameters.H>
#include <eos.H>
#include <network.H>
#include <burner.H>
#include <fstream>
#include <iostream>
#include <nse_solver.H>

AMREX_INLINE
void burn_cell_c()
{

    burn_t state;


    // // Echo initial conditions at burn and fill burn state input

    state.T = unit_test_rp::temperature;
    state.rho = unit_test_rp::density;
    state.y_e = unit_test_rp::ye;

    // set initial chemical potential of proton and neutron

    state.mu_p = unit_test_rp::mu_p;
    state.mu_n = unit_test_rp::mu_n;

    // set a reference cell size.

    state.dx = 1.0e6_rt;

    std::cout << "chemical potential of proton is " << unit_test_rp::mu_p << std::endl;
    std::cout << "chemical potential of neutron is " << unit_test_rp::mu_n << std::endl;

    // find the nse state using (rho, T, Ye)

    bool ye_is_valid = true;
    amrex::Real eps = 1.0e-10_rt;
    auto NSE_STATE = get_actual_nse_state(nse_input_rty, state, eps, ye_is_valid);

    std::cout << "State Density (g/cm^3): " << state.rho << std::endl;
    std::cout << "State Temperature (K): " << state.T << std::endl;
    std::cout << "electron fraction is " << state.y_e << std::endl;

    std::cout << "NSE state: " << std::endl;
    for (int n = 0; n < NumSpec; ++n) {
      std::cout << short_spec_names_cxx[n] << " : " << NSE_STATE.xn[n] << std::endl;
    }

    // Let state.y equal to nse_state.y to make sure its in nse state

    for (int n = 0; n < NumSpec; ++n){
      state.y[SFS+n] = NSE_STATE.y[SFS+n];
    }

    if (in_nse(state)){
      std::cout << "We're in NSE. " << std::endl;
    }
    else{
      std::cout << "We're not in NSE. " << std::endl;
    }
    std::cout << "-----------------------------------------------" << std::endl;

    // now test the EOS+NSE function.  The basic idea is the following:
    //
    // 1. find the e corresponding to the current NSE state (eos_input_rt)
    // 2. perturb e
    // 3. call the NSE function using input (rho, e, Ye) to get the X's
    //    and T that is consistent with e.
    // 4. call the EOS with this new T (eos_input_rt) and updated NSE
    //    composition and see if we recover the same perturbed e

    eos_t eos_state;
    eos_state.T = state.T;
    eos_state.rho = state.rho;
    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = state.y[SFS+n] / state.rho;
    }

    // get the initial e and abar corresponding to this T

    eos(eos_input_rt, eos_state);

    // now perturb e and find NSE using (rho, e, Ye)

    amrex::Real e_start = eos_state.e * 1.05;
    amrex::Real T_start{eos_state.T};
    amrex::Real abar_start{eos_state.abar};

    // Initialize the internal energy

    state.e = e_start;

    // Find the new NSE state using input (rho, e, Ye)

    // We will test results from two modes of solving for T using e.

    // First using mode == 1
    {
    solve_nse_e_mode = 1;
    auto nse_state = get_actual_nse_state(nse_input_rey, state, 1.e-10, true);

    // now try calling the EOS with updated T and the NSE composition and see if
    // we get back the same energy and see how abar changed.

    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = nse_state.xn[n];
    }

    eos_state.T = nse_state.T;

    eos(eos_input_rt, eos_state);

    // Print to see the change in T and abar due to perturbed e
    // and see if we recover the same energy using eos.

    std::cout << "Using nse.solve_nse_e_mode == 1:" << std::endl;
    std::cout << "change in T: " << T_start << " " << eos_state.T <<  std::endl;
    std::cout << "change in abar: " << abar_start << " " << eos_state.abar <<  std::endl;
    std::cout << "recovered energy: " << e_start << " " << eos_state.e << std::endl;

    // See the result of constraint eq and jacobian using the solved input

    amrex::Array1D<amrex::Real, 1, 3> x {nse_state.mu_p, nse_state.mu_n, nse_state.T};
    amrex::Array1D<amrex::Real, 1, 3> f;
    amrex::Array2D<amrex::Real, 1, 3, 1, 3> jac;
    int flag = 0;

    nse_solver_data<burn_t> state_data = {nse_state, {0.0_rt}};

    nse_fcn(x, f, state_data, flag);
    nse_jcn(x, jac, state_data, flag);

    std::cout << "Constraint Eqns using solved input are: ";
    for (int q = 1; q <= 3; ++q) {
        std::cout << f(q) << " ";
    }
    std::cout << std::endl;

    std::cout << "Constraint Jacobian using solved input are: ";
    for (int q = 1; q <= 3; ++q) {
        for (int p = 1; p <= 3; ++p) {
            std::cout << jac(q, p) << " ";
        }
    }
    std::cout << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;
    }

    // Now do the same but with mode == 2
    {
    solve_nse_e_mode = 2;
    auto nse_state = get_actual_nse_state(nse_input_rey, state, 1.e-10, true);

    // now try calling the EOS with updated T and the NSE composition and see if
    // we get back the same energy and see how abar changed.

    for (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = nse_state.xn[n];
    }

    eos_state.T = nse_state.T;

    eos(eos_input_rt, eos_state);

    // Print to see the change in T and abar due to perturbed e
    // and see if we recover the same energy using eos.

    std::cout << "Using nse.solve_nse_e_mode == 2:" << std::endl;
    std::cout << "change in T: " << T_start << " " << eos_state.T <<  std::endl;
    std::cout << "change in abar: " << abar_start << " " << eos_state.abar <<  std::endl;
    std::cout << "recovered energy: " << e_start << " " << eos_state.e << std::endl;

    // See the result of constraint eq and jacobian using the solved input

    amrex::Array1D<amrex::Real, 1, 3> x {nse_state.mu_p, nse_state.mu_n, nse_state.T};
    amrex::Array1D<amrex::Real, 1, 3> f;
    amrex::Array2D<amrex::Real, 1, 3, 1, 3> jac;
    int flag = 0;

    nse_solver_data<burn_t> state_data = {nse_state, {0.0_rt}};

    nse_fcn(x, f, state_data, flag);
    nse_jcn(x, jac, state_data, flag);

    std::cout << "Constraint Eqns using solved input are: ";
    for (int q = 1; q <= 3; ++q) {
        std::cout << f(q) << " ";
    }
    std::cout << std::endl;

    std::cout << "Constraint Jacobian using solved input are: ";
    for (int q = 1; q <= 3; ++q) {
        for (int p = 1; p <= 3; ++p) {
            std::cout << jac(q, p) << " ";
        }
    }
    std::cout << std::endl;
    }
}
#endif
