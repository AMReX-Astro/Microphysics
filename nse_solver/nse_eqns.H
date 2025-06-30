#ifndef NSE_EQNS_H
#define NSE_EQNS_H

#include <AMReX_REAL.H>
#include <AMReX_Algorithm.H>
#include <nse_solver.H>
#include <extern_parameters.H>

///
/// NSE constraint equations:
///
/// 1) Mass Conservation: Sum_k X_k - 1 = 0
/// 2) Charge Conservation: Sum_k Z_k X_k / A_k - Ye = 0
/// 3) Energy Conservation: e(rho, X_k^NSE, Ye) - e_in = 0
///
/// When the input is (rho, T, Ye), we solve for the first two Eqns
/// When the input is (rho, e, Ye), we solve for all three Eqns.
/// the energy constraint is used to solve for Temperature.
///

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_fcn(amrex::Array1D<amrex::Real, 1, 2>& x,
             amrex::Array1D<amrex::Real, 1, 2>& fvec,
             const T& state_data, int& iflag) {
    // Constraint Equations for (rho, T, Ye) input.

    // here state is the nse_state from get_nonexponent_nse_state

    amrex::ignore_unused(iflag);

    auto nse_state = state_data.state;
    nse_state.mu_p = x(1);
    nse_state.mu_n = x(2);

    // Apply exponent component for calculating nse mass fractions

    apply_nse_exponent(nse_state, state_data.u_c);

    fvec(1) = -1.0_rt;

    for (int n = 0; n < NumSpec; ++n) {
#ifdef NEW_NETWORK_IMPLEMENTATION
        if (n == NSE_INDEX::H1_index) {
            continue;
        }
#endif
        // constraint equation 1, mass fraction sum to 1

        fvec(1) += nse_state.xn[n];
    }

    // constraint equation 2, electron fraction should be the same

    fvec(2) = nse_state.y_e - state_data.state.y_e;

}


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_fcn(amrex::Array1D<amrex::Real, 1, 3>& x,
             amrex::Array1D<amrex::Real, 1, 3>& fvec,
             const T& state_data, int& iflag) {
    // Constraint Equations for (rho, e, Ye) input.

    amrex::ignore_unused(iflag);

    // Initialize burn state

    burn_t burn_state;
    burn_state.rho = state_data.state.rho;
    burn_state.mu_p = x(1);
    burn_state.mu_n = x(2);
    burn_state.T    = x(3);

    amrex::Real e_in = state_data.state.e;

    // Find NSE X_k

    // Find nonexponent piece, note that nse_state is eos_t
    // so that we can do eos on it.

    auto nse_state = get_nonexponent_nse_state<burn_t>(burn_state);

    // Find screening term

#if SCREEN_METHOD == SCREEN_METHOD_chabrier1998
    compute_coulomb_contribution(state_data.u_c, nse_state);
#endif

    // Apply exponent for calculating NSE mass fractions

    apply_nse_exponent(nse_state, state_data.u_c);

    fvec(1) = -1.0_rt;

    for (int n = 0; n < NumSpec; ++n) {
#ifdef NEW_NETWORK_IMPLEMENTATION
        if (n == NSE_INDEX::H1_index) {
            continue;
        }
#endif
        // constraint equation 1, mass fraction sum to 1

        fvec(1) += nse_state.xn[n];
    }

    // constraint equation 2, electron fraction should be the same
    // note that nse_state.y_e is already updated in apply_nse_exponent()

    fvec(2) = nse_state.y_e - state_data.state.y_e;

    // constraint equation 3, e(rho, T, X^NSE) - e_in = 0

    // Create eos state to find the corresponding internal energy

    eos_extra_t eos_state;
    for  (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = nse_state.xn[n];
    }
    eos_state.rho = nse_state.rho;
    eos_state.T = nse_state.T;
    eos(eos_input_rt, eos_state);

    fvec(3) = eos_state.e - e_in;
}

///
/// NSE constraint Jacobian: ∂F_i/∂x_j
///
/// When the input is (rho, T, Ye), x_j = {μ_p, μ_n}
/// When the input is (rho, e, Ye), x_j = {μ_p, μ_n, T}
///

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_jcn(amrex::Array1D<amrex::Real, 1, 2>& x,
             amrex::Array2D<amrex::Real, 1, 2, 1, 2>& fjac,
             const T& state_data, int& iflag) {
    // 2x2 Constraint Jacobians for (rho, T, Ye) input.

    // here state is the nse_state from get_nonexponent_nse_state

    amrex::ignore_unused(iflag);

    auto nse_state = state_data.state;
    nse_state.mu_p = x(1);
    nse_state.mu_n = x(2);

    apply_nse_exponent(nse_state, state_data.u_c);

    // evaluate jacobian of the constraint

    fjac(1, 1) = 0.0_rt;
    fjac(1, 2) = 0.0_rt;
    fjac(2, 1) = 0.0_rt;
    fjac(2, 2) = 0.0_rt;

    // if we are doing drive_initial_convection, we want to use
    // the temperature that comes in through T_fixed

    amrex::Real T_in = nse_state.T_fixed > 0.0_rt ? nse_state.T_fixed : nse_state.T;
    amrex::Real ikTMeV = C::MeV2erg / (C::k_B * T_in);

    for (int n = 0; n < NumSpec; ++n) {
#ifdef NEW_NETWORK_IMPLEMENTATION
        if (n == NSE_INDEX::H1_index) {
            continue;
        }
#endif
        fjac(1, 1) += nse_state.xn[n] * zion[n] * ikTMeV;
        fjac(1, 2) += nse_state.xn[n] * (aion[n] - zion[n]) * ikTMeV;
        fjac(2, 1) += nse_state.xn[n] * zion[n] * zion[n] * aion_inv[n] * ikTMeV;
        fjac(2, 2) += nse_state.xn[n] * zion[n] * (aion[n] - zion[n]) * aion_inv[n] * ikTMeV;
    }

}


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_jcn(amrex::Array1D<amrex::Real, 1, 3>& x,
             amrex::Array2D<amrex::Real, 1, 3, 1, 3>& fjac,
             const T& state_data, int& iflag) {
    // 3x3 Constraint Jacobians for (rho, e, Ye) input.

    amrex::ignore_unused(iflag);

    // Initialize burn state

    burn_t burn_state;
    burn_state.rho = state_data.state.rho;
    burn_state.mu_p = x(1);
    burn_state.mu_n = x(2);
    burn_state.T    = x(3);

    // Find NSE X_k

    // Find nonexponent piece.

    auto nse_state = get_nonexponent_nse_state(burn_state);

    // Find screening term

#if SCREEN_METHOD == SCREEN_METHOD_chabrier1998
    compute_coulomb_contribution(state_data.u_c, nse_state);
#endif

    // Apply exponent for calculating NSE mass fractions

    apply_nse_exponent(nse_state, state_data.u_c);

    // evaluate jacobian of the constraint

    for (int i = 1; i <= 3; ++i) {
        for (int j = 1; j <= 3; ++j) {
            fjac(i, j) = 0.0_rt;
        }
    }

    // if we are doing drive_initial_convection, we want to use
    // the temperature that comes in through T_fixed

    amrex::Real T_in = nse_state.T_fixed > 0.0_rt ? nse_state.T_fixed : nse_state.T;
    amrex::Real ikTMeV = C::MeV2erg / (C::k_B * T_in);

    // Create eos state to find relevant dedZ, dedA and abar via eos

    eos_extra_t eos_state;
    for  (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = nse_state.xn[n];
    }
    eos_state.rho = nse_state.rho;
    eos_state.T = nse_state.T;
    eos(eos_input_rt, eos_state);

    for (int n = 0; n < NumSpec; ++n) {
#ifdef NEW_NETWORK_IMPLEMENTATION
        if (n == NSE_INDEX::H1_index) {
            continue;
        }
#endif
        // Pre-compute some relevant derivatives

        amrex::Real dXdmu_p = nse_state.xn[n] * zion[n] * ikTMeV;
        amrex::Real dXdmu_n = nse_state.xn[n] * (aion[n] - zion[n]) * ikTMeV;

        // fill partition function

        amrex::Real pf = 1.0_rt;
        amrex::Real dpf_dT = 0.0_rt;

#ifndef NEW_NETWORK_IMPLEMENTATION
        get_partition_function(n+1, tfactors, pf, dpf_dT);
#endif

        // There is another screening term I need to add if we have screening.
        amrex::Real dXdT = nse_state.xn[n] *
            (dpf_dT / pf + (1.5_rt - zion[n] * nse_state.mu_p - (aion[n] - zion[n]) * nse_state.mu_n -
                            network::bion(n+1)) / nse_state.T);

        fjac(1, 1) += dXdmu_p;
        fjac(1, 2) += dXdmu_n;
        fjac(1, 3) += dXdT;
        fjac(2, 1) += zion[n] * aion_inv[n] * dXdmu_p;
        fjac(2, 2) += zion[n] * aion_inv[n] * dXdmu_n;
        fjac(2, 3) += zion[n] * aion_inv[n] * dXdT;
        fjac(3, 1) += aion_inv[n] * dXdmu_p;
        fjac(3, 2) += aion_inv[n] * dXdmu_n;
        fjac(3, 3) += aion_inv[n] * dXdT;
    }

    // To save time, in the future, calculate X_k individually along with the jacobian,
    // then call eos after getting relevant X_k to find abar, dedA and dedZ

    fjac(3, 1) = - nse_state.abar * nse_state.abar *
        (nse_state.dedA + state_data.state.y_e * nse_state.dedZ) * fjac(3, 1);

    fjac(3, 2) = - nse_state.abar * nse_state.abar *
        (nse_state.dedA + state_data.state.y_e * nse_state.dedZ) * fjac(3, 2);

    fjac(3, 3) = nse_state.dedT - nse_state.abar * nse_state.abar *
        (nse_state.dedA + state_data.state.y_e * nse_state.dedZ) * fjac(3, 3);
}

#endif
