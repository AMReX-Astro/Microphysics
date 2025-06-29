#ifndef NSE_EQNS_H
#define NSE_EQNS_H

#include <AMReX_REAL.H>
#include <AMReX_Algorithm.H>

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
        // since we're using rhoX, make it sum to rho

        fvec(1) += nse_state.xn[n];
    }

    // constraint equation 2, electron fraction should be the same

    fvec(2) = nse_state.y_e - state_data.state.y_e;

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

#endif
