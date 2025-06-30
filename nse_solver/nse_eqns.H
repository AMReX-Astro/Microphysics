#ifndef NSE_EQNS_H
#define NSE_EQNS_H

#include <AMReX_REAL.H>
#include <AMReX_Algorithm.H>

#include <extern_parameters.H>

///
/// This file contains routines needed to solve:
/// NSE Composition Equations
///
/// X_k = m_k / ρ (2J_k + 1) G_k (m_k kB T / (2 π ħ^2))^{3/2}
///       × exp { Z_k μ_p + N_k μ_n + B_k - μ^c_k / (kB T)}
///

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
T get_nonexponent_nse_state(const T& state) {
    // This function calculates a portion of the nse_state that is not dependent
    // on the chemical potential. This way we don't have to calculate the same
    // component repeatedly.

    // This function returns an incomplete nse_state where the xn lacks the exponent term.

    T nse_state = state;

    // set partition function and spin

    amrex::Real pf = 1.0_rt;
    [[maybe_unused]] amrex::Real dpf_dT;
    amrex::Real spin = 1.0_rt;

    // if we are doing drive_initial_convection, we want to use
    // the temperature that comes in through T_fixed

    amrex::Real T_in = state.T_fixed > 0.0_rt ? state.T_fixed : state.T;
    constexpr amrex::Real i2PiHbar2 = 1.0_rt / (2.0_rt * M_PI * C::hbar * C::hbar);

#ifndef NEW_NETWORK_IMPLEMENTATION
    auto tfactors = evaluate_tfactors(T_in);
#endif

    for (int n = 0; n < NumSpec; ++n) {
#ifdef NEW_NETWORK_IMPLEMENTATION
        if (n == NSE_INDEX::H1_index) {
            nse_state.xn[n] = 0.0_rt;
            continue;
        }
#endif

        // fill partition function and get spin

#ifndef NEW_NETWORK_IMPLEMENTATION
        spin = get_spin_state(n+1);
        get_partition_function(n+1, tfactors, pf, dpf_dT);
#endif

        // find nse mass frac without the exponent term.

        amrex::Real power = network::mion(n+1) * C::k_B * T_in * i2PiHbar2;

        nse_state.xn[n] = network::mion(n+1) * pf * spin / state.rho *
            std::sqrt(amrex::Math::powi<3>(power));
    }

    return nse_state;
}


#if SCREEN_METHOD == SCREEN_METHOD_chabrier1998
template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void compute_coulomb_contribution(amrex::Array1D<amrex::Real, 1, NumSpec>& u_c,
                                const T& state) {
    // This function computes the coulomb contribution 1D array

    // if we use chabrier1998 screening
    // Get the required terms to calculate coulomb correction term, u_c

    amrex::Real T_in = state.T_fixed > 0.0_rt ? state.T_fixed : state.T;

    //
    // Find n_e for original state;
    // Note that y_e depends on the mass fraction,
    // but we use the coulomb correction to compute the mass fraction
    // So here y_e is simply the actual y_e we want to achieve.
    // so we just treat u_c as a constant.
    //

    const amrex::Real n_e = state.rho * state.y_e * C::n_A;
    const amrex::Real Gamma_e = C::q_e * C::q_e *
        std::cbrt(1.333333333333_rt * M_PI * n_e) / (C::k_B * T_in);
    amrex::Real gamma;

    for (int n = 0; n < NumSpec; ++n) {
#ifdef NEW_NETWORK_IMPLEMENTATION
        if (n == NSE_INDEX::H1_index) {
            continue;
        }
#endif
        // term for calculating u_c

        gamma = std::cbrt(amrex::Math::powi<5>(zion[n])) * Gamma_e;

        // chemical potential for coulomb correction
        // see appendix of Calder 2007, doi:10.1086/510709 for more detail

        // reuse existing implementation from screening routine
        auto f = chabrier1998_helmholtz_F<amrex::Real>(gamma);

        //
        // Here u_c is a dimensionless quantity.
        // Otherwise:
        // u_c = C::k_B * T_in / C::MeV2erg * f;
        //

        u_c(n+1) = f;
    }
}
#endif


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void apply_nse_exponent(T& nse_state,
                        const amrex::Array1D<amrex::Real, 1, NumSpec>& u_c) {
    // This function applies the nse exponent which depends on
    // the chemical potential.

    // if we are doing drive_initial_convection, we want to use
    // the temperature that comes in through T_fixed

    amrex::Real T_in = nse_state.T_fixed > 0.0_rt ? nse_state.T_fixed : nse_state.T;
    amrex::Real ikTMeV = C::MeV2erg / (C::k_B * T_in);
    amrex::Real exponent;

    for (int n = 0; n < NumSpec; ++n) {
#ifdef NEW_NETWORK_IMPLEMENTATION
        if (n == NSE_INDEX::H1_index) {
            nse_state.xn[n] = 0.0_rt;
            continue;
        }
#endif

        // find nse mass frac

        // prevent an overflow on exp by capping the exponent -- we hope that a subsequent
        // iteration will make it happy again

        exponent = amrex::min(500.0_rt,
                              (zion[n] * nse_state.mu_p + (aion[n] - zion[n]) *
                               nse_state.mu_n + network::bion(n+1)) *
                              ikTMeV - u_c(n+1));

        nse_state.xn[n] *= std::exp(exponent);
    }

    nse_state.y_e = 0.0_rt;

    for (int n = 0; n < NumSpec; ++n) {
#ifdef NEW_NETWORK_IMPLEMENTATION
        if (n == NSE_INDEX::H1_index) {
            continue;
        }
#endif
        // constraint equation 1, mass fraction sum to 1

        nse_state.y_e += nse_state.xn[n] * zion[n] * aion_inv[n];
    }
}

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
