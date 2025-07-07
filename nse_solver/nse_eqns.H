#ifndef NSE_EQNS_H
#define NSE_EQNS_H

#include <AMReX_REAL.H>
#include <AMReX_Algorithm.H>
#include <nse_solver.H>
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
T get_T_independent_nonexponent_nse_state(const T& state) {
    // This function calculates a portion of the nse_state that is not dependent
    // on Temperature. This way we don't have to calculate the same
    // component repeatedly.

    T nse_state = state;

    // set spin

    amrex::Real spin = 1.0_rt;

    for (int n = 0; n < NumSpec; ++n) {
#ifdef NEW_NETWORK_IMPLEMENTATION
        if (n == NSE_INDEX::H1_index) {
            nse_state.xn[n] = 0.0_rt;
            continue;
        }
#endif

        // fill spin

#ifndef NEW_NETWORK_IMPLEMENTATION
        spin = get_spin_state(n+1);
#endif

        // find temperature independent parts

        nse_state.xn[n] = network::mion(n+1) * spin / state.rho;
    }

    return nse_state;
}


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void apply_T_dependent_nonexponent_nse_state(T& nse_state) {
    // This function applies the temperature dependent piece
    // that does not dependent on the chemical potential.

    // The input nse_state should be the state obtained from
    // get_T_independent_nse_state()

    amrex::Real T_in = nse_state.T_fixed > 0.0_rt ? nse_state.T_fixed : nse_state.T;
    constexpr amrex::Real kBi2PiHbar2 = C::k_B / (2.0_rt * M_PI * C::hbar * C::hbar);

    // Initialize temperature-dependent nuclear partition function

    amrex::Real pf = 1.0_rt;
    [[maybe_unused]] amrex::Real dpf_dT;
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

#ifndef NEW_NETWORK_IMPLEMENTATION
        get_partition_function(n+1, tfactors, pf, dpf_dT);
#endif

        // Apply temperature-dependent part

        amrex::Real power = network::mion(n+1) * T_in * kBi2PiHbar2;
        nse_state.xn[n] *= pf * std::sqrt(amrex::Math::powi<3>(power));
    }
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
    amrex::Real Gamma;

    for (int n = 0; n < NumSpec; ++n) {
#ifdef NEW_NETWORK_IMPLEMENTATION
        if (n == NSE_INDEX::H1_index) {
            continue;
        }
#endif
        // term for calculating u_c

        Gamma = std::cbrt(amrex::Math::powi<5>(zion[n])) * Gamma_e;

        // chemical potential for coulomb correction
        // see appendix of Calder 2007, doi:10.1086/510709 for more detail

        // reuse existing implementation from screening routine

        //
        // Here u_c is a dimensionless quantity.
        // i.e. normalized by kb T
        //

        u_c(n+1) = chabrier1998_helmholtz_F<amrex::Real>(Gamma);
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

    for (int n = 0; n < NumSpec; ++n) {
#ifdef NEW_NETWORK_IMPLEMENTATION
        if (n == NSE_INDEX::H1_index) {
            nse_state.xn[n] = 0.0_rt;
            continue;
        }
#endif

        // Apply the exponent part

        // prevent an overflow on exp by capping the exponent -- we hope that a subsequent
        // iteration will make it happy again

        amrex::Real exponent = amrex::min(500.0_rt,
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
    // Constraint Equations for (rho, T, Ye) input.

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

    // here state is the nse_state from get_T_independent_nonexponent_nse_state

    auto nse_state = state_data.state;
    nse_state.mu_p = x(1);
    nse_state.mu_n = x(2);
    nse_state.T    = x(3);

    // Finish evaluating NSE X_k

    // Apply the temperature dependent bit

    apply_T_dependent_nonexponent_nse_state(nse_state);

    // Find screening term

    amrex::Array1D<amrex::Real, 1, NumSpec> u_c{0.0_rt};

#if SCREEN_METHOD == SCREEN_METHOD_chabrier1998
    compute_coulomb_contribution(u_c, nse_state);
#endif

    // Apply exponent for calculating NSE mass fractions
    // This also updates y_e.

    apply_nse_exponent(nse_state, u_c);

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

    eos_t eos_state;
    for  (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = nse_state.xn[n];
    }
    eos_state.rho = nse_state.rho;
    eos_state.T = nse_state.T;
    eos(eos_input_rt, eos_state);

    // scale the energy constraint by e_in so that its on order 1
    // similar to other constraint eqs.

    fvec(3) = eos_state.e / state_data.state.e - 1.0_rt;
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
        // Pre-compute some relevant derivatives

        amrex::Real dX_dmu_p = nse_state.xn[n] * zion[n] * ikTMeV;
        amrex::Real dX_dmu_n = nse_state.xn[n] * (aion[n] - zion[n]) * ikTMeV;

        fjac(1, 1) += dX_dmu_p;
        fjac(1, 2) += dX_dmu_n;
        fjac(2, 1) += zion[n] * aion_inv[n] * dX_dmu_p;
        fjac(2, 2) += zion[n] * aion_inv[n] * dX_dmu_n;
    }

}


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_jcn(amrex::Array1D<amrex::Real, 1, 3>& x,
             amrex::Array2D<amrex::Real, 1, 3, 1, 3>& fjac,
             const T& state_data, int& iflag) {
    // 3x3 Constraint Jacobians for (rho, e, Ye) input.

    amrex::ignore_unused(iflag);

    // evaluate jacobian of the constraint eqns

    for (int i = 1; i <= 3; ++i) {
        for (int j = 1; j <= 3; ++j) {
            fjac(i, j) = 0.0_rt;
        }
    }

    // Initialize some variables

    auto nse_state = state_data.state;
    nse_state.mu_p = x(1);
    nse_state.mu_n = x(2);
    nse_state.T    = x(3);

    const amrex::Real T_in = nse_state.T;
    const amrex::Real e_in = state_data.state.e;
    const amrex::Real Ye = state_data.state.y_e;
    const amrex::Real ikTMeV = C::MeV2erg / (C::k_B * T_in);
    constexpr amrex::Real kBi2PiHbar2 = C::k_B / (2.0_rt * M_PI * C::hbar * C::hbar);

#if SCREEN_METHOD == SCREEN_METHOD_chabrier1998
    const amrex::Real n_e = state_data.state.rho * Ye * C::n_A;
    const amrex::Real Gamma_e = C::q_e * C::q_e *
        std::cbrt(1.333333333333_rt * M_PI * n_e) / (C::k_B * T_in);

    constexpr amrex::Real A_1 = -0.9052_rt;
    constexpr amrex::Real A_2 = 0.6322_rt;
    constexpr amrex::Real sqrt_A2 = gcem::sqrt(A_2);
    constexpr amrex::Real A_3 = -0.5_rt * gcem::sqrt(3.0_rt) - A_1 / sqrt_A2;
#endif

#ifndef NEW_NETWORK_IMPLEMENTATION
    auto tfactors = evaluate_tfactors(T_in);
#endif

    //
    // Evaluate the X_k and Jacobian at the same time to save runtime
    //

    for (int n = 0; n < NumSpec; ++n) {
#ifdef NEW_NETWORK_IMPLEMENTATION
        if (n == NSE_INDEX::H1_index) {
            continue;
        }
#endif
        // Find X_k for each individual nuclei

        // fill partition function and its temperature derivative

        amrex::Real pf = 1.0_rt;
        amrex::Real dpf_dT = 0.0_rt;

#ifndef NEW_NETWORK_IMPLEMENTATION
        get_partition_function(n+1, tfactors, pf, dpf_dT);
#endif

        // Apply temperature dependent part

        amrex::Real power = network::mion(n+1) * T_in * kBi2PiHbar2;
        nse_state.xn[n] *= pf * std::sqrt(amrex::Math::powi<3>(power));

        // Find coulomb potential and derivative

        amrex::Real u_c = 0.0_rt;
        amrex::Real du_c_dT = 0.0_rt;

#if SCREEN_METHOD == SCREEN_METHOD_chabrier1998

        // Find dΓ / dT = - Γ / T

        amrex::Real Gamma = std::cbrt(amrex::Math::powi<5>(zion[n])) * Gamma_e;
        amrex::Real dGamma_dT = - Gamma / T_in;

        // Compute coulomb potential

        u_c = chabrier1998_helmholtz_F<amrex::Real>(Gamma);

        // Find temperature derivative of the coulomb potential term: d (μ^c/kbT) / dT
        // Perhaps switch to using autodiff, but the expression seems simple enough.

        du_c_dT = dGamma_dT * std::sqrt(Gamma) * (A_3 / (1.0_rt + Gamma) +
                                                  A_1 / std::sqrt(A_2 + Gamma));
#endif

        // Apply exponent part to xn

        amrex::Real exponent = amrex::min(500.0_rt,
                                          (zion[n] * nse_state.mu_p + (aion[n] - zion[n]) *
                                           nse_state.mu_n + network::bion(n+1)) *
                                          ikTMeV - u_c);

        nse_state.xn[n] *= std::exp(exponent);

        // Now that NSE X's are computed, find Jacobian:

        // Pre-compute some relevant derivatives

        amrex::Real dX_dmu_p = nse_state.xn[n] * zion[n] * ikTMeV;
        amrex::Real dX_dmu_n = nse_state.xn[n] * (aion[n] - zion[n]) * ikTMeV;
        amrex::Real dX_dT = nse_state.xn[n] *
            (dpf_dT / pf - du_c_dT + (1.5_rt - (zion[n] * nse_state.mu_p +
                                                (aion[n] - zion[n]) * nse_state.mu_n +
                                                network::bion(n+1)) * ikTMeV) / T_in);

        fjac(1, 1) += dX_dmu_p;
        fjac(1, 2) += dX_dmu_n;
        fjac(1, 3) += dX_dT;
        fjac(2, 1) += zion[n] * aion_inv[n] * dX_dmu_p;
        fjac(2, 2) += zion[n] * aion_inv[n] * dX_dmu_n;
        fjac(2, 3) += zion[n] * aion_inv[n] * dX_dT;
        fjac(3, 1) += aion_inv[n] * dX_dmu_p;
        fjac(3, 2) += aion_inv[n] * dX_dmu_n;
        fjac(3, 3) += aion_inv[n] * dX_dT;
    }

    //
    // Now that NSE X's are all computed,
    // Create eos state to find relevant dedZ, dedA and abar via eos
    // This completes the Jacobian.
    //

    eos_extra_t eos_state;
    for  (int n = 0; n < NumSpec; ++n) {
        eos_state.xn[n] = nse_state.xn[n];
    }
    eos_state.rho = nse_state.rho;
    eos_state.T = T_in;
    eos(eos_input_rt, eos_state);

    fjac(3, 1) = - eos_state.abar * eos_state.abar *
        (eos_state.dedA + Ye * eos_state.dedZ) * fjac(3, 1) / e_in;

    fjac(3, 2) = - eos_state.abar * eos_state.abar *
        (eos_state.dedA + Ye * eos_state.dedZ) * fjac(3, 2) / e_in;

    fjac(3, 3) = (eos_state.dedT - eos_state.abar * eos_state.abar *
                  (eos_state.dedA + Ye * eos_state.dedZ) * fjac(3, 3)) / e_in;
}

#endif
