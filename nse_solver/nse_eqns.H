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
T get_T_independent_nse_state(const T& state) {
    // Computes (2J + 1) m^{5/2} (kB 1e9 / 2 π ħ^2)^{3/2} / ρ
    // a portion of the nse_state that is not dependent on Temperature.

    T nse_state = state;
    constexpr amrex::Real kBGK_2PiHbar2 = (C::k_B * 1.0e9_rt) / (2.0_rt * M_PI * C::hbar * C::hbar);

    amrex::constexpr_for<0, NumSpec>([&] (auto n)
    {
        // Note, using gcem::sqrt instead of std, would get underflow and get 0 here
        nse_state.xn[n] = std::sqrt(amrex::Math::powi<5>(network::mion<n+1>())) *
            get_spin_state<n+1>() * gcem::sqrt(amrex::Math::powi<3>(kBGK_2PiHbar2)) / state.rho;
    });

    return nse_state;
}


#if SCREEN_METHOD == SCREEN_METHOD_chabrier1998
template <int do_T9_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void compute_coulomb_contribution(amrex::Array2D<amrex::Real, 1, NumSpec, 1, 2, amrex::Order::C>& u_c,
                                  const T& state) {
    // This function computes the coulomb contribution and its T9 derivative

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

    for (int n = 0; n < NumSpec; ++n) {

        // term for calculating u_c

        amrex::Real Gamma = std::cbrt(amrex::Math::powi<5>(zion[n])) * Gamma_e;

        // chemical potential for coulomb correction
        // see appendix of Calder 2007, doi:10.1086/510709 for more detail

        // reuse existing implementation from screening routine

        //
        // Here u_c is a dimensionless quantity.
        // i.e. normalized by kb T
        //

        u_c(n+1, 1) = chabrier1998_helmholtz_F<amrex::Real>(Gamma);

        if constexpr (do_T9_derivatives) {
            constexpr amrex::Real A_1 = -0.9052_rt;
            constexpr amrex::Real A_2 = 0.6322_rt;
            constexpr amrex::Real A_3 = -0.5_rt * gcem::sqrt(3.0_rt) - A_1 / gcem::sqrt(A_2);

            // Find dΓ / dT9 = - Γ / T9

            amrex::Real dGamma_dT9 = - Gamma * 1.0e-9_rt / T_in;

            // Find temperature derivative of the coulomb potential term: d (μ^c/kbT) / dT
            // Perhaps switch to using autodiff, but the expression seems simple enough.

            u_c(n+1, 2) = dGamma_dT9 * std::sqrt(Gamma) *
                (A_3 / (1.0_rt + Gamma) + A_1 / std::sqrt(A_2 + Gamma));
        }
    }
}
#endif


template <int do_T9_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void evaluate_T_dependent_nse_state(amrex::Array2D<amrex::Real, 1, NumSpec, 1, 2, amrex::Order::C>& log_data,
                                    const T& state) {
    // Evaluates the temperature dependent piece and its T9 derivative.
    // i.e. (B - μ^c)/kBT + 1.5 lnT9 + ln(pf)
    // This is the term will be used inside the exp()

    amrex::Real T_in = state.T_fixed > 0.0_rt ? state.T_fixed : state.T;
    amrex::Real ikTMeV = 1.0_rt / (C::k_MeV * T_in);

    auto tfactors = evaluate_tfactors(T_in);

    // Evaluate the coulomb potential: μ^c / kBT and potentially find its T9 derivative

    amrex::Array2D<amrex::Real, 1, NumSpec, 1, 2, amrex::Order::C> u_c {0.0_rt};
#if SCREEN_METHOD == SCREEN_METHOD_chabrier1998
    compute_coulomb_contribution<do_T9_derivatives>(u_c, state);
#endif

    part_fun::pf_cache_t pf_cache;

    for (int n = 1; n <= NumSpec; ++n) {

        // Evaluate internal nuclear partition function, in log space

        amrex::Real logpf;
        [[maybe_unused]] amrex::Real dlogpf_dT9;
        get_partition_function_cached(n, tfactors, pf_cache, logpf, dlogpf_dT9);

        // Compute the log term

        log_data(n, 1) = network::bion(n) * ikTMeV - u_c(n, 1) +
            1.5_rt * tfactors.lnT9 + logpf;

        if constexpr (do_T9_derivatives) {
            log_data(n, 2) = (1.5_rt - network::bion(n) * ikTMeV) * tfactors.T9i +
                dlogpf_dT9 - u_c(n, 2);
        }
    }
}


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void apply_nse_exponent(T& nse_state,
                        const amrex::Array2D<amrex::Real, 1, NumSpec, 1, 2, amrex::Order::C>& log_data) {
    // This function applies the nse exponent which depends on
    // the chemical potential and the precomputed Temperature term
    // from evaluate_T_dependent_exponent_nse_state

    // if we are doing drive_initial_convection, we want to use
    // the temperature that comes in through T_fixed

    amrex::Real T_in = nse_state.T_fixed > 0.0_rt ? nse_state.T_fixed : nse_state.T;
    amrex::Real ikTMeV = 1.0_rt / (C::k_MeV * T_in);

    for (int n = 0; n < NumSpec; ++n) {

        // Apply the exponent part

        // prevent an overflow on exp by capping the exponent -- we hope that a subsequent
        // iteration will make it happy again

        amrex::Real exponent = amrex::min(500.0_rt,
                                          (zion[n] * nse_state.mu_p +
                                           (aion[n] - zion[n]) * nse_state.mu_n) *
                                          ikTMeV + log_data(n+1, 1));

        nse_state.xn[n] *= std::exp(exponent);
    }

    // Compute the corresponding electron fraction, Ye.

    nse_state.y_e = 0.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
        nse_state.y_e += nse_state.xn[n] * zion[n] * aion_inv[n];
    }
}

///
/// NSE constraint equations:
///
/// 1) Mass Conservation: Sum_k X_k - 1 = 0
/// 2) Charge Conservation: Sum_k Z_k X_k / A_k - Ye = 0
/// 3) Energy Conservation: e(rho, X_k^NSE, Ye) / e_in - 1 = 0
///
/// When the input is (rho, T, Ye), x = (μ_p, μ_n), solve for the first two Eqns
/// When the input is (rho, e, Ye), x = (μ_p, μ_n, T9), solve for all three Eqns.
/// the energy constraint is used to solve for Temperature.
///

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void nse_fcn(amrex::Array1D<amrex::Real, 1, 2>& x,
             amrex::Array1D<amrex::Real, 1, 2>& fvec,
             const T& state_data, int& iflag) {
    // Constraint Equations for (rho, T, Ye) input.
    // x = (μ_p, μ_n)

    amrex::ignore_unused(iflag);

    auto nse_state = state_data.state;
    nse_state.mu_p = x(1);
    nse_state.mu_n = x(2);

    // Apply exponent component for calculating nse mass fractions

    apply_nse_exponent(nse_state, state_data.log_data);

    // constraint equation 1, mass fraction sum to 1

    fvec(1) = -1.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
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
    // x = (μ_p, μ_n, T9)

    amrex::ignore_unused(iflag);

    // here state is the nse_state from get_T_independent_nse_state

    auto nse_state = state_data.state;
    nse_state.mu_p = x(1);
    nse_state.mu_n = x(2);
    nse_state.T    = x(3) * 1.0e9_rt;

    // Finish evaluating NSE X_k

    // Evaluate the temperature dependent piece, including screening.

    amrex::Array2D<amrex::Real, 1, NumSpec, 1, 2, amrex::Order::C> log_data{0.0_rt};
    constexpr int do_T9_derivatives = 0;
    evaluate_T_dependent_nse_state<do_T9_derivatives>(log_data, nse_state);

    // Apply exponent for calculating NSE mass fractions
    // This also updates y_e.

    apply_nse_exponent(nse_state, log_data);

    // constraint equation 1, mass fraction sum to 1

    fvec(1) = -1.0_rt;
    for (int n = 0; n < NumSpec; ++n) {
        fvec(1) += nse_state.xn[n];
    }

    // constraint equation 2, electron fraction should be the same
    // note that nse_state.y_e is already updated in apply_nse_exponent()

    fvec(2) = nse_state.y_e - state_data.state.y_e;

    // constraint equation 3, e(rho, T, X^NSE) / e_in - 1 = 0

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
/// When the input is (rho, e, Ye), x_j = {μ_p, μ_n, T9}
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

    apply_nse_exponent(nse_state, state_data.log_data);

    // evaluate jacobian of the constraint

    fjac(1, 1) = 0.0_rt;
    fjac(1, 2) = 0.0_rt;
    fjac(2, 1) = 0.0_rt;
    fjac(2, 2) = 0.0_rt;

    // if we are doing drive_initial_convection, we want to use
    // the temperature that comes in through T_fixed

    amrex::Real T_in = nse_state.T_fixed > 0.0_rt ? nse_state.T_fixed : nse_state.T;
    amrex::Real ikTMeV = 1.0_rt / (C::k_MeV * T_in);

    for (int n = 0; n < NumSpec; ++n) {

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
    nse_state.T    = x(3) * 1.0e9_rt;

    const amrex::Real T_in = nse_state.T;
    const amrex::Real T9i = 1.0_rt / x(3);
    const amrex::Real e_in = state_data.state.e;
    const amrex::Real Ye = state_data.state.y_e;
    const amrex::Real ikTMeV = 1.0_rt / (C::k_MeV * T_in);

    // Evaluate the temperature dependent piece and its T9 derivative

    amrex::Array2D<amrex::Real, 1, NumSpec, 1, 2, amrex::Order::C> log_data {0.0_rt};
    constexpr int do_T9_derivatives = 1;
    evaluate_T_dependent_nse_state<do_T9_derivatives>(log_data, nse_state);

    //
    // Evaluate the X_k and Jacobian at the same time to save runtime
    //

    for (int n = 0; n < NumSpec; ++n) {

        // Find X_k for each individual nuclei

        // Apply exponent part to xn

        amrex::Real mu_kBT = (zion[n] * nse_state.mu_p + (aion[n] - zion[n]) * nse_state.mu_n) * ikTMeV;
        amrex::Real exponent = amrex::min(500.0_rt, log_data(n+1,1) + mu_kBT);

        nse_state.xn[n] *= std::exp(exponent);

        // Now that NSE X's are computed, find Jacobian:

        // Pre-compute some relevant derivatives

        amrex::Real dX_dmu_p = nse_state.xn[n] * zion[n] * ikTMeV;
        amrex::Real dX_dmu_n = nse_state.xn[n] * (aion[n] - zion[n]) * ikTMeV;
        amrex::Real dX_dT9 = nse_state.xn[n] * (log_data(n+1, 2) - mu_kBT * T9i);

        fjac(1, 1) += dX_dmu_p;
        fjac(1, 2) += dX_dmu_n;
        fjac(1, 3) += dX_dT9;
        fjac(2, 1) += zion[n] * aion_inv[n] * dX_dmu_p;
        fjac(2, 2) += zion[n] * aion_inv[n] * dX_dmu_n;
        fjac(2, 3) += zion[n] * aion_inv[n] * dX_dT9;
        fjac(3, 1) += aion_inv[n] * dX_dmu_p;
        fjac(3, 2) += aion_inv[n] * dX_dmu_n;
        fjac(3, 3) += aion_inv[n] * dX_dT9;
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

    amrex::Real de_dT9 = 1.0e9_rt * eos_state.dedT;

    fjac(3, 1) = - eos_state.abar * eos_state.abar *
        (eos_state.dedA + Ye * eos_state.dedZ) * fjac(3, 1) / e_in;

    fjac(3, 2) = - eos_state.abar * eos_state.abar *
        (eos_state.dedA + Ye * eos_state.dedZ) * fjac(3, 2) / e_in;

    fjac(3, 3) = (de_dT9 - eos_state.abar * eos_state.abar *
                  (eos_state.dedA + Ye * eos_state.dedZ) * fjac(3, 3)) / e_in;
}

#endif
