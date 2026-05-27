#ifndef NSE_EQNS_H
#define NSE_EQNS_H

#include <AMReX_REAL.H>
#include <AMReX_Algorithm.H>
#include <extern_parameters.H>
#include <actual_network.H>
#include <partition_functions.H>

///
/// This file contains routines needed to solve:
/// NSE Composition Equations
///
/// X_k = m_k / ρ (2J_k + 1) G_k (m_k kB T / (2 π ħ^2))^{3/2}
///       × exp { Z_k μ_p + N_k μ_n + B_k - μ^c_k / (kB T)}
///

template <typename T>
struct nse_solver_data
{
    T state;

    // Holds the Temperature dependent log terms
    // and its derivative w.r.t T9, i.e. dlogX/dT9.
    amrex::Array1D<amrex::Real, 1, NumSpec> log_data;
    amrex::Array1D<amrex::Real, 1, NumSpec> dlogX_dT9;
};


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
T get_T_independent_nse_state(const T& state) {
    // Computes (2J + 1) m^{5/2} (kB 1e9 / 2 π ħ^2)^{3/2} / ρ
    // a portion of the nse_state that is not dependent on Temperature.

    T nse_state = state;
    constexpr amrex::Real kBGK_2PiHbar2 = (C::k_B * 1.0e9_rt) / (2.0_rt * M_PI * C::hbar * C::hbar);
    const amrex::Real rho_inv = 1.0_rt / state.rho;

    for (int n = 0; n < NumSpec; ++n) {
        // Absorb mion into kBGK_2PiHbar2 to prevent underflow at compile time.
        nse_state.xn[n] = network::mion(n+1) * get_spin_state(n+1) *
            std::sqrt(amrex::Math::powi<3>(network::mion(n+1) * kBGK_2PiHbar2)) * rho_inv;
    }

    return nse_state;
}


#if SCREEN_METHOD == SCREEN_METHOD_chabrier1998
template <int do_T9_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void compute_coulomb_contribution(amrex::Array1D<amrex::Real, 1, NumSpec>& u_c,
                                  amrex::Array1D<amrex::Real, 1, NumSpec>& du_c_dT9,
                                  const T& state) {
    // This function computes the coulomb contribution and its T9 derivative

    // if we use chabrier1998 screening
    // Get the required terms to calculate coulomb correction term, u_c

    amrex::Real T_in = state.T_fixed > 0.0_rt ? state.T_fixed : state.T;

    //
    // Find n_e for original state;
    // Note that y_e depends on the mass fraction,
    // but we use the coulomb correction to compute the mass fraction
    // But the one of the constraint is on input y_e
    // So here y_e is simply the actual y_e we want to achieve.
    // so u_c only depends on temperature.
    //

    const amrex::Real T9i = 1.0e9_rt / T_in;
    const amrex::Real n_e = state.rho * state.y_e * C::n_A;
    constexpr amrex::Real gamma_e_constants =
        C::q_e * C::q_e * gcem::pow(4.0_rt * M_PI / 3.0_rt, 1.0_rt / 3.0_rt) / (C::k_B * 1.0e9_rt);
    const amrex::Real Gamma_e = gamma_e_constants * std::cbrt(n_e) * T9i;

    for (int n = 1; n <= NumSpec; ++n) {

        // term for calculating u_c
        amrex::Real Z_53 = amrex::Math::powi<5>(std::cbrt(NetworkProperties::zion(n)));
        amrex::Real Gamma = Z_53 * Gamma_e;

        // chemical potential for coulomb correction
        // see appendix of Calder 2007, doi:10.1086/510709 for more detail

        // reuse existing implementation from screening routine

        //
        // Here u_c is a dimensionless quantity.
        // i.e. normalized by kB T
        //

        u_c(n) = chabrier1998_helmholtz_F<amrex::Real>(Gamma);

        if constexpr (do_T9_derivatives) {
            constexpr amrex::Real A_1 = -0.9052_rt;
            constexpr amrex::Real A_2 = 0.6322_rt;
            constexpr amrex::Real A_3 = -0.5_rt * gcem::sqrt(3.0_rt) - A_1 / gcem::sqrt(A_2);

            // Find dΓ / dT9 = - Γ / T9

            amrex::Real dGamma_dT9 = - Gamma * T9i;

            // Find temperature derivative of the coulomb potential term: d (μ^c/kBT) / dT
            // Perhaps switch to using autodiff, but the expression seems simple enough.

            du_c_dT9(n) = dGamma_dT9 * std::sqrt(Gamma) *
                (A_3 / (1.0_rt + Gamma) + A_1 / std::sqrt(A_2 + Gamma));
        }
    }
}
#endif


template <int do_T9_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_T_dependent_nse_state(nse_solver_data<T>& state_data) {
    // It does two things:
    // 1. Applies T9^{3/2} to the input state,
    //    where the input state should come from get_T_independent_nse_state()
    // 2. Evaluates the temperature dependent piece in the exponent
    //    i.e. - μ^c/kBT  + ln(pf)
    //    This is the term will be used inside the exp()
    //    It also evaluates the dlog/T9 derivative, including the effect of T9^{3/2}
    //    i.e. dlog/dT9 = 1.5 / T9 - du^c/dT9 + dlogpf/dT9

    auto& state = state_data.state;
    const amrex::Real T_in = state.T_fixed > 0.0_rt ? state.T_fixed : state.T;
    const amrex::Real T9 = T_in * 1.0e-9_rt;
    const amrex::Real T9_32 = std::sqrt(amrex::Math::powi<3>(T9));
    [[maybe_unused]] const amrex::Real T9i = 1.0_rt / T9;

    // Evaluate the coulomb potential: μ^c / kBT and potentially find its T9 derivative

    amrex::Array1D<amrex::Real, 1, NumSpec> u_c {0.0_rt};
    amrex::Array1D<amrex::Real, 1, NumSpec> du_c_dT9 {0.0_rt};
#if SCREEN_METHOD == SCREEN_METHOD_chabrier1998
    compute_coulomb_contribution<do_T9_derivatives>(u_c, du_c_dT9, state);
#endif

    // Create pf_cache for caching interpolation index.
    part_fun::pf_cache_t pf_cache;

    for (int n = 0; n < NumSpec; ++n) {

        // Evaluate internal nuclear partition function, in log space

        amrex::Real logpf;
        [[maybe_unused]] amrex::Real dlogpf_dT9;
        get_partition_function_cached(n+1, T9, pf_cache, logpf, dlogpf_dT9);

        // Compute the log term

        state_data.log_data(n+1) = logpf - u_c(n+1);

        // Apply T9^{3/2} to the state

        state.xn[n] *= T9_32;

        if constexpr (do_T9_derivatives) {
            state_data.dlogX_dT9(n+1) = 1.5_rt * T9i + dlogpf_dT9 - du_c_dT9(n+1);
        }
    }
}


template <int do_T9_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void apply_nse_exponent(nse_solver_data<T>& state_data) {
    // Last step to compute the NSE mass fractions.
    // It applies the nse exponent which depends on
    // the chemical potential and the precomputed Temperature term
    // from get_T_dependent_exponent_nse_state

    // This also completes dlogX_dT9 with the contribution from:
    // (Z_k μ_p + N_k μ_n + B_k) / (kB T)

    // if we are doing drive_initial_convection, we want to use
    // the temperature that comes in through T_fixed

    auto& nse_state = state_data.state;

    const amrex::Real T_in = nse_state.T_fixed > 0.0_rt ? nse_state.T_fixed : nse_state.T;
    const amrex::Real T9i = 1.0e9_rt / T_in;
    constexpr amrex::Real ikMeV_GK = 1.0e-9_rt / C::k_MeV;
    const amrex::Real ikTMeV = ikMeV_GK * T9i;

    for (int n = 0; n < NumSpec; ++n) {

        // Apply the exponent part

        // prevent an overflow on exp by capping the exponent -- we hope that a subsequent
        // iteration will make it happy again

        amrex::Real mu_kBT = (zion[n] * nse_state.mu_p +
                              (aion[n] - zion[n]) * nse_state.mu_n +
                              network::bion(n+1)) * ikTMeV;

        amrex::Real exponent = amrex::min(500.0_rt, state_data.log_data(n+1) + mu_kBT);
        nse_state.xn[n] *= std::exp(exponent);

        if constexpr (do_T9_derivatives) {
            state_data.dlogX_dT9(n+1) += -mu_kBT * T9i;
        }
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

    auto nse_state_data = state_data;
    auto& nse_state = nse_state_data.state;

    nse_state.mu_p = x(1);
    nse_state.mu_n = x(2);

    // Apply exponent component for calculating nse mass fractions
    constexpr int do_T9_derivatives = 0;
    apply_nse_exponent<do_T9_derivatives>(nse_state_data);

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

    // here nse_state from get_T_independent_nse_state
    auto nse_state_data = state_data;
    auto& nse_state = nse_state_data.state;

    nse_state.mu_p = x(1);
    nse_state.mu_n = x(2);
    nse_state.T    = x(3) * 1.0e9_rt;

    // Finish evaluating NSE X_k

    // get the temperature dependent piece, including screening.

    constexpr int do_T9_derivatives = 0;
    get_T_dependent_nse_state<do_T9_derivatives>(nse_state_data);

    // Apply exponent for calculating NSE mass fractions
    // This also updates y_e.

    apply_nse_exponent<do_T9_derivatives>(nse_state_data);

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

    auto nse_state_data = state_data;
    auto& nse_state = nse_state_data.state;

    nse_state.mu_p = x(1);
    nse_state.mu_n = x(2);

    constexpr int do_T9_derivatives = 0;
    apply_nse_exponent<do_T9_derivatives>(nse_state_data);

    // evaluate jacobian of the constraint

    fjac(1, 1) = 0.0_rt;
    fjac(1, 2) = 0.0_rt;
    fjac(2, 1) = 0.0_rt;
    fjac(2, 2) = 0.0_rt;

    // if we are doing drive_initial_convection, we want to use
    // the temperature that comes in through T_fixed

    const amrex::Real T_in = nse_state.T_fixed > 0.0_rt ? nse_state.T_fixed : nse_state.T;
    const amrex::Real ikTMeV = 1.0_rt / (C::k_MeV * T_in);

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

    auto nse_state_data = state_data;
    auto& nse_state = nse_state_data.state;

    nse_state.mu_p = x(1);
    nse_state.mu_n = x(2);
    nse_state.T    = x(3) * 1.0e9_rt;

    const amrex::Real T_in = nse_state.T;
    const amrex::Real T9i = 1.0_rt / x(3);
    const amrex::Real e_inv = 1.0_rt / state_data.state.e;
    const amrex::Real Ye = state_data.state.y_e;
    constexpr amrex::Real ikMeV_GK = 1.0e-9_rt / C::k_MeV;
    const amrex::Real ikTMeV = ikMeV_GK * T9i;

    // Evaluate the temperature dependent piece and its T9 derivative

    constexpr int do_T9_derivatives = 1;
    get_T_dependent_nse_state<do_T9_derivatives>(nse_state_data);

    // Apply exponent for calculating NSE mass fractions and
    // and completes dlog/dT9 derivative. This also updates y_e.

    apply_nse_exponent<do_T9_derivatives>(nse_state_data);

    //
    // Evaluate the Jacobian
    //

    for (int n = 0; n < NumSpec; ++n) {

        // Pre-compute some relevant derivatives

        amrex::Real dX_dmu_p = nse_state.xn[n] * zion[n] * ikTMeV;
        amrex::Real dX_dmu_n = nse_state.xn[n] * (aion[n] - zion[n]) * ikTMeV;
        amrex::Real dX_dT9 = nse_state.xn[n] * nse_state_data.dlogX_dT9(n+1);
        amrex::Real Z_over_A = zion[n] * aion_inv[n];

        fjac(1, 1) += dX_dmu_p;
        fjac(1, 2) += dX_dmu_n;
        fjac(1, 3) += dX_dT9;
        fjac(2, 1) += Z_over_A * dX_dmu_p;
        fjac(2, 2) += Z_over_A * dX_dmu_n;
        fjac(2, 3) += Z_over_A * dX_dT9;
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

    // de/dT9 normalize by input energy.
    amrex::Real de_dT9_norm = 1.0e9_rt * eos_state.dedT * e_inv;

    // -abar^2 * (de/dabar + Ye de/dzbar) normalized by input energy
    amrex::Real common_term = -eos_state.abar * eos_state.abar *
        (eos_state.dedA + Ye * eos_state.dedZ) * e_inv;

    fjac(3, 1) *= common_term;
    fjac(3, 2) *= common_term;
    fjac(3, 3) = de_dT9_norm + common_term * fjac(3, 3);
}

#endif
