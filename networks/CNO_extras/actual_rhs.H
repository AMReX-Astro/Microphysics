#ifndef ACTUAL_RHS_H
#define ACTUAL_RHS_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <extern_parameters.H>
#include <actual_network.H>
#include <burn_type.H>
#include <jacobian_utilities.H>
#ifdef SCREENING
#include <screen.H>
#endif
#include <microphysics_autodiff.H>
#ifdef NEUTRINOS
#include <neutrino.H>
#endif
#include <rate_type.H>
#include <reaclib_rates.H>
#include <approximate_rates.H>
#include <modified_rates.H>
#include <branched_rates.H>
#include <temperature_table_rates.H>
#include <tabular_weak_rates.H>
#include <derived_rates.H>

using namespace amrex::literals;

template<class T>
AMREX_GPU_DEVICE AMREX_INLINE
void ener_gener_rate(T const& dydt, amrex::Real& enuc)
{

    // Computes the instantaneous energy generation rate (from the nuclei)

    // This is basically e = m c**2

    enuc = 0.0_rt;

    amrex::constexpr_for<1, NumSpec+1>([&] (auto n)
    {
        enuc += dydt(n) * network::mion<n>();
    });

    enuc *= C::enuc_conv2;
}


#ifdef SCREENING
template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void evaluate_screening(const burn_t& state,
                        const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
                        T& rate_eval) {

    using namespace Rates;

    // Evaluate screening terms

    using number_t = std::conditional_t<do_T_derivatives, autodiff::dual, amrex::Real>;
    number_t temp = state.T;
    if constexpr (do_T_derivatives) {
        // seed the dual number for temperature before calculating anything with it
        autodiff::seed(temp);
    }

    plasma_state_t<number_t> pstate{};
    fill_plasma_state(pstate, temp, state.rho, Y);

    amrex::Real log_scor, dlog_scor_dT;

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 10.0_rt, 18.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Ne18) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Ne18) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 9.0_rt, 18.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_F18) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_F18) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 6.0_rt, 12.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_C12) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_C12) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 11.0_rt, 21.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Na21) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Na21) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 8.0_rt, 16.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_O16) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_O16) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 3.0_rt, 2.0_rt, 3.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He3_He3) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He3_He3) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 9.0_rt, 17.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_F17) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_F17) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 7.0_rt, 14.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_N14) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_N14) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 3.0_rt, 2.0_rt, 4.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He3_He4) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He3_He4) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 9.0_rt, 19.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_F19) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_F19) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 6.0_rt, 13.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_C13) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_C13) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 6.0_rt, 12.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_C12) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_C12) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 1.0_rt, 1.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_p) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_p) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 8.0_rt, 17.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_O17) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_O17) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 7.0_rt, 15.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_N15) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_N15) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 10.0_rt, 20.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Ne20) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Ne20) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 8.0_rt, 18.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_O18) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_O18) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 2.0_rt, 3.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_He3) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_He3) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 7.0_rt, 15.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_N15) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_N15) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 10.0_rt, 20.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Ne20) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Ne20) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 7.0_rt, 13.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_N13) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_N13) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 8.0_rt, 14.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_O14) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_O14) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 2.0_rt, 4.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_He4) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_He4) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 4.0_rt, 8.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Be8) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_Be8) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 8.0_rt, 16.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_O16) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_O16) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 7.0_rt, 14.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_N14) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_N14) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 9.0_rt, 17.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_F17) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_F17) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 10.0_rt, 19.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Ne19) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_p_Ne19) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 8.0_rt, 15.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_O15) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_O15) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 7.0_rt, 13.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_N13) = log_scor;
        if constexpr (do_T_derivatives) {
            rate_eval.dlog_screen_dT(k_He4_N13) = dlog_scor_dT;
        }
    }

}
#endif


template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void evaluate_rates(const burn_t& state,
                    const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
                    T& rate_eval) {

    using namespace Rates;

    // here Y is consistent with state.xn

    [[maybe_unused]] const amrex::Real rhoy = state.rho * state.y_e;

    const tf_t tfactors = evaluate_tfactors(state.T);

    // Precompute screening terms

#ifdef SCREENING
    evaluate_screening<do_T_derivatives, T>(state, Y, rate_eval);
#endif

    // Fill in different rates

    reaclib_rates::fill_rates<do_T_derivatives, T>(tfactors, rate_eval);

    tabular_weak_rates::fill_rates(state.T, rhoy, Y, rate_eval);

    temp_tabular::fill_rates<do_T_derivatives, T>(tfactors, rate_eval);

    // fill modified rates next -- these can have ReacLib or
    // TemperatureTabular/StarLib rates as the original rate
    modified_rates::fill_rates<do_T_derivatives, T>(tfactors, rate_eval);

    // fill branched rates -- these can have ReacLib or
    // TemperatureTabular/StarLib as the underlying rate or branches.
    branched_rates::fill_rates<do_T_derivatives, T>(tfactors, rate_eval);

    // Calculate Derived Rates next. This should go last but before
    // approx rates, since those may have ReacLibRate,
    // TemperatureTabularRate / StarLibRate, or DerivedRate as the
    // underlying rate.
    derived_rates::fill_rates<do_T_derivatives, T>(tfactors, rate_eval);



    // Fill approximate rates
    approximate_rates::fill_rates<do_T_derivatives, T>(tfactors, state.rho, Y, rate_eval);

}


// backwards-compatible wrapper
template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void evaluate_rates(const burn_t& state,
                    T& rate_eval) {

    // create molar fractions
    amrex::Array1D<amrex::Real, 1, NumSpec> Y;
    for (int n = 1; n <= NumSpec; ++n) {
        Y(n) = state.xn[n-1] * aion_inv[n-1];
    }

    evaluate_rates<do_T_derivatives>(state, Y, rate_eval);
}


#ifdef NSE_NET
AMREX_GPU_DEVICE AMREX_INLINE
void get_ydot_weak(const burn_t& state,
                   amrex::Array1D<amrex::Real, 1, neqs>& ydot_nuc,
                   amrex::Real& enuc_weak,
                   [[maybe_unused]] const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {
    ///
    /// Calculate Ydots contribute only from weak reactions.
    /// This is used to calculate dyedt and energy generation from
    /// weak reactions for self-consistent NSE
    ///

    using namespace modified_rates;
    using namespace temp_tabular;

    // initialize ydot_nuc to 0

    for (int i = 1; i <= neqs; ++i) {
        ydot_nuc(i) = 0.0_rt;
    }

    rate_t rate_eval;

    constexpr int do_T_derivatives = 0;
    [[maybe_unused]] amrex::Real rate, drate_dT, edot_nu, edot_gamma;
    [[maybe_unused]] amrex::Real rhoy = state.rho * state.y_e;

    rate_eval.enuc_weak = 0.0_rt;

    // Compute all weak rates and get ydot_weak

#ifdef SCREENING
    {
        plasma_state_t<amrex::Real> pstate{};
        fill_plasma_state(pstate, state.T, state.rho, Y);
        amrex::Real log_scor, dlog_scor_dT;
        {
            constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 3.0_rt, 2.0_rt, 4.0_rt);
            static_assert(scn_fac.z1 == 2.0_rt);
            actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
            rate_eval.log_screen(k_He3_He4) = log_scor;
        }

        {
            constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 10.0_rt, 19.0_rt);
            static_assert(scn_fac.z1 == 1.0_rt);
            actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
            rate_eval.log_screen(k_p_Ne19) = log_scor;
        }

        {
            constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 2.0_rt, 3.0_rt);
            static_assert(scn_fac.z1 == 1.0_rt);
            actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
            rate_eval.log_screen(k_p_He3) = log_scor;
        }

        {
            constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 1.0_rt, 1.0_rt);
            static_assert(scn_fac.z1 == 1.0_rt);
            actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
            rate_eval.log_screen(k_p_p) = log_scor;
        }

    }
#endif

    const tf_t tfactors = evaluate_tfactors(state.T);

    {
        // N13_to_C13_beta_pos_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        reaclib_rates::rate_N13_to_C13_beta_pos_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_N13_to_C13_beta_pos_reaclib) = rate;
    }

    {
        // O14_to_N14_beta_pos_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        reaclib_rates::rate_O14_to_N14_beta_pos_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_O14_to_N14_beta_pos_reaclib) = rate;
    }

    {
        // O15_to_N15_beta_pos_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        reaclib_rates::rate_O15_to_N15_beta_pos_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_O15_to_N15_beta_pos_reaclib) = rate;
    }

    {
        // p_p_to_d_beta_pos_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_p);
#endif
        reaclib_rates::rate_p_p_to_d_beta_pos_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_p_to_d_beta_pos_reaclib) = rate;
    }

    {
        // p_p_to_d_electron_capture_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_p);
#endif
        reaclib_rates::rate_p_p_to_d_electron_capture_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_p_to_d_electron_capture_reaclib) = rate;
    }

    {
        // He3_p_to_He4_beta_pos_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_He3);
#endif
        reaclib_rates::rate_He3_p_to_He4_beta_pos_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He3_p_to_He4_beta_pos_reaclib) = rate;
    }

    {
        // Ne19_3p_to_Mg22_beta_pos_modified

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Ne19);
#endif
        modified_rates::rate_Ne19_3p_to_Mg22_beta_pos_modified<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_Ne19_3p_to_Mg22_beta_pos_modified) = rate;
    }

    {
        // p_p_p_to_He3_beta_pos_modified

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_p);
#endif
        modified_rates::rate_p_p_p_to_He3_beta_pos_modified<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_p_p_to_He3_beta_pos_modified) = rate;
    }

    {
        // p_p_p_to_He3_electron_capture_modified

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_p);
#endif
        modified_rates::rate_p_p_p_to_He3_electron_capture_modified<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_p_p_to_He3_electron_capture_modified) = rate;
    }

    {
        // He3_He4_p_to_He4_He4_beta_pos_modified

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He3_He4);
#endif
        modified_rates::rate_He3_He4_p_to_He4_He4_beta_pos_modified<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He3_He4_p_to_He4_He4_beta_pos_modified) = rate;
    }

    tabular_weak_rates::fill_rates(state.T, rhoy, Y, rate_eval);

    const auto& screened_rates = rate_eval.screened_rates;

    ydot_nuc(H1) =
        -3.0*screened_rates(k_Ne19_3p_to_Mg22_beta_pos_modified)*Y(Ne19)*Y(H1)*state.rho +
        -1.5*screened_rates(k_p_p_p_to_He3_beta_pos_modified)*amrex::Math::powi<2>(Y(H1))*state.rho +
        -1.5*screened_rates(k_p_p_p_to_He3_electron_capture_modified)*amrex::Math::powi<2>(Y(H1))*state.rho +
        -screened_rates(k_He3_p_to_He4_beta_pos_reaclib)*Y(He3)*Y(H1)*state.rho +
        -screened_rates(k_He3_He4_p_to_He4_He4_beta_pos_modified)*Y(He3)*Y(He4)*state.rho;

    ydot_nuc(He3) =
        0.5*screened_rates(k_p_p_p_to_He3_beta_pos_modified)*amrex::Math::powi<2>(Y(H1))*state.rho +
        0.5*screened_rates(k_p_p_p_to_He3_electron_capture_modified)*amrex::Math::powi<2>(Y(H1))*state.rho +
        -screened_rates(k_He3_p_to_He4_beta_pos_reaclib)*Y(He3)*Y(H1)*state.rho +
        -screened_rates(k_He3_He4_p_to_He4_He4_beta_pos_modified)*Y(He3)*Y(He4)*state.rho;

    ydot_nuc(He4) =
        screened_rates(k_He3_p_to_He4_beta_pos_reaclib)*Y(He3)*Y(H1)*state.rho +
        screened_rates(k_He3_He4_p_to_He4_He4_beta_pos_modified)*Y(He3)*Y(He4)*state.rho;

    ydot_nuc(C12) = 0.0_rt;

    ydot_nuc(C13) =
        screened_rates(k_N13_to_C13_beta_pos_reaclib)*Y(N13);

    ydot_nuc(N13) =
        -screened_rates(k_N13_to_C13_beta_pos_reaclib)*Y(N13);

    ydot_nuc(N14) =
        screened_rates(k_O14_to_N14_beta_pos_reaclib)*Y(O14);

    ydot_nuc(N15) =
        screened_rates(k_O15_to_N15_beta_pos_reaclib)*Y(O15);

    ydot_nuc(O14) =
        -screened_rates(k_O14_to_N14_beta_pos_reaclib)*Y(O14);

    ydot_nuc(O15) =
        -screened_rates(k_O15_to_N15_beta_pos_reaclib)*Y(O15);

    ydot_nuc(O16) = 0.0_rt;

    ydot_nuc(O17) =
        (screened_rates(k_F17_to_O17_electron_capture_weaktab)*Y(F17) + -screened_rates(k_O17_to_F17_beta_neg_weaktab)*Y(O17));

    ydot_nuc(O18) =
        (screened_rates(k_F18_to_O18_electron_capture_weaktab)*Y(F18) + -screened_rates(k_O18_to_F18_beta_neg_weaktab)*Y(O18));

    ydot_nuc(F17) =
        (-screened_rates(k_F17_to_O17_electron_capture_weaktab)*Y(F17) + screened_rates(k_O17_to_F17_beta_neg_weaktab)*Y(O17));

    ydot_nuc(F18) =
        (-screened_rates(k_F18_to_O18_electron_capture_weaktab)*Y(F18) + screened_rates(k_O18_to_F18_beta_neg_weaktab)*Y(O18)) +
        (screened_rates(k_Ne18_to_F18_electron_capture_weaktab)*Y(Ne18) + -screened_rates(k_F18_to_Ne18_beta_neg_weaktab)*Y(F18));

    ydot_nuc(F19) =
        (screened_rates(k_Ne19_to_F19_electron_capture_weaktab)*Y(Ne19) + -screened_rates(k_F19_to_Ne19_beta_neg_weaktab)*Y(F19));

    ydot_nuc(Ne18) =
        (-screened_rates(k_Ne18_to_F18_electron_capture_weaktab)*Y(Ne18) + screened_rates(k_F18_to_Ne18_beta_neg_weaktab)*Y(F18));

    ydot_nuc(Ne19) =
        -screened_rates(k_Ne19_3p_to_Mg22_beta_pos_modified)*Y(Ne19)*Y(H1)*state.rho +
        (-screened_rates(k_Ne19_to_F19_electron_capture_weaktab)*Y(Ne19) + screened_rates(k_F19_to_Ne19_beta_neg_weaktab)*Y(F19));

    ydot_nuc(Ne20) = 0.0_rt;

    ydot_nuc(Mg22) =
        screened_rates(k_Ne19_3p_to_Mg22_beta_pos_modified)*Y(Ne19)*Y(H1)*state.rho;

    ydot_nuc(Mg24) = 0.0_rt;

    ydot_nuc(Fe56) = 0.0_rt;

    enuc_weak = rate_eval.enuc_weak;
}
#endif


AMREX_GPU_DEVICE AMREX_INLINE
void rhs_nuc(const burn_t& state,
             amrex::Array1D<amrex::Real, 1, neqs>& ydot_nuc,
             const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
             const amrex::Array1D<amrex::Real, 1, Rates::NumRates>& screened_rates) {

    using namespace Rates;

    ydot_nuc(H1) =
        (-screened_rates(k_C12_p_to_N13_reaclib)*Y(C12)*Y(H1)*state.rho + screened_rates(k_N13_to_p_C12_derived)*Y(N13)) +
        (-screened_rates(k_C13_p_to_N14_reaclib)*Y(C13)*Y(H1)*state.rho + screened_rates(k_N14_to_p_C13_derived)*Y(N14)) +
        (-screened_rates(k_N13_p_to_O14_reaclib)*Y(N13)*Y(H1)*state.rho + screened_rates(k_O14_to_p_N13_derived)*Y(O14)) +
        (-screened_rates(k_N14_p_to_O15_reaclib)*Y(N14)*Y(H1)*state.rho + screened_rates(k_O15_to_p_N14_derived)*Y(O15)) +
        (-screened_rates(k_N15_p_to_O16_reaclib)*Y(N15)*Y(H1)*state.rho + screened_rates(k_O16_to_p_N15_derived)*Y(O16)) +
        (-screened_rates(k_O16_p_to_F17_reaclib)*Y(O16)*Y(H1)*state.rho + screened_rates(k_F17_to_p_O16_derived)*Y(F17)) +
        (-screened_rates(k_O17_p_to_F18_reaclib)*Y(O17)*Y(H1)*state.rho + screened_rates(k_F18_to_p_O17_derived)*Y(F18)) +
        (-screened_rates(k_O18_p_to_F19_reaclib)*Y(O18)*Y(H1)*state.rho + screened_rates(k_F19_to_p_O18_derived)*Y(F19)) +
        (-screened_rates(k_F17_p_to_Ne18_reaclib)*Y(F17)*Y(H1)*state.rho + screened_rates(k_Ne18_to_p_F17_derived)*Y(Ne18)) +
        (-screened_rates(k_F18_p_to_Ne19_reaclib)*Y(F18)*Y(H1)*state.rho + screened_rates(k_Ne19_to_p_F18_derived)*Y(Ne19)) +
        (-screened_rates(k_F19_p_to_Ne20_reaclib)*Y(F19)*Y(H1)*state.rho + screened_rates(k_Ne20_to_p_F19_derived)*Y(Ne20)) +
        (screened_rates(k_N13_He4_to_p_O16_reaclib)*Y(He4)*Y(N13)*state.rho + -screened_rates(k_O16_p_to_He4_N13_derived)*Y(O16)*Y(H1)*state.rho) +
        (-screened_rates(k_N15_p_to_He4_C12_reaclib)*Y(N15)*Y(H1)*state.rho + screened_rates(k_C12_He4_to_p_N15_derived)*Y(C12)*Y(He4)*state.rho) +
        (screened_rates(k_O14_He4_to_p_F17_reaclib)*Y(He4)*Y(O14)*state.rho + -screened_rates(k_F17_p_to_He4_O14_derived)*Y(F17)*Y(H1)*state.rho) +
        (-screened_rates(k_O17_p_to_He4_N14_reaclib)*Y(O17)*Y(H1)*state.rho + screened_rates(k_N14_He4_to_p_O17_derived)*Y(He4)*Y(N14)*state.rho) +
        (-screened_rates(k_O18_p_to_He4_N15_reaclib)*Y(O18)*Y(H1)*state.rho + screened_rates(k_N15_He4_to_p_O18_derived)*Y(He4)*Y(N15)*state.rho) +
        (-screened_rates(k_F18_p_to_He4_O15_reaclib)*Y(F18)*Y(H1)*state.rho + screened_rates(k_O15_He4_to_p_F18_derived)*Y(He4)*Y(O15)*state.rho) +
        (-screened_rates(k_F19_p_to_He4_O16_reaclib)*Y(F19)*Y(H1)*state.rho + screened_rates(k_O16_He4_to_p_F19_derived)*Y(He4)*Y(O16)*state.rho) +
        (screened_rates(k_F17_He4_to_p_Ne20_derived)*Y(F17)*Y(He4)*state.rho + -screened_rates(k_Ne20_p_to_He4_F17_reaclib)*Y(Ne20)*Y(H1)*state.rho) +
        -3.0*screened_rates(k_Ne19_3p_to_Mg22_beta_pos_modified)*Y(Ne19)*Y(H1)*state.rho +
        -2.0*screened_rates(k_Ne20_p_p_to_Mg22_modified)*Y(Ne20)*Y(H1)*state.rho +
        -1.5*screened_rates(k_p_p_p_to_He3_beta_pos_modified)*amrex::Math::powi<2>(Y(H1))*state.rho +
        -1.5*screened_rates(k_p_p_p_to_He3_electron_capture_modified)*amrex::Math::powi<2>(Y(H1))*state.rho +
        screened_rates(k_He3_He3_to_p_p_He4_reaclib)*amrex::Math::powi<2>(Y(He3))*state.rho +
        -screened_rates(k_He3_p_to_He4_beta_pos_reaclib)*Y(He3)*Y(H1)*state.rho +
        -screened_rates(k_He3_He4_p_to_He4_He4_beta_pos_modified)*Y(He3)*Y(He4)*state.rho;

    ydot_nuc(He3) =
        0.5*screened_rates(k_p_p_p_to_He3_beta_pos_modified)*amrex::Math::powi<2>(Y(H1))*state.rho +
        0.5*screened_rates(k_p_p_p_to_He3_electron_capture_modified)*amrex::Math::powi<2>(Y(H1))*state.rho +
        -screened_rates(k_He3_He3_to_p_p_He4_reaclib)*amrex::Math::powi<2>(Y(He3))*state.rho +
        -screened_rates(k_He3_p_to_He4_beta_pos_reaclib)*Y(He3)*Y(H1)*state.rho +
        -screened_rates(k_He3_He4_p_to_He4_He4_beta_pos_modified)*Y(He3)*Y(He4)*state.rho;

    ydot_nuc(He4) =
        (-screened_rates(k_C12_He4_to_O16_reaclib)*Y(C12)*Y(He4)*state.rho + screened_rates(k_O16_to_He4_C12_derived)*Y(O16)) +
        (-screened_rates(k_N14_He4_to_F18_reaclib)*Y(He4)*Y(N14)*state.rho + screened_rates(k_F18_to_He4_N14_derived)*Y(F18)) +
        (-screened_rates(k_N15_He4_to_F19_reaclib)*Y(He4)*Y(N15)*state.rho + screened_rates(k_F19_to_He4_N15_derived)*Y(F19)) +
        (-screened_rates(k_O14_He4_to_Ne18_reaclib)*Y(He4)*Y(O14)*state.rho + screened_rates(k_Ne18_to_He4_O14_derived)*Y(Ne18)) +
        (-screened_rates(k_O15_He4_to_Ne19_reaclib)*Y(He4)*Y(O15)*state.rho + screened_rates(k_Ne19_to_He4_O15_derived)*Y(Ne19)) +
        (-screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(He4)*Y(O16)*state.rho + screened_rates(k_Ne20_to_He4_O16_derived)*Y(Ne20)) +
        (-screened_rates(k_Ne20_He4_to_Mg24_reaclib)*Y(He4)*Y(Ne20)*state.rho + screened_rates(k_Mg24_to_He4_Ne20_derived)*Y(Mg24)) +
        (-screened_rates(k_N13_He4_to_p_O16_reaclib)*Y(He4)*Y(N13)*state.rho + screened_rates(k_O16_p_to_He4_N13_derived)*Y(O16)*Y(H1)*state.rho) +
        (screened_rates(k_N15_p_to_He4_C12_reaclib)*Y(N15)*Y(H1)*state.rho + -screened_rates(k_C12_He4_to_p_N15_derived)*Y(C12)*Y(He4)*state.rho) +
        (-screened_rates(k_O14_He4_to_p_F17_reaclib)*Y(He4)*Y(O14)*state.rho + screened_rates(k_F17_p_to_He4_O14_derived)*Y(F17)*Y(H1)*state.rho) +
        (screened_rates(k_O17_p_to_He4_N14_reaclib)*Y(O17)*Y(H1)*state.rho + -screened_rates(k_N14_He4_to_p_O17_derived)*Y(He4)*Y(N14)*state.rho) +
        (screened_rates(k_O18_p_to_He4_N15_reaclib)*Y(O18)*Y(H1)*state.rho + -screened_rates(k_N15_He4_to_p_O18_derived)*Y(He4)*Y(N15)*state.rho) +
        (screened_rates(k_F18_p_to_He4_O15_reaclib)*Y(F18)*Y(H1)*state.rho + -screened_rates(k_O15_He4_to_p_F18_derived)*Y(He4)*Y(O15)*state.rho) +
        (screened_rates(k_F19_p_to_He4_O16_reaclib)*Y(F19)*Y(H1)*state.rho + -screened_rates(k_O16_He4_to_p_F19_derived)*Y(He4)*Y(O16)*state.rho) +
        (-0.5*screened_rates(k_He4_He4_He4_to_C12_reaclib)*amrex::Math::powi<3>(Y(He4))*amrex::Math::powi<2>(state.rho) + 3.0*screened_rates(k_C12_to_3He4_derived)*Y(C12)) +
        (-screened_rates(k_F17_He4_to_p_Ne20_derived)*Y(F17)*Y(He4)*state.rho + screened_rates(k_Ne20_p_to_He4_F17_reaclib)*Y(Ne20)*Y(H1)*state.rho) +
        (-screened_rates(k_Ne18_He4_to_Mg22_approx)*Y(He4)*Y(Ne18)*state.rho + screened_rates(k_Mg22_to_He4_Ne18_derived)*Y(Mg22)) +
        0.5*screened_rates(k_He3_He3_to_p_p_He4_reaclib)*amrex::Math::powi<2>(Y(He3))*state.rho +
        screened_rates(k_He3_p_to_He4_beta_pos_reaclib)*Y(He3)*Y(H1)*state.rho +
        screened_rates(k_He3_He4_p_to_He4_He4_beta_pos_modified)*Y(He3)*Y(He4)*state.rho;

    ydot_nuc(C12) =
        (-screened_rates(k_C12_p_to_N13_reaclib)*Y(C12)*Y(H1)*state.rho + screened_rates(k_N13_to_p_C12_derived)*Y(N13)) +
        (-screened_rates(k_C12_He4_to_O16_reaclib)*Y(C12)*Y(He4)*state.rho + screened_rates(k_O16_to_He4_C12_derived)*Y(O16)) +
        (screened_rates(k_N15_p_to_He4_C12_reaclib)*Y(N15)*Y(H1)*state.rho + -screened_rates(k_C12_He4_to_p_N15_derived)*Y(C12)*Y(He4)*state.rho) +
        (0.16666666666666667*screened_rates(k_He4_He4_He4_to_C12_reaclib)*amrex::Math::powi<3>(Y(He4))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_C12_to_3He4_derived)*Y(C12));

    ydot_nuc(C13) =
        screened_rates(k_N13_to_C13_beta_pos_reaclib)*Y(N13) +
        (-screened_rates(k_C13_p_to_N14_reaclib)*Y(C13)*Y(H1)*state.rho + screened_rates(k_N14_to_p_C13_derived)*Y(N14));

    ydot_nuc(N13) =
        -screened_rates(k_N13_to_C13_beta_pos_reaclib)*Y(N13) +
        (screened_rates(k_C12_p_to_N13_reaclib)*Y(C12)*Y(H1)*state.rho + -screened_rates(k_N13_to_p_C12_derived)*Y(N13)) +
        (-screened_rates(k_N13_p_to_O14_reaclib)*Y(N13)*Y(H1)*state.rho + screened_rates(k_O14_to_p_N13_derived)*Y(O14)) +
        (-screened_rates(k_N13_He4_to_p_O16_reaclib)*Y(He4)*Y(N13)*state.rho + screened_rates(k_O16_p_to_He4_N13_derived)*Y(O16)*Y(H1)*state.rho);

    ydot_nuc(N14) =
        screened_rates(k_O14_to_N14_beta_pos_reaclib)*Y(O14) +
        (screened_rates(k_C13_p_to_N14_reaclib)*Y(C13)*Y(H1)*state.rho + -screened_rates(k_N14_to_p_C13_derived)*Y(N14)) +
        (-screened_rates(k_N14_p_to_O15_reaclib)*Y(N14)*Y(H1)*state.rho + screened_rates(k_O15_to_p_N14_derived)*Y(O15)) +
        (-screened_rates(k_N14_He4_to_F18_reaclib)*Y(He4)*Y(N14)*state.rho + screened_rates(k_F18_to_He4_N14_derived)*Y(F18)) +
        (screened_rates(k_O17_p_to_He4_N14_reaclib)*Y(O17)*Y(H1)*state.rho + -screened_rates(k_N14_He4_to_p_O17_derived)*Y(He4)*Y(N14)*state.rho);

    ydot_nuc(N15) =
        screened_rates(k_O15_to_N15_beta_pos_reaclib)*Y(O15) +
        (-screened_rates(k_N15_p_to_O16_reaclib)*Y(N15)*Y(H1)*state.rho + screened_rates(k_O16_to_p_N15_derived)*Y(O16)) +
        (-screened_rates(k_N15_He4_to_F19_reaclib)*Y(He4)*Y(N15)*state.rho + screened_rates(k_F19_to_He4_N15_derived)*Y(F19)) +
        (-screened_rates(k_N15_p_to_He4_C12_reaclib)*Y(N15)*Y(H1)*state.rho + screened_rates(k_C12_He4_to_p_N15_derived)*Y(C12)*Y(He4)*state.rho) +
        (screened_rates(k_O18_p_to_He4_N15_reaclib)*Y(O18)*Y(H1)*state.rho + -screened_rates(k_N15_He4_to_p_O18_derived)*Y(He4)*Y(N15)*state.rho);

    ydot_nuc(O14) =
        -screened_rates(k_O14_to_N14_beta_pos_reaclib)*Y(O14) +
        (screened_rates(k_N13_p_to_O14_reaclib)*Y(N13)*Y(H1)*state.rho + -screened_rates(k_O14_to_p_N13_derived)*Y(O14)) +
        (-screened_rates(k_O14_He4_to_Ne18_reaclib)*Y(He4)*Y(O14)*state.rho + screened_rates(k_Ne18_to_He4_O14_derived)*Y(Ne18)) +
        (-screened_rates(k_O14_He4_to_p_F17_reaclib)*Y(He4)*Y(O14)*state.rho + screened_rates(k_F17_p_to_He4_O14_derived)*Y(F17)*Y(H1)*state.rho);

    ydot_nuc(O15) =
        -screened_rates(k_O15_to_N15_beta_pos_reaclib)*Y(O15) +
        (screened_rates(k_N14_p_to_O15_reaclib)*Y(N14)*Y(H1)*state.rho + -screened_rates(k_O15_to_p_N14_derived)*Y(O15)) +
        (-screened_rates(k_O15_He4_to_Ne19_reaclib)*Y(He4)*Y(O15)*state.rho + screened_rates(k_Ne19_to_He4_O15_derived)*Y(Ne19)) +
        (screened_rates(k_F18_p_to_He4_O15_reaclib)*Y(F18)*Y(H1)*state.rho + -screened_rates(k_O15_He4_to_p_F18_derived)*Y(He4)*Y(O15)*state.rho);

    ydot_nuc(O16) =
        (screened_rates(k_C12_He4_to_O16_reaclib)*Y(C12)*Y(He4)*state.rho + -screened_rates(k_O16_to_He4_C12_derived)*Y(O16)) +
        (screened_rates(k_N15_p_to_O16_reaclib)*Y(N15)*Y(H1)*state.rho + -screened_rates(k_O16_to_p_N15_derived)*Y(O16)) +
        (-screened_rates(k_O16_p_to_F17_reaclib)*Y(O16)*Y(H1)*state.rho + screened_rates(k_F17_to_p_O16_derived)*Y(F17)) +
        (-screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(He4)*Y(O16)*state.rho + screened_rates(k_Ne20_to_He4_O16_derived)*Y(Ne20)) +
        (screened_rates(k_N13_He4_to_p_O16_reaclib)*Y(He4)*Y(N13)*state.rho + -screened_rates(k_O16_p_to_He4_N13_derived)*Y(O16)*Y(H1)*state.rho) +
        (screened_rates(k_F19_p_to_He4_O16_reaclib)*Y(F19)*Y(H1)*state.rho + -screened_rates(k_O16_He4_to_p_F19_derived)*Y(He4)*Y(O16)*state.rho);

    ydot_nuc(O17) =
        (-screened_rates(k_O17_p_to_F18_reaclib)*Y(O17)*Y(H1)*state.rho + screened_rates(k_F18_to_p_O17_derived)*Y(F18)) +
        (-screened_rates(k_O17_p_to_He4_N14_reaclib)*Y(O17)*Y(H1)*state.rho + screened_rates(k_N14_He4_to_p_O17_derived)*Y(He4)*Y(N14)*state.rho) +
        (screened_rates(k_F17_to_O17_electron_capture_weaktab)*Y(F17) + -screened_rates(k_O17_to_F17_beta_neg_weaktab)*Y(O17));

    ydot_nuc(O18) =
        (-screened_rates(k_O18_p_to_F19_reaclib)*Y(O18)*Y(H1)*state.rho + screened_rates(k_F19_to_p_O18_derived)*Y(F19)) +
        (-screened_rates(k_O18_p_to_He4_N15_reaclib)*Y(O18)*Y(H1)*state.rho + screened_rates(k_N15_He4_to_p_O18_derived)*Y(He4)*Y(N15)*state.rho) +
        (screened_rates(k_F18_to_O18_electron_capture_weaktab)*Y(F18) + -screened_rates(k_O18_to_F18_beta_neg_weaktab)*Y(O18));

    ydot_nuc(F17) =
        (screened_rates(k_O16_p_to_F17_reaclib)*Y(O16)*Y(H1)*state.rho + -screened_rates(k_F17_to_p_O16_derived)*Y(F17)) +
        (-screened_rates(k_F17_p_to_Ne18_reaclib)*Y(F17)*Y(H1)*state.rho + screened_rates(k_Ne18_to_p_F17_derived)*Y(Ne18)) +
        (screened_rates(k_O14_He4_to_p_F17_reaclib)*Y(He4)*Y(O14)*state.rho + -screened_rates(k_F17_p_to_He4_O14_derived)*Y(F17)*Y(H1)*state.rho) +
        (-screened_rates(k_F17_He4_to_p_Ne20_derived)*Y(F17)*Y(He4)*state.rho + screened_rates(k_Ne20_p_to_He4_F17_reaclib)*Y(Ne20)*Y(H1)*state.rho) +
        (-screened_rates(k_F17_to_O17_electron_capture_weaktab)*Y(F17) + screened_rates(k_O17_to_F17_beta_neg_weaktab)*Y(O17));

    ydot_nuc(F18) =
        (screened_rates(k_N14_He4_to_F18_reaclib)*Y(He4)*Y(N14)*state.rho + -screened_rates(k_F18_to_He4_N14_derived)*Y(F18)) +
        (screened_rates(k_O17_p_to_F18_reaclib)*Y(O17)*Y(H1)*state.rho + -screened_rates(k_F18_to_p_O17_derived)*Y(F18)) +
        (-screened_rates(k_F18_p_to_Ne19_reaclib)*Y(F18)*Y(H1)*state.rho + screened_rates(k_Ne19_to_p_F18_derived)*Y(Ne19)) +
        (-screened_rates(k_F18_p_to_He4_O15_reaclib)*Y(F18)*Y(H1)*state.rho + screened_rates(k_O15_He4_to_p_F18_derived)*Y(He4)*Y(O15)*state.rho) +
        (-screened_rates(k_F18_to_O18_electron_capture_weaktab)*Y(F18) + screened_rates(k_O18_to_F18_beta_neg_weaktab)*Y(O18)) +
        (screened_rates(k_Ne18_to_F18_electron_capture_weaktab)*Y(Ne18) + -screened_rates(k_F18_to_Ne18_beta_neg_weaktab)*Y(F18));

    ydot_nuc(F19) =
        (screened_rates(k_N15_He4_to_F19_reaclib)*Y(He4)*Y(N15)*state.rho + -screened_rates(k_F19_to_He4_N15_derived)*Y(F19)) +
        (screened_rates(k_O18_p_to_F19_reaclib)*Y(O18)*Y(H1)*state.rho + -screened_rates(k_F19_to_p_O18_derived)*Y(F19)) +
        (-screened_rates(k_F19_p_to_Ne20_reaclib)*Y(F19)*Y(H1)*state.rho + screened_rates(k_Ne20_to_p_F19_derived)*Y(Ne20)) +
        (-screened_rates(k_F19_p_to_He4_O16_reaclib)*Y(F19)*Y(H1)*state.rho + screened_rates(k_O16_He4_to_p_F19_derived)*Y(He4)*Y(O16)*state.rho) +
        (screened_rates(k_Ne19_to_F19_electron_capture_weaktab)*Y(Ne19) + -screened_rates(k_F19_to_Ne19_beta_neg_weaktab)*Y(F19));

    ydot_nuc(Ne18) =
        (screened_rates(k_O14_He4_to_Ne18_reaclib)*Y(He4)*Y(O14)*state.rho + -screened_rates(k_Ne18_to_He4_O14_derived)*Y(Ne18)) +
        (screened_rates(k_F17_p_to_Ne18_reaclib)*Y(F17)*Y(H1)*state.rho + -screened_rates(k_Ne18_to_p_F17_derived)*Y(Ne18)) +
        (-screened_rates(k_Ne18_He4_to_Mg22_approx)*Y(He4)*Y(Ne18)*state.rho + screened_rates(k_Mg22_to_He4_Ne18_derived)*Y(Mg22)) +
        (-screened_rates(k_Ne18_to_F18_electron_capture_weaktab)*Y(Ne18) + screened_rates(k_F18_to_Ne18_beta_neg_weaktab)*Y(F18));

    ydot_nuc(Ne19) =
        (screened_rates(k_O15_He4_to_Ne19_reaclib)*Y(He4)*Y(O15)*state.rho + -screened_rates(k_Ne19_to_He4_O15_derived)*Y(Ne19)) +
        (screened_rates(k_F18_p_to_Ne19_reaclib)*Y(F18)*Y(H1)*state.rho + -screened_rates(k_Ne19_to_p_F18_derived)*Y(Ne19)) +
        -screened_rates(k_Ne19_3p_to_Mg22_beta_pos_modified)*Y(Ne19)*Y(H1)*state.rho +
        (-screened_rates(k_Ne19_to_F19_electron_capture_weaktab)*Y(Ne19) + screened_rates(k_F19_to_Ne19_beta_neg_weaktab)*Y(F19));

    ydot_nuc(Ne20) =
        (screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(He4)*Y(O16)*state.rho + -screened_rates(k_Ne20_to_He4_O16_derived)*Y(Ne20)) +
        (screened_rates(k_F19_p_to_Ne20_reaclib)*Y(F19)*Y(H1)*state.rho + -screened_rates(k_Ne20_to_p_F19_derived)*Y(Ne20)) +
        (-screened_rates(k_Ne20_He4_to_Mg24_reaclib)*Y(He4)*Y(Ne20)*state.rho + screened_rates(k_Mg24_to_He4_Ne20_derived)*Y(Mg24)) +
        (screened_rates(k_F17_He4_to_p_Ne20_derived)*Y(F17)*Y(He4)*state.rho + -screened_rates(k_Ne20_p_to_He4_F17_reaclib)*Y(Ne20)*Y(H1)*state.rho) +
        -screened_rates(k_Ne20_p_p_to_Mg22_modified)*Y(Ne20)*Y(H1)*state.rho;

    ydot_nuc(Mg22) =
        (screened_rates(k_Ne18_He4_to_Mg22_approx)*Y(He4)*Y(Ne18)*state.rho + -screened_rates(k_Mg22_to_He4_Ne18_derived)*Y(Mg22)) +
        screened_rates(k_Ne19_3p_to_Mg22_beta_pos_modified)*Y(Ne19)*Y(H1)*state.rho +
        screened_rates(k_Ne20_p_p_to_Mg22_modified)*Y(Ne20)*Y(H1)*state.rho;

    ydot_nuc(Mg24) =
        (screened_rates(k_Ne20_He4_to_Mg24_reaclib)*Y(He4)*Y(Ne20)*state.rho + -screened_rates(k_Mg24_to_He4_Ne20_derived)*Y(Mg24));

    ydot_nuc(Fe56) = 0.0_rt;

}


AMREX_GPU_DEVICE AMREX_INLINE
void actual_rhs (const burn_t& state, amrex::Array1D<amrex::Real, 1, neqs>& ydot)
{

    // Set molar abundances
    amrex::Array1D<amrex::Real, 1, NumSpec> Y;
    for (int i = 1; i <= NumSpec; ++i) {
        Y(i) = state.xn[i-1] * aion_inv[i-1];
    }

    // build the rates

    rate_t rate_eval;

    constexpr int do_T_derivatives = 0;

    evaluate_rates<do_T_derivatives, rate_t>(state, Y, rate_eval);

    rhs_nuc(state, ydot, Y, rate_eval.screened_rates);

    // ion binding energy contributions

    amrex::Real enuc;
    ener_gener_rate(ydot, enuc);

    // include any weak rate neutrino losses
    enuc += rate_eval.enuc_weak;

    // Get the thermal neutrino losses

    amrex::Real sneut{};
#ifdef NEUTRINOS
    constexpr int do_derivatives{0};
    amrex::Real dsneutdt{}, dsneutdd{}, dsnuda{}, dsnudz{};

    neutrino_cooling<do_derivatives>(state.T, state.rho, state.abar, state.zbar,
                                     sneut, dsneutdt, dsneutdd, dsnuda, dsnudz);
#endif

    // Append the energy equation (this is erg/g/s)

    ydot(net_ienuc) = enuc - sneut;

}


template<class MatrixType>
AMREX_GPU_DEVICE AMREX_INLINE
void jac_nuc(const burn_t& state,
             MatrixType& jac,
             const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
             const amrex::Array1D<amrex::Real, 1, Rates::NumRates>& screened_rates)
{

    using namespace Rates;

    amrex::Real scratch;

    scratch = -screened_rates(k_C12_p_to_N13_reaclib)*Y(C12)*state.rho - screened_rates(k_C13_p_to_N14_reaclib)*Y(C13)*state.rho - screened_rates(k_F17_p_to_He4_O14_derived)*Y(F17)*state.rho - screened_rates(k_F17_p_to_Ne18_reaclib)*Y(F17)*state.rho - screened_rates(k_F18_p_to_He4_O15_reaclib)*Y(F18)*state.rho - screened_rates(k_F18_p_to_Ne19_reaclib)*Y(F18)*state.rho - screened_rates(k_F19_p_to_He4_O16_reaclib)*Y(F19)*state.rho - screened_rates(k_F19_p_to_Ne20_reaclib)*Y(F19)*state.rho - screened_rates(k_He3_p_to_He4_beta_pos_reaclib)*Y(He3)*state.rho - screened_rates(k_N13_p_to_O14_reaclib)*Y(N13)*state.rho - screened_rates(k_N14_p_to_O15_reaclib)*Y(N14)*state.rho - screened_rates(k_N15_p_to_He4_C12_reaclib)*Y(N15)*state.rho - screened_rates(k_N15_p_to_O16_reaclib)*Y(N15)*state.rho - 3.0*screened_rates(k_Ne19_3p_to_Mg22_beta_pos_modified)*Y(Ne19)*state.rho - 2.0*screened_rates(k_Ne20_p_p_to_Mg22_modified)*Y(Ne20)*state.rho - screened_rates(k_Ne20_p_to_He4_F17_reaclib)*Y(Ne20)*state.rho - screened_rates(k_O16_p_to_F17_reaclib)*Y(O16)*state.rho - screened_rates(k_O16_p_to_He4_N13_derived)*Y(O16)*state.rho - screened_rates(k_O17_p_to_F18_reaclib)*Y(O17)*state.rho - screened_rates(k_O17_p_to_He4_N14_reaclib)*Y(O17)*state.rho - screened_rates(k_O18_p_to_F19_reaclib)*Y(O18)*state.rho - screened_rates(k_O18_p_to_He4_N15_reaclib)*Y(O18)*state.rho - 3.0*screened_rates(k_p_p_p_to_He3_beta_pos_modified)*Y(H1)*state.rho - 3.0*screened_rates(k_p_p_p_to_He3_electron_capture_modified)*Y(H1)*state.rho;
    jac.set(H1, H1, scratch);

    scratch = 2.0*screened_rates(k_He3_He3_to_p_p_He4_reaclib)*Y(He3)*state.rho - screened_rates(k_He3_He4_p_to_He4_He4_beta_pos_modified)*Y(He4)*state.rho - screened_rates(k_He3_p_to_He4_beta_pos_reaclib)*Y(H1)*state.rho;
    jac.set(H1, He3, scratch);

    scratch = screened_rates(k_C12_He4_to_p_N15_derived)*Y(C12)*state.rho + screened_rates(k_F17_He4_to_p_Ne20_derived)*Y(F17)*state.rho - screened_rates(k_He3_He4_p_to_He4_He4_beta_pos_modified)*Y(He3)*state.rho + screened_rates(k_N13_He4_to_p_O16_reaclib)*Y(N13)*state.rho + screened_rates(k_N14_He4_to_p_O17_derived)*Y(N14)*state.rho + screened_rates(k_N15_He4_to_p_O18_derived)*Y(N15)*state.rho + screened_rates(k_O14_He4_to_p_F17_reaclib)*Y(O14)*state.rho + screened_rates(k_O15_He4_to_p_F18_derived)*Y(O15)*state.rho + screened_rates(k_O16_He4_to_p_F19_derived)*Y(O16)*state.rho;
    jac.set(H1, He4, scratch);

    scratch = screened_rates(k_C12_He4_to_p_N15_derived)*Y(He4)*state.rho - screened_rates(k_C12_p_to_N13_reaclib)*Y(H1)*state.rho;
    jac.set(H1, C12, scratch);

    scratch = -screened_rates(k_C13_p_to_N14_reaclib)*Y(H1)*state.rho;
    jac.set(H1, C13, scratch);

    scratch = screened_rates(k_N13_He4_to_p_O16_reaclib)*Y(He4)*state.rho - screened_rates(k_N13_p_to_O14_reaclib)*Y(H1)*state.rho + screened_rates(k_N13_to_p_C12_derived);
    jac.set(H1, N13, scratch);

    scratch = screened_rates(k_N14_He4_to_p_O17_derived)*Y(He4)*state.rho - screened_rates(k_N14_p_to_O15_reaclib)*Y(H1)*state.rho + screened_rates(k_N14_to_p_C13_derived);
    jac.set(H1, N14, scratch);

    scratch = screened_rates(k_N15_He4_to_p_O18_derived)*Y(He4)*state.rho - screened_rates(k_N15_p_to_He4_C12_reaclib)*Y(H1)*state.rho - screened_rates(k_N15_p_to_O16_reaclib)*Y(H1)*state.rho;
    jac.set(H1, N15, scratch);

    scratch = screened_rates(k_O14_He4_to_p_F17_reaclib)*Y(He4)*state.rho + screened_rates(k_O14_to_p_N13_derived);
    jac.set(H1, O14, scratch);

    scratch = screened_rates(k_O15_He4_to_p_F18_derived)*Y(He4)*state.rho + screened_rates(k_O15_to_p_N14_derived);
    jac.set(H1, O15, scratch);

    scratch = screened_rates(k_O16_He4_to_p_F19_derived)*Y(He4)*state.rho - screened_rates(k_O16_p_to_F17_reaclib)*Y(H1)*state.rho - screened_rates(k_O16_p_to_He4_N13_derived)*Y(H1)*state.rho + screened_rates(k_O16_to_p_N15_derived);
    jac.set(H1, O16, scratch);

    scratch = -screened_rates(k_O17_p_to_F18_reaclib)*Y(H1)*state.rho - screened_rates(k_O17_p_to_He4_N14_reaclib)*Y(H1)*state.rho;
    jac.set(H1, O17, scratch);

    scratch = -screened_rates(k_O18_p_to_F19_reaclib)*Y(H1)*state.rho - screened_rates(k_O18_p_to_He4_N15_reaclib)*Y(H1)*state.rho;
    jac.set(H1, O18, scratch);

    scratch = screened_rates(k_F17_He4_to_p_Ne20_derived)*Y(He4)*state.rho - screened_rates(k_F17_p_to_He4_O14_derived)*Y(H1)*state.rho - screened_rates(k_F17_p_to_Ne18_reaclib)*Y(H1)*state.rho + screened_rates(k_F17_to_p_O16_derived);
    jac.set(H1, F17, scratch);

    scratch = -screened_rates(k_F18_p_to_He4_O15_reaclib)*Y(H1)*state.rho - screened_rates(k_F18_p_to_Ne19_reaclib)*Y(H1)*state.rho + screened_rates(k_F18_to_p_O17_derived);
    jac.set(H1, F18, scratch);

    scratch = -screened_rates(k_F19_p_to_He4_O16_reaclib)*Y(H1)*state.rho - screened_rates(k_F19_p_to_Ne20_reaclib)*Y(H1)*state.rho + screened_rates(k_F19_to_p_O18_derived);
    jac.set(H1, F19, scratch);

    scratch = screened_rates(k_Ne18_to_p_F17_derived);
    jac.set(H1, Ne18, scratch);

    scratch = -3.0*screened_rates(k_Ne19_3p_to_Mg22_beta_pos_modified)*Y(H1)*state.rho + screened_rates(k_Ne19_to_p_F18_derived);
    jac.set(H1, Ne19, scratch);

    scratch = -2.0*screened_rates(k_Ne20_p_p_to_Mg22_modified)*Y(H1)*state.rho - screened_rates(k_Ne20_p_to_He4_F17_reaclib)*Y(H1)*state.rho + screened_rates(k_Ne20_to_p_F19_derived);
    jac.set(H1, Ne20, scratch);

    jac.set(H1, Mg22, 0.0);

    jac.set(H1, Mg24, 0.0);

    jac.set(H1, Fe56, 0.0);

    scratch = -screened_rates(k_He3_p_to_He4_beta_pos_reaclib)*Y(He3)*state.rho + 1.0*screened_rates(k_p_p_p_to_He3_beta_pos_modified)*Y(H1)*state.rho + 1.0*screened_rates(k_p_p_p_to_He3_electron_capture_modified)*Y(H1)*state.rho;
    jac.set(He3, H1, scratch);

    scratch = -2.0*screened_rates(k_He3_He3_to_p_p_He4_reaclib)*Y(He3)*state.rho - screened_rates(k_He3_He4_p_to_He4_He4_beta_pos_modified)*Y(He4)*state.rho - screened_rates(k_He3_p_to_He4_beta_pos_reaclib)*Y(H1)*state.rho;
    jac.set(He3, He3, scratch);

    scratch = -screened_rates(k_He3_He4_p_to_He4_He4_beta_pos_modified)*Y(He3)*state.rho;
    jac.set(He3, He4, scratch);

    jac.set(He3, C12, 0.0);

    jac.set(He3, C13, 0.0);

    jac.set(He3, N13, 0.0);

    jac.set(He3, N14, 0.0);

    jac.set(He3, N15, 0.0);

    jac.set(He3, O14, 0.0);

    jac.set(He3, O15, 0.0);

    jac.set(He3, O16, 0.0);

    jac.set(He3, O17, 0.0);

    jac.set(He3, O18, 0.0);

    jac.set(He3, F17, 0.0);

    jac.set(He3, F18, 0.0);

    jac.set(He3, F19, 0.0);

    jac.set(He3, Ne18, 0.0);

    jac.set(He3, Ne19, 0.0);

    jac.set(He3, Ne20, 0.0);

    jac.set(He3, Mg22, 0.0);

    jac.set(He3, Mg24, 0.0);

    jac.set(He3, Fe56, 0.0);

    scratch = screened_rates(k_F17_p_to_He4_O14_derived)*Y(F17)*state.rho + screened_rates(k_F18_p_to_He4_O15_reaclib)*Y(F18)*state.rho + screened_rates(k_F19_p_to_He4_O16_reaclib)*Y(F19)*state.rho + screened_rates(k_He3_p_to_He4_beta_pos_reaclib)*Y(He3)*state.rho + screened_rates(k_N15_p_to_He4_C12_reaclib)*Y(N15)*state.rho + screened_rates(k_Ne20_p_to_He4_F17_reaclib)*Y(Ne20)*state.rho + screened_rates(k_O16_p_to_He4_N13_derived)*Y(O16)*state.rho + screened_rates(k_O17_p_to_He4_N14_reaclib)*Y(O17)*state.rho + screened_rates(k_O18_p_to_He4_N15_reaclib)*Y(O18)*state.rho;
    jac.set(He4, H1, scratch);

    scratch = 1.0*screened_rates(k_He3_He3_to_p_p_He4_reaclib)*Y(He3)*state.rho + screened_rates(k_He3_He4_p_to_He4_He4_beta_pos_modified)*Y(He4)*state.rho + screened_rates(k_He3_p_to_He4_beta_pos_reaclib)*Y(H1)*state.rho;
    jac.set(He4, He3, scratch);

    scratch = -screened_rates(k_C12_He4_to_O16_reaclib)*Y(C12)*state.rho - screened_rates(k_C12_He4_to_p_N15_derived)*Y(C12)*state.rho - screened_rates(k_F17_He4_to_p_Ne20_derived)*Y(F17)*state.rho + screened_rates(k_He3_He4_p_to_He4_He4_beta_pos_modified)*Y(He3)*state.rho - 1.5*screened_rates(k_He4_He4_He4_to_C12_reaclib)*amrex::Math::powi<2>(Y(He4))*amrex::Math::powi<2>(state.rho) - screened_rates(k_N13_He4_to_p_O16_reaclib)*Y(N13)*state.rho - screened_rates(k_N14_He4_to_F18_reaclib)*Y(N14)*state.rho - screened_rates(k_N14_He4_to_p_O17_derived)*Y(N14)*state.rho - screened_rates(k_N15_He4_to_F19_reaclib)*Y(N15)*state.rho - screened_rates(k_N15_He4_to_p_O18_derived)*Y(N15)*state.rho - screened_rates(k_Ne18_He4_to_Mg22_approx)*Y(Ne18)*state.rho - screened_rates(k_Ne20_He4_to_Mg24_reaclib)*Y(Ne20)*state.rho - screened_rates(k_O14_He4_to_Ne18_reaclib)*Y(O14)*state.rho - screened_rates(k_O14_He4_to_p_F17_reaclib)*Y(O14)*state.rho - screened_rates(k_O15_He4_to_Ne19_reaclib)*Y(O15)*state.rho - screened_rates(k_O15_He4_to_p_F18_derived)*Y(O15)*state.rho - screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(O16)*state.rho - screened_rates(k_O16_He4_to_p_F19_derived)*Y(O16)*state.rho;
    jac.set(He4, He4, scratch);

    scratch = -screened_rates(k_C12_He4_to_O16_reaclib)*Y(He4)*state.rho - screened_rates(k_C12_He4_to_p_N15_derived)*Y(He4)*state.rho + 3.0*screened_rates(k_C12_to_3He4_derived);
    jac.set(He4, C12, scratch);

    jac.set(He4, C13, 0.0);

    scratch = -screened_rates(k_N13_He4_to_p_O16_reaclib)*Y(He4)*state.rho;
    jac.set(He4, N13, scratch);

    scratch = -screened_rates(k_N14_He4_to_F18_reaclib)*Y(He4)*state.rho - screened_rates(k_N14_He4_to_p_O17_derived)*Y(He4)*state.rho;
    jac.set(He4, N14, scratch);

    scratch = -screened_rates(k_N15_He4_to_F19_reaclib)*Y(He4)*state.rho - screened_rates(k_N15_He4_to_p_O18_derived)*Y(He4)*state.rho + screened_rates(k_N15_p_to_He4_C12_reaclib)*Y(H1)*state.rho;
    jac.set(He4, N15, scratch);

    scratch = -screened_rates(k_O14_He4_to_Ne18_reaclib)*Y(He4)*state.rho - screened_rates(k_O14_He4_to_p_F17_reaclib)*Y(He4)*state.rho;
    jac.set(He4, O14, scratch);

    scratch = -screened_rates(k_O15_He4_to_Ne19_reaclib)*Y(He4)*state.rho - screened_rates(k_O15_He4_to_p_F18_derived)*Y(He4)*state.rho;
    jac.set(He4, O15, scratch);

    scratch = -screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(He4)*state.rho - screened_rates(k_O16_He4_to_p_F19_derived)*Y(He4)*state.rho + screened_rates(k_O16_p_to_He4_N13_derived)*Y(H1)*state.rho + screened_rates(k_O16_to_He4_C12_derived);
    jac.set(He4, O16, scratch);

    scratch = screened_rates(k_O17_p_to_He4_N14_reaclib)*Y(H1)*state.rho;
    jac.set(He4, O17, scratch);

    scratch = screened_rates(k_O18_p_to_He4_N15_reaclib)*Y(H1)*state.rho;
    jac.set(He4, O18, scratch);

    scratch = -screened_rates(k_F17_He4_to_p_Ne20_derived)*Y(He4)*state.rho + screened_rates(k_F17_p_to_He4_O14_derived)*Y(H1)*state.rho;
    jac.set(He4, F17, scratch);

    scratch = screened_rates(k_F18_p_to_He4_O15_reaclib)*Y(H1)*state.rho + screened_rates(k_F18_to_He4_N14_derived);
    jac.set(He4, F18, scratch);

    scratch = screened_rates(k_F19_p_to_He4_O16_reaclib)*Y(H1)*state.rho + screened_rates(k_F19_to_He4_N15_derived);
    jac.set(He4, F19, scratch);

    scratch = -screened_rates(k_Ne18_He4_to_Mg22_approx)*Y(He4)*state.rho + screened_rates(k_Ne18_to_He4_O14_derived);
    jac.set(He4, Ne18, scratch);

    scratch = screened_rates(k_Ne19_to_He4_O15_derived);
    jac.set(He4, Ne19, scratch);

    scratch = -screened_rates(k_Ne20_He4_to_Mg24_reaclib)*Y(He4)*state.rho + screened_rates(k_Ne20_p_to_He4_F17_reaclib)*Y(H1)*state.rho + screened_rates(k_Ne20_to_He4_O16_derived);
    jac.set(He4, Ne20, scratch);

    scratch = screened_rates(k_Mg22_to_He4_Ne18_derived);
    jac.set(He4, Mg22, scratch);

    scratch = screened_rates(k_Mg24_to_He4_Ne20_derived);
    jac.set(He4, Mg24, scratch);

    jac.set(He4, Fe56, 0.0);

    scratch = -screened_rates(k_C12_p_to_N13_reaclib)*Y(C12)*state.rho + screened_rates(k_N15_p_to_He4_C12_reaclib)*Y(N15)*state.rho;
    jac.set(C12, H1, scratch);

    jac.set(C12, He3, 0.0);

    scratch = -screened_rates(k_C12_He4_to_O16_reaclib)*Y(C12)*state.rho - screened_rates(k_C12_He4_to_p_N15_derived)*Y(C12)*state.rho + 0.5*screened_rates(k_He4_He4_He4_to_C12_reaclib)*amrex::Math::powi<2>(Y(He4))*amrex::Math::powi<2>(state.rho);
    jac.set(C12, He4, scratch);

    scratch = -screened_rates(k_C12_He4_to_O16_reaclib)*Y(He4)*state.rho - screened_rates(k_C12_He4_to_p_N15_derived)*Y(He4)*state.rho - screened_rates(k_C12_p_to_N13_reaclib)*Y(H1)*state.rho - screened_rates(k_C12_to_3He4_derived);
    jac.set(C12, C12, scratch);

    jac.set(C12, C13, 0.0);

    scratch = screened_rates(k_N13_to_p_C12_derived);
    jac.set(C12, N13, scratch);

    jac.set(C12, N14, 0.0);

    scratch = screened_rates(k_N15_p_to_He4_C12_reaclib)*Y(H1)*state.rho;
    jac.set(C12, N15, scratch);

    jac.set(C12, O14, 0.0);

    jac.set(C12, O15, 0.0);

    scratch = screened_rates(k_O16_to_He4_C12_derived);
    jac.set(C12, O16, scratch);

    jac.set(C12, O17, 0.0);

    jac.set(C12, O18, 0.0);

    jac.set(C12, F17, 0.0);

    jac.set(C12, F18, 0.0);

    jac.set(C12, F19, 0.0);

    jac.set(C12, Ne18, 0.0);

    jac.set(C12, Ne19, 0.0);

    jac.set(C12, Ne20, 0.0);

    jac.set(C12, Mg22, 0.0);

    jac.set(C12, Mg24, 0.0);

    jac.set(C12, Fe56, 0.0);

    scratch = -screened_rates(k_C13_p_to_N14_reaclib)*Y(C13)*state.rho;
    jac.set(C13, H1, scratch);

    jac.set(C13, He3, 0.0);

    jac.set(C13, He4, 0.0);

    jac.set(C13, C12, 0.0);

    scratch = -screened_rates(k_C13_p_to_N14_reaclib)*Y(H1)*state.rho;
    jac.set(C13, C13, scratch);

    scratch = screened_rates(k_N13_to_C13_beta_pos_reaclib);
    jac.set(C13, N13, scratch);

    scratch = screened_rates(k_N14_to_p_C13_derived);
    jac.set(C13, N14, scratch);

    jac.set(C13, N15, 0.0);

    jac.set(C13, O14, 0.0);

    jac.set(C13, O15, 0.0);

    jac.set(C13, O16, 0.0);

    jac.set(C13, O17, 0.0);

    jac.set(C13, O18, 0.0);

    jac.set(C13, F17, 0.0);

    jac.set(C13, F18, 0.0);

    jac.set(C13, F19, 0.0);

    jac.set(C13, Ne18, 0.0);

    jac.set(C13, Ne19, 0.0);

    jac.set(C13, Ne20, 0.0);

    jac.set(C13, Mg22, 0.0);

    jac.set(C13, Mg24, 0.0);

    jac.set(C13, Fe56, 0.0);

    scratch = screened_rates(k_C12_p_to_N13_reaclib)*Y(C12)*state.rho - screened_rates(k_N13_p_to_O14_reaclib)*Y(N13)*state.rho + screened_rates(k_O16_p_to_He4_N13_derived)*Y(O16)*state.rho;
    jac.set(N13, H1, scratch);

    jac.set(N13, He3, 0.0);

    scratch = -screened_rates(k_N13_He4_to_p_O16_reaclib)*Y(N13)*state.rho;
    jac.set(N13, He4, scratch);

    scratch = screened_rates(k_C12_p_to_N13_reaclib)*Y(H1)*state.rho;
    jac.set(N13, C12, scratch);

    jac.set(N13, C13, 0.0);

    scratch = -screened_rates(k_N13_He4_to_p_O16_reaclib)*Y(He4)*state.rho - screened_rates(k_N13_p_to_O14_reaclib)*Y(H1)*state.rho - screened_rates(k_N13_to_C13_beta_pos_reaclib) - screened_rates(k_N13_to_p_C12_derived);
    jac.set(N13, N13, scratch);

    jac.set(N13, N14, 0.0);

    jac.set(N13, N15, 0.0);

    scratch = screened_rates(k_O14_to_p_N13_derived);
    jac.set(N13, O14, scratch);

    jac.set(N13, O15, 0.0);

    scratch = screened_rates(k_O16_p_to_He4_N13_derived)*Y(H1)*state.rho;
    jac.set(N13, O16, scratch);

    jac.set(N13, O17, 0.0);

    jac.set(N13, O18, 0.0);

    jac.set(N13, F17, 0.0);

    jac.set(N13, F18, 0.0);

    jac.set(N13, F19, 0.0);

    jac.set(N13, Ne18, 0.0);

    jac.set(N13, Ne19, 0.0);

    jac.set(N13, Ne20, 0.0);

    jac.set(N13, Mg22, 0.0);

    jac.set(N13, Mg24, 0.0);

    jac.set(N13, Fe56, 0.0);

    scratch = screened_rates(k_C13_p_to_N14_reaclib)*Y(C13)*state.rho - screened_rates(k_N14_p_to_O15_reaclib)*Y(N14)*state.rho + screened_rates(k_O17_p_to_He4_N14_reaclib)*Y(O17)*state.rho;
    jac.set(N14, H1, scratch);

    jac.set(N14, He3, 0.0);

    scratch = -screened_rates(k_N14_He4_to_F18_reaclib)*Y(N14)*state.rho - screened_rates(k_N14_He4_to_p_O17_derived)*Y(N14)*state.rho;
    jac.set(N14, He4, scratch);

    jac.set(N14, C12, 0.0);

    scratch = screened_rates(k_C13_p_to_N14_reaclib)*Y(H1)*state.rho;
    jac.set(N14, C13, scratch);

    jac.set(N14, N13, 0.0);

    scratch = -screened_rates(k_N14_He4_to_F18_reaclib)*Y(He4)*state.rho - screened_rates(k_N14_He4_to_p_O17_derived)*Y(He4)*state.rho - screened_rates(k_N14_p_to_O15_reaclib)*Y(H1)*state.rho - screened_rates(k_N14_to_p_C13_derived);
    jac.set(N14, N14, scratch);

    jac.set(N14, N15, 0.0);

    scratch = screened_rates(k_O14_to_N14_beta_pos_reaclib);
    jac.set(N14, O14, scratch);

    scratch = screened_rates(k_O15_to_p_N14_derived);
    jac.set(N14, O15, scratch);

    jac.set(N14, O16, 0.0);

    scratch = screened_rates(k_O17_p_to_He4_N14_reaclib)*Y(H1)*state.rho;
    jac.set(N14, O17, scratch);

    jac.set(N14, O18, 0.0);

    jac.set(N14, F17, 0.0);

    scratch = screened_rates(k_F18_to_He4_N14_derived);
    jac.set(N14, F18, scratch);

    jac.set(N14, F19, 0.0);

    jac.set(N14, Ne18, 0.0);

    jac.set(N14, Ne19, 0.0);

    jac.set(N14, Ne20, 0.0);

    jac.set(N14, Mg22, 0.0);

    jac.set(N14, Mg24, 0.0);

    jac.set(N14, Fe56, 0.0);

    scratch = -screened_rates(k_N15_p_to_He4_C12_reaclib)*Y(N15)*state.rho - screened_rates(k_N15_p_to_O16_reaclib)*Y(N15)*state.rho + screened_rates(k_O18_p_to_He4_N15_reaclib)*Y(O18)*state.rho;
    jac.set(N15, H1, scratch);

    jac.set(N15, He3, 0.0);

    scratch = screened_rates(k_C12_He4_to_p_N15_derived)*Y(C12)*state.rho - screened_rates(k_N15_He4_to_F19_reaclib)*Y(N15)*state.rho - screened_rates(k_N15_He4_to_p_O18_derived)*Y(N15)*state.rho;
    jac.set(N15, He4, scratch);

    scratch = screened_rates(k_C12_He4_to_p_N15_derived)*Y(He4)*state.rho;
    jac.set(N15, C12, scratch);

    jac.set(N15, C13, 0.0);

    jac.set(N15, N13, 0.0);

    jac.set(N15, N14, 0.0);

    scratch = -screened_rates(k_N15_He4_to_F19_reaclib)*Y(He4)*state.rho - screened_rates(k_N15_He4_to_p_O18_derived)*Y(He4)*state.rho - screened_rates(k_N15_p_to_He4_C12_reaclib)*Y(H1)*state.rho - screened_rates(k_N15_p_to_O16_reaclib)*Y(H1)*state.rho;
    jac.set(N15, N15, scratch);

    jac.set(N15, O14, 0.0);

    scratch = screened_rates(k_O15_to_N15_beta_pos_reaclib);
    jac.set(N15, O15, scratch);

    scratch = screened_rates(k_O16_to_p_N15_derived);
    jac.set(N15, O16, scratch);

    jac.set(N15, O17, 0.0);

    scratch = screened_rates(k_O18_p_to_He4_N15_reaclib)*Y(H1)*state.rho;
    jac.set(N15, O18, scratch);

    jac.set(N15, F17, 0.0);

    jac.set(N15, F18, 0.0);

    scratch = screened_rates(k_F19_to_He4_N15_derived);
    jac.set(N15, F19, scratch);

    jac.set(N15, Ne18, 0.0);

    jac.set(N15, Ne19, 0.0);

    jac.set(N15, Ne20, 0.0);

    jac.set(N15, Mg22, 0.0);

    jac.set(N15, Mg24, 0.0);

    jac.set(N15, Fe56, 0.0);

    scratch = screened_rates(k_F17_p_to_He4_O14_derived)*Y(F17)*state.rho + screened_rates(k_N13_p_to_O14_reaclib)*Y(N13)*state.rho;
    jac.set(O14, H1, scratch);

    jac.set(O14, He3, 0.0);

    scratch = -screened_rates(k_O14_He4_to_Ne18_reaclib)*Y(O14)*state.rho - screened_rates(k_O14_He4_to_p_F17_reaclib)*Y(O14)*state.rho;
    jac.set(O14, He4, scratch);

    jac.set(O14, C12, 0.0);

    jac.set(O14, C13, 0.0);

    scratch = screened_rates(k_N13_p_to_O14_reaclib)*Y(H1)*state.rho;
    jac.set(O14, N13, scratch);

    jac.set(O14, N14, 0.0);

    jac.set(O14, N15, 0.0);

    scratch = -screened_rates(k_O14_He4_to_Ne18_reaclib)*Y(He4)*state.rho - screened_rates(k_O14_He4_to_p_F17_reaclib)*Y(He4)*state.rho - screened_rates(k_O14_to_N14_beta_pos_reaclib) - screened_rates(k_O14_to_p_N13_derived);
    jac.set(O14, O14, scratch);

    jac.set(O14, O15, 0.0);

    jac.set(O14, O16, 0.0);

    jac.set(O14, O17, 0.0);

    jac.set(O14, O18, 0.0);

    scratch = screened_rates(k_F17_p_to_He4_O14_derived)*Y(H1)*state.rho;
    jac.set(O14, F17, scratch);

    jac.set(O14, F18, 0.0);

    jac.set(O14, F19, 0.0);

    scratch = screened_rates(k_Ne18_to_He4_O14_derived);
    jac.set(O14, Ne18, scratch);

    jac.set(O14, Ne19, 0.0);

    jac.set(O14, Ne20, 0.0);

    jac.set(O14, Mg22, 0.0);

    jac.set(O14, Mg24, 0.0);

    jac.set(O14, Fe56, 0.0);

    scratch = screened_rates(k_F18_p_to_He4_O15_reaclib)*Y(F18)*state.rho + screened_rates(k_N14_p_to_O15_reaclib)*Y(N14)*state.rho;
    jac.set(O15, H1, scratch);

    jac.set(O15, He3, 0.0);

    scratch = -screened_rates(k_O15_He4_to_Ne19_reaclib)*Y(O15)*state.rho - screened_rates(k_O15_He4_to_p_F18_derived)*Y(O15)*state.rho;
    jac.set(O15, He4, scratch);

    jac.set(O15, C12, 0.0);

    jac.set(O15, C13, 0.0);

    jac.set(O15, N13, 0.0);

    scratch = screened_rates(k_N14_p_to_O15_reaclib)*Y(H1)*state.rho;
    jac.set(O15, N14, scratch);

    jac.set(O15, N15, 0.0);

    jac.set(O15, O14, 0.0);

    scratch = -screened_rates(k_O15_He4_to_Ne19_reaclib)*Y(He4)*state.rho - screened_rates(k_O15_He4_to_p_F18_derived)*Y(He4)*state.rho - screened_rates(k_O15_to_N15_beta_pos_reaclib) - screened_rates(k_O15_to_p_N14_derived);
    jac.set(O15, O15, scratch);

    jac.set(O15, O16, 0.0);

    jac.set(O15, O17, 0.0);

    jac.set(O15, O18, 0.0);

    jac.set(O15, F17, 0.0);

    scratch = screened_rates(k_F18_p_to_He4_O15_reaclib)*Y(H1)*state.rho;
    jac.set(O15, F18, scratch);

    jac.set(O15, F19, 0.0);

    jac.set(O15, Ne18, 0.0);

    scratch = screened_rates(k_Ne19_to_He4_O15_derived);
    jac.set(O15, Ne19, scratch);

    jac.set(O15, Ne20, 0.0);

    jac.set(O15, Mg22, 0.0);

    jac.set(O15, Mg24, 0.0);

    jac.set(O15, Fe56, 0.0);

    scratch = screened_rates(k_F19_p_to_He4_O16_reaclib)*Y(F19)*state.rho + screened_rates(k_N15_p_to_O16_reaclib)*Y(N15)*state.rho - screened_rates(k_O16_p_to_F17_reaclib)*Y(O16)*state.rho - screened_rates(k_O16_p_to_He4_N13_derived)*Y(O16)*state.rho;
    jac.set(O16, H1, scratch);

    jac.set(O16, He3, 0.0);

    scratch = screened_rates(k_C12_He4_to_O16_reaclib)*Y(C12)*state.rho + screened_rates(k_N13_He4_to_p_O16_reaclib)*Y(N13)*state.rho - screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(O16)*state.rho - screened_rates(k_O16_He4_to_p_F19_derived)*Y(O16)*state.rho;
    jac.set(O16, He4, scratch);

    scratch = screened_rates(k_C12_He4_to_O16_reaclib)*Y(He4)*state.rho;
    jac.set(O16, C12, scratch);

    jac.set(O16, C13, 0.0);

    scratch = screened_rates(k_N13_He4_to_p_O16_reaclib)*Y(He4)*state.rho;
    jac.set(O16, N13, scratch);

    jac.set(O16, N14, 0.0);

    scratch = screened_rates(k_N15_p_to_O16_reaclib)*Y(H1)*state.rho;
    jac.set(O16, N15, scratch);

    jac.set(O16, O14, 0.0);

    jac.set(O16, O15, 0.0);

    scratch = -screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(He4)*state.rho - screened_rates(k_O16_He4_to_p_F19_derived)*Y(He4)*state.rho - screened_rates(k_O16_p_to_F17_reaclib)*Y(H1)*state.rho - screened_rates(k_O16_p_to_He4_N13_derived)*Y(H1)*state.rho - screened_rates(k_O16_to_He4_C12_derived) - screened_rates(k_O16_to_p_N15_derived);
    jac.set(O16, O16, scratch);

    jac.set(O16, O17, 0.0);

    jac.set(O16, O18, 0.0);

    scratch = screened_rates(k_F17_to_p_O16_derived);
    jac.set(O16, F17, scratch);

    jac.set(O16, F18, 0.0);

    scratch = screened_rates(k_F19_p_to_He4_O16_reaclib)*Y(H1)*state.rho;
    jac.set(O16, F19, scratch);

    jac.set(O16, Ne18, 0.0);

    jac.set(O16, Ne19, 0.0);

    scratch = screened_rates(k_Ne20_to_He4_O16_derived);
    jac.set(O16, Ne20, scratch);

    jac.set(O16, Mg22, 0.0);

    jac.set(O16, Mg24, 0.0);

    jac.set(O16, Fe56, 0.0);

    scratch = -screened_rates(k_O17_p_to_F18_reaclib)*Y(O17)*state.rho - screened_rates(k_O17_p_to_He4_N14_reaclib)*Y(O17)*state.rho;
    jac.set(O17, H1, scratch);

    jac.set(O17, He3, 0.0);

    scratch = screened_rates(k_N14_He4_to_p_O17_derived)*Y(N14)*state.rho;
    jac.set(O17, He4, scratch);

    jac.set(O17, C12, 0.0);

    jac.set(O17, C13, 0.0);

    jac.set(O17, N13, 0.0);

    scratch = screened_rates(k_N14_He4_to_p_O17_derived)*Y(He4)*state.rho;
    jac.set(O17, N14, scratch);

    jac.set(O17, N15, 0.0);

    jac.set(O17, O14, 0.0);

    jac.set(O17, O15, 0.0);

    jac.set(O17, O16, 0.0);

    scratch = -screened_rates(k_O17_p_to_F18_reaclib)*Y(H1)*state.rho - screened_rates(k_O17_p_to_He4_N14_reaclib)*Y(H1)*state.rho - screened_rates(k_O17_to_F17_beta_neg_weaktab);
    jac.set(O17, O17, scratch);

    jac.set(O17, O18, 0.0);

    scratch = screened_rates(k_F17_to_O17_electron_capture_weaktab);
    jac.set(O17, F17, scratch);

    scratch = screened_rates(k_F18_to_p_O17_derived);
    jac.set(O17, F18, scratch);

    jac.set(O17, F19, 0.0);

    jac.set(O17, Ne18, 0.0);

    jac.set(O17, Ne19, 0.0);

    jac.set(O17, Ne20, 0.0);

    jac.set(O17, Mg22, 0.0);

    jac.set(O17, Mg24, 0.0);

    jac.set(O17, Fe56, 0.0);

    scratch = -screened_rates(k_O18_p_to_F19_reaclib)*Y(O18)*state.rho - screened_rates(k_O18_p_to_He4_N15_reaclib)*Y(O18)*state.rho;
    jac.set(O18, H1, scratch);

    jac.set(O18, He3, 0.0);

    scratch = screened_rates(k_N15_He4_to_p_O18_derived)*Y(N15)*state.rho;
    jac.set(O18, He4, scratch);

    jac.set(O18, C12, 0.0);

    jac.set(O18, C13, 0.0);

    jac.set(O18, N13, 0.0);

    jac.set(O18, N14, 0.0);

    scratch = screened_rates(k_N15_He4_to_p_O18_derived)*Y(He4)*state.rho;
    jac.set(O18, N15, scratch);

    jac.set(O18, O14, 0.0);

    jac.set(O18, O15, 0.0);

    jac.set(O18, O16, 0.0);

    jac.set(O18, O17, 0.0);

    scratch = -screened_rates(k_O18_p_to_F19_reaclib)*Y(H1)*state.rho - screened_rates(k_O18_p_to_He4_N15_reaclib)*Y(H1)*state.rho - screened_rates(k_O18_to_F18_beta_neg_weaktab);
    jac.set(O18, O18, scratch);

    jac.set(O18, F17, 0.0);

    scratch = screened_rates(k_F18_to_O18_electron_capture_weaktab);
    jac.set(O18, F18, scratch);

    scratch = screened_rates(k_F19_to_p_O18_derived);
    jac.set(O18, F19, scratch);

    jac.set(O18, Ne18, 0.0);

    jac.set(O18, Ne19, 0.0);

    jac.set(O18, Ne20, 0.0);

    jac.set(O18, Mg22, 0.0);

    jac.set(O18, Mg24, 0.0);

    jac.set(O18, Fe56, 0.0);

    scratch = -screened_rates(k_F17_p_to_He4_O14_derived)*Y(F17)*state.rho - screened_rates(k_F17_p_to_Ne18_reaclib)*Y(F17)*state.rho + screened_rates(k_Ne20_p_to_He4_F17_reaclib)*Y(Ne20)*state.rho + screened_rates(k_O16_p_to_F17_reaclib)*Y(O16)*state.rho;
    jac.set(F17, H1, scratch);

    jac.set(F17, He3, 0.0);

    scratch = -screened_rates(k_F17_He4_to_p_Ne20_derived)*Y(F17)*state.rho + screened_rates(k_O14_He4_to_p_F17_reaclib)*Y(O14)*state.rho;
    jac.set(F17, He4, scratch);

    jac.set(F17, C12, 0.0);

    jac.set(F17, C13, 0.0);

    jac.set(F17, N13, 0.0);

    jac.set(F17, N14, 0.0);

    jac.set(F17, N15, 0.0);

    scratch = screened_rates(k_O14_He4_to_p_F17_reaclib)*Y(He4)*state.rho;
    jac.set(F17, O14, scratch);

    jac.set(F17, O15, 0.0);

    scratch = screened_rates(k_O16_p_to_F17_reaclib)*Y(H1)*state.rho;
    jac.set(F17, O16, scratch);

    scratch = screened_rates(k_O17_to_F17_beta_neg_weaktab);
    jac.set(F17, O17, scratch);

    jac.set(F17, O18, 0.0);

    scratch = -screened_rates(k_F17_He4_to_p_Ne20_derived)*Y(He4)*state.rho - screened_rates(k_F17_p_to_He4_O14_derived)*Y(H1)*state.rho - screened_rates(k_F17_p_to_Ne18_reaclib)*Y(H1)*state.rho - screened_rates(k_F17_to_O17_electron_capture_weaktab) - screened_rates(k_F17_to_p_O16_derived);
    jac.set(F17, F17, scratch);

    jac.set(F17, F18, 0.0);

    jac.set(F17, F19, 0.0);

    scratch = screened_rates(k_Ne18_to_p_F17_derived);
    jac.set(F17, Ne18, scratch);

    jac.set(F17, Ne19, 0.0);

    scratch = screened_rates(k_Ne20_p_to_He4_F17_reaclib)*Y(H1)*state.rho;
    jac.set(F17, Ne20, scratch);

    jac.set(F17, Mg22, 0.0);

    jac.set(F17, Mg24, 0.0);

    jac.set(F17, Fe56, 0.0);

    scratch = -screened_rates(k_F18_p_to_He4_O15_reaclib)*Y(F18)*state.rho - screened_rates(k_F18_p_to_Ne19_reaclib)*Y(F18)*state.rho + screened_rates(k_O17_p_to_F18_reaclib)*Y(O17)*state.rho;
    jac.set(F18, H1, scratch);

    jac.set(F18, He3, 0.0);

    scratch = screened_rates(k_N14_He4_to_F18_reaclib)*Y(N14)*state.rho + screened_rates(k_O15_He4_to_p_F18_derived)*Y(O15)*state.rho;
    jac.set(F18, He4, scratch);

    jac.set(F18, C12, 0.0);

    jac.set(F18, C13, 0.0);

    jac.set(F18, N13, 0.0);

    scratch = screened_rates(k_N14_He4_to_F18_reaclib)*Y(He4)*state.rho;
    jac.set(F18, N14, scratch);

    jac.set(F18, N15, 0.0);

    jac.set(F18, O14, 0.0);

    scratch = screened_rates(k_O15_He4_to_p_F18_derived)*Y(He4)*state.rho;
    jac.set(F18, O15, scratch);

    jac.set(F18, O16, 0.0);

    scratch = screened_rates(k_O17_p_to_F18_reaclib)*Y(H1)*state.rho;
    jac.set(F18, O17, scratch);

    scratch = screened_rates(k_O18_to_F18_beta_neg_weaktab);
    jac.set(F18, O18, scratch);

    jac.set(F18, F17, 0.0);

    scratch = -screened_rates(k_F18_p_to_He4_O15_reaclib)*Y(H1)*state.rho - screened_rates(k_F18_p_to_Ne19_reaclib)*Y(H1)*state.rho - screened_rates(k_F18_to_He4_N14_derived) - screened_rates(k_F18_to_Ne18_beta_neg_weaktab) - screened_rates(k_F18_to_O18_electron_capture_weaktab) - screened_rates(k_F18_to_p_O17_derived);
    jac.set(F18, F18, scratch);

    jac.set(F18, F19, 0.0);

    scratch = screened_rates(k_Ne18_to_F18_electron_capture_weaktab);
    jac.set(F18, Ne18, scratch);

    scratch = screened_rates(k_Ne19_to_p_F18_derived);
    jac.set(F18, Ne19, scratch);

    jac.set(F18, Ne20, 0.0);

    jac.set(F18, Mg22, 0.0);

    jac.set(F18, Mg24, 0.0);

    jac.set(F18, Fe56, 0.0);

    scratch = -screened_rates(k_F19_p_to_He4_O16_reaclib)*Y(F19)*state.rho - screened_rates(k_F19_p_to_Ne20_reaclib)*Y(F19)*state.rho + screened_rates(k_O18_p_to_F19_reaclib)*Y(O18)*state.rho;
    jac.set(F19, H1, scratch);

    jac.set(F19, He3, 0.0);

    scratch = screened_rates(k_N15_He4_to_F19_reaclib)*Y(N15)*state.rho + screened_rates(k_O16_He4_to_p_F19_derived)*Y(O16)*state.rho;
    jac.set(F19, He4, scratch);

    jac.set(F19, C12, 0.0);

    jac.set(F19, C13, 0.0);

    jac.set(F19, N13, 0.0);

    jac.set(F19, N14, 0.0);

    scratch = screened_rates(k_N15_He4_to_F19_reaclib)*Y(He4)*state.rho;
    jac.set(F19, N15, scratch);

    jac.set(F19, O14, 0.0);

    jac.set(F19, O15, 0.0);

    scratch = screened_rates(k_O16_He4_to_p_F19_derived)*Y(He4)*state.rho;
    jac.set(F19, O16, scratch);

    jac.set(F19, O17, 0.0);

    scratch = screened_rates(k_O18_p_to_F19_reaclib)*Y(H1)*state.rho;
    jac.set(F19, O18, scratch);

    jac.set(F19, F17, 0.0);

    jac.set(F19, F18, 0.0);

    scratch = -screened_rates(k_F19_p_to_He4_O16_reaclib)*Y(H1)*state.rho - screened_rates(k_F19_p_to_Ne20_reaclib)*Y(H1)*state.rho - screened_rates(k_F19_to_He4_N15_derived) - screened_rates(k_F19_to_Ne19_beta_neg_weaktab) - screened_rates(k_F19_to_p_O18_derived);
    jac.set(F19, F19, scratch);

    jac.set(F19, Ne18, 0.0);

    scratch = screened_rates(k_Ne19_to_F19_electron_capture_weaktab);
    jac.set(F19, Ne19, scratch);

    scratch = screened_rates(k_Ne20_to_p_F19_derived);
    jac.set(F19, Ne20, scratch);

    jac.set(F19, Mg22, 0.0);

    jac.set(F19, Mg24, 0.0);

    jac.set(F19, Fe56, 0.0);

    scratch = screened_rates(k_F17_p_to_Ne18_reaclib)*Y(F17)*state.rho;
    jac.set(Ne18, H1, scratch);

    jac.set(Ne18, He3, 0.0);

    scratch = -screened_rates(k_Ne18_He4_to_Mg22_approx)*Y(Ne18)*state.rho + screened_rates(k_O14_He4_to_Ne18_reaclib)*Y(O14)*state.rho;
    jac.set(Ne18, He4, scratch);

    jac.set(Ne18, C12, 0.0);

    jac.set(Ne18, C13, 0.0);

    jac.set(Ne18, N13, 0.0);

    jac.set(Ne18, N14, 0.0);

    jac.set(Ne18, N15, 0.0);

    scratch = screened_rates(k_O14_He4_to_Ne18_reaclib)*Y(He4)*state.rho;
    jac.set(Ne18, O14, scratch);

    jac.set(Ne18, O15, 0.0);

    jac.set(Ne18, O16, 0.0);

    jac.set(Ne18, O17, 0.0);

    jac.set(Ne18, O18, 0.0);

    scratch = screened_rates(k_F17_p_to_Ne18_reaclib)*Y(H1)*state.rho;
    jac.set(Ne18, F17, scratch);

    scratch = screened_rates(k_F18_to_Ne18_beta_neg_weaktab);
    jac.set(Ne18, F18, scratch);

    jac.set(Ne18, F19, 0.0);

    scratch = -screened_rates(k_Ne18_He4_to_Mg22_approx)*Y(He4)*state.rho - screened_rates(k_Ne18_to_F18_electron_capture_weaktab) - screened_rates(k_Ne18_to_He4_O14_derived) - screened_rates(k_Ne18_to_p_F17_derived);
    jac.set(Ne18, Ne18, scratch);

    jac.set(Ne18, Ne19, 0.0);

    jac.set(Ne18, Ne20, 0.0);

    scratch = screened_rates(k_Mg22_to_He4_Ne18_derived);
    jac.set(Ne18, Mg22, scratch);

    jac.set(Ne18, Mg24, 0.0);

    jac.set(Ne18, Fe56, 0.0);

    scratch = screened_rates(k_F18_p_to_Ne19_reaclib)*Y(F18)*state.rho - screened_rates(k_Ne19_3p_to_Mg22_beta_pos_modified)*Y(Ne19)*state.rho;
    jac.set(Ne19, H1, scratch);

    jac.set(Ne19, He3, 0.0);

    scratch = screened_rates(k_O15_He4_to_Ne19_reaclib)*Y(O15)*state.rho;
    jac.set(Ne19, He4, scratch);

    jac.set(Ne19, C12, 0.0);

    jac.set(Ne19, C13, 0.0);

    jac.set(Ne19, N13, 0.0);

    jac.set(Ne19, N14, 0.0);

    jac.set(Ne19, N15, 0.0);

    jac.set(Ne19, O14, 0.0);

    scratch = screened_rates(k_O15_He4_to_Ne19_reaclib)*Y(He4)*state.rho;
    jac.set(Ne19, O15, scratch);

    jac.set(Ne19, O16, 0.0);

    jac.set(Ne19, O17, 0.0);

    jac.set(Ne19, O18, 0.0);

    jac.set(Ne19, F17, 0.0);

    scratch = screened_rates(k_F18_p_to_Ne19_reaclib)*Y(H1)*state.rho;
    jac.set(Ne19, F18, scratch);

    scratch = screened_rates(k_F19_to_Ne19_beta_neg_weaktab);
    jac.set(Ne19, F19, scratch);

    jac.set(Ne19, Ne18, 0.0);

    scratch = -screened_rates(k_Ne19_3p_to_Mg22_beta_pos_modified)*Y(H1)*state.rho - screened_rates(k_Ne19_to_F19_electron_capture_weaktab) - screened_rates(k_Ne19_to_He4_O15_derived) - screened_rates(k_Ne19_to_p_F18_derived);
    jac.set(Ne19, Ne19, scratch);

    jac.set(Ne19, Ne20, 0.0);

    jac.set(Ne19, Mg22, 0.0);

    jac.set(Ne19, Mg24, 0.0);

    jac.set(Ne19, Fe56, 0.0);

    scratch = screened_rates(k_F19_p_to_Ne20_reaclib)*Y(F19)*state.rho - screened_rates(k_Ne20_p_p_to_Mg22_modified)*Y(Ne20)*state.rho - screened_rates(k_Ne20_p_to_He4_F17_reaclib)*Y(Ne20)*state.rho;
    jac.set(Ne20, H1, scratch);

    jac.set(Ne20, He3, 0.0);

    scratch = screened_rates(k_F17_He4_to_p_Ne20_derived)*Y(F17)*state.rho - screened_rates(k_Ne20_He4_to_Mg24_reaclib)*Y(Ne20)*state.rho + screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(O16)*state.rho;
    jac.set(Ne20, He4, scratch);

    jac.set(Ne20, C12, 0.0);

    jac.set(Ne20, C13, 0.0);

    jac.set(Ne20, N13, 0.0);

    jac.set(Ne20, N14, 0.0);

    jac.set(Ne20, N15, 0.0);

    jac.set(Ne20, O14, 0.0);

    jac.set(Ne20, O15, 0.0);

    scratch = screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(He4)*state.rho;
    jac.set(Ne20, O16, scratch);

    jac.set(Ne20, O17, 0.0);

    jac.set(Ne20, O18, 0.0);

    scratch = screened_rates(k_F17_He4_to_p_Ne20_derived)*Y(He4)*state.rho;
    jac.set(Ne20, F17, scratch);

    jac.set(Ne20, F18, 0.0);

    scratch = screened_rates(k_F19_p_to_Ne20_reaclib)*Y(H1)*state.rho;
    jac.set(Ne20, F19, scratch);

    jac.set(Ne20, Ne18, 0.0);

    jac.set(Ne20, Ne19, 0.0);

    scratch = -screened_rates(k_Ne20_He4_to_Mg24_reaclib)*Y(He4)*state.rho - screened_rates(k_Ne20_p_p_to_Mg22_modified)*Y(H1)*state.rho - screened_rates(k_Ne20_p_to_He4_F17_reaclib)*Y(H1)*state.rho - screened_rates(k_Ne20_to_He4_O16_derived) - screened_rates(k_Ne20_to_p_F19_derived);
    jac.set(Ne20, Ne20, scratch);

    jac.set(Ne20, Mg22, 0.0);

    scratch = screened_rates(k_Mg24_to_He4_Ne20_derived);
    jac.set(Ne20, Mg24, scratch);

    jac.set(Ne20, Fe56, 0.0);

    scratch = screened_rates(k_Ne19_3p_to_Mg22_beta_pos_modified)*Y(Ne19)*state.rho + screened_rates(k_Ne20_p_p_to_Mg22_modified)*Y(Ne20)*state.rho;
    jac.set(Mg22, H1, scratch);

    jac.set(Mg22, He3, 0.0);

    scratch = screened_rates(k_Ne18_He4_to_Mg22_approx)*Y(Ne18)*state.rho;
    jac.set(Mg22, He4, scratch);

    jac.set(Mg22, C12, 0.0);

    jac.set(Mg22, C13, 0.0);

    jac.set(Mg22, N13, 0.0);

    jac.set(Mg22, N14, 0.0);

    jac.set(Mg22, N15, 0.0);

    jac.set(Mg22, O14, 0.0);

    jac.set(Mg22, O15, 0.0);

    jac.set(Mg22, O16, 0.0);

    jac.set(Mg22, O17, 0.0);

    jac.set(Mg22, O18, 0.0);

    jac.set(Mg22, F17, 0.0);

    jac.set(Mg22, F18, 0.0);

    jac.set(Mg22, F19, 0.0);

    scratch = screened_rates(k_Ne18_He4_to_Mg22_approx)*Y(He4)*state.rho;
    jac.set(Mg22, Ne18, scratch);

    scratch = screened_rates(k_Ne19_3p_to_Mg22_beta_pos_modified)*Y(H1)*state.rho;
    jac.set(Mg22, Ne19, scratch);

    scratch = screened_rates(k_Ne20_p_p_to_Mg22_modified)*Y(H1)*state.rho;
    jac.set(Mg22, Ne20, scratch);

    scratch = -screened_rates(k_Mg22_to_He4_Ne18_derived);
    jac.set(Mg22, Mg22, scratch);

    jac.set(Mg22, Mg24, 0.0);

    jac.set(Mg22, Fe56, 0.0);

    jac.set(Mg24, H1, 0.0);

    jac.set(Mg24, He3, 0.0);

    scratch = screened_rates(k_Ne20_He4_to_Mg24_reaclib)*Y(Ne20)*state.rho;
    jac.set(Mg24, He4, scratch);

    jac.set(Mg24, C12, 0.0);

    jac.set(Mg24, C13, 0.0);

    jac.set(Mg24, N13, 0.0);

    jac.set(Mg24, N14, 0.0);

    jac.set(Mg24, N15, 0.0);

    jac.set(Mg24, O14, 0.0);

    jac.set(Mg24, O15, 0.0);

    jac.set(Mg24, O16, 0.0);

    jac.set(Mg24, O17, 0.0);

    jac.set(Mg24, O18, 0.0);

    jac.set(Mg24, F17, 0.0);

    jac.set(Mg24, F18, 0.0);

    jac.set(Mg24, F19, 0.0);

    jac.set(Mg24, Ne18, 0.0);

    jac.set(Mg24, Ne19, 0.0);

    scratch = screened_rates(k_Ne20_He4_to_Mg24_reaclib)*Y(He4)*state.rho;
    jac.set(Mg24, Ne20, scratch);

    jac.set(Mg24, Mg22, 0.0);

    scratch = -screened_rates(k_Mg24_to_He4_Ne20_derived);
    jac.set(Mg24, Mg24, scratch);

    jac.set(Mg24, Fe56, 0.0);

    jac.set(Fe56, H1, 0.0);

    jac.set(Fe56, He3, 0.0);

    jac.set(Fe56, He4, 0.0);

    jac.set(Fe56, C12, 0.0);

    jac.set(Fe56, C13, 0.0);

    jac.set(Fe56, N13, 0.0);

    jac.set(Fe56, N14, 0.0);

    jac.set(Fe56, N15, 0.0);

    jac.set(Fe56, O14, 0.0);

    jac.set(Fe56, O15, 0.0);

    jac.set(Fe56, O16, 0.0);

    jac.set(Fe56, O17, 0.0);

    jac.set(Fe56, O18, 0.0);

    jac.set(Fe56, F17, 0.0);

    jac.set(Fe56, F18, 0.0);

    jac.set(Fe56, F19, 0.0);

    jac.set(Fe56, Ne18, 0.0);

    jac.set(Fe56, Ne19, 0.0);

    jac.set(Fe56, Ne20, 0.0);

    jac.set(Fe56, Mg22, 0.0);

    jac.set(Fe56, Mg24, 0.0);

    jac.set(Fe56, Fe56, 0.0);


}



template<class MatrixType>
AMREX_GPU_DEVICE AMREX_INLINE
void actual_jac(const burn_t& state, MatrixType& jac)
{

    // Set molar abundances
    amrex::Array1D<amrex::Real, 1, NumSpec> Y;
    for (int i = 1; i <= NumSpec; ++i) {
        Y(i) = state.xn[i-1] * aion_inv[i-1];
    }


    rate_derivs_t rate_eval;

    constexpr int do_T_derivatives = 1;

    evaluate_rates<do_T_derivatives, rate_derivs_t>(state, Y, rate_eval);

    // Species Jacobian elements with respect to other species

    jac_nuc(state, jac, Y, rate_eval.screened_rates);

    // Energy generation rate Jacobian elements with respect to species

    for (int j = 1; j <= NumSpec; ++j) {
        auto jac_slice_2 = [&](int i) -> amrex::Real { return jac.get(i, j); };
        amrex::Real tmp{};
        ener_gener_rate(jac_slice_2, tmp);
        jac(net_ienuc,j) = tmp;
    }

    // Account for the thermal neutrino losses

    amrex::Real dsneutdt{};
#ifdef NEUTRINOS
    amrex::Real sneut, dsneutdd, dsnuda, dsnudz;
    constexpr int do_derivatives{1};
    neutrino_cooling<do_derivatives>(state.T, state.rho, state.abar, state.zbar,
                                     sneut, dsneutdt, dsneutdd, dsnuda, dsnudz);

    for (int j = 1; j <= NumSpec; ++j) {
       amrex::Real b1 = (-state.abar * state.abar * dsnuda + (zion[j-1] - state.zbar) * state.abar * dsnudz);
       jac.add(net_ienuc, j, -b1);
    }
#endif

    // Evaluate the Jacobian elements with respect to energy by
    // calling the RHS using d(rate) / dT and then transform them
    // to our energy integration variable.

    amrex::Array1D<amrex::Real, 1, neqs>  yderivs;

    rhs_nuc(state, yderivs, Y, rate_eval.dscreened_rates_dT);

    for (int k = 1; k <= NumSpec; k++) {
        jac.set(k, net_ienuc, temperature_to_energy_jacobian(state, yderivs(k)));
    }


    // finally, d(de/dt)/de

    amrex::Real jac_e_T;
    ener_gener_rate(yderivs, jac_e_T);
    jac_e_T -= dsneutdt;
    jac.set(net_ienuc, net_ienuc, temperature_to_energy_jacobian(state, jac_e_T));

}


AMREX_INLINE
void actual_rhs_init () {


}


#endif
