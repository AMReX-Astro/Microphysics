#ifndef actual_rhs_H
#define actual_rhs_H

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
#include <reaclib_rates.H>
#include <table_rates.H>
#include <temperature_table_rates.H>
#include <derived_rates.H>

using namespace amrex::literals;
using namespace ArrayUtil;

using namespace Species;
using namespace Rates;

using namespace rate_tables;


template<class T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
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


template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void evaluate_rates(const burn_t& state,
                    const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
                    T& rate_eval) {

    // here Y is consistent with state.xn

    [[maybe_unused]] amrex::Real rhoy = state.rho * state.y_e;

    // Calculate Reaclib rates

    using number_t = std::conditional_t<do_T_derivatives, autodiff::dual, amrex::Real>;
    number_t temp = state.T;
    if constexpr (do_T_derivatives) {
        // seed the dual number for temperature before calculating anything with it
        autodiff::seed(temp);
    }

#ifdef SCREENING
    plasma_state_t<number_t> pstate{};
    fill_plasma_state(pstate, temp, state.rho, Y);
#endif

    const tf_t tfactors = evaluate_tfactors(state.T);

    fill_reaclib_rates<do_T_derivatives, T>(tfactors, rate_eval);

    temp_tabular::fill_rates<do_T_derivatives, T>(tfactors, rate_eval);

    // Calculate Derived Rates. This should go last but before approx rates.
    fill_derived_rates<do_T_derivatives, T>(tfactors, rate_eval);


#ifdef SCREENING
    // Evaluate screening factors

    amrex::Real ratraw, dratraw_dT;
    amrex::Real scor, dscor_dt;
    [[maybe_unused]] amrex::Real scor2, dscor2_dt;


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 6.0_rt, 12.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_p_C12_to_N13_reaclib);
    rate_eval.screened_rates(k_p_C12_to_N13_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_C12_to_N13_reaclib);
        rate_eval.dscreened_rates_dT(k_p_C12_to_N13_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 6.0_rt, 12.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_He4_C12_to_O16_reaclib);
    rate_eval.screened_rates(k_He4_C12_to_O16_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_C12_to_O16_reaclib);
        rate_eval.dscreened_rates_dT(k_He4_C12_to_O16_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_C12_to_p_N15_derived);
    rate_eval.screened_rates(k_He4_C12_to_p_N15_derived) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_C12_to_p_N15_derived);
        rate_eval.dscreened_rates_dT(k_He4_C12_to_p_N15_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 6.0_rt, 13.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_p_C13_to_N14_reaclib);
    rate_eval.screened_rates(k_p_C13_to_N14_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_C13_to_N14_reaclib);
        rate_eval.dscreened_rates_dT(k_p_C13_to_N14_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 7.0_rt, 13.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_p_N13_to_O14_reaclib);
    rate_eval.screened_rates(k_p_N13_to_O14_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_N13_to_O14_reaclib);
        rate_eval.dscreened_rates_dT(k_p_N13_to_O14_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 7.0_rt, 14.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_p_N14_to_O15_reaclib);
    rate_eval.screened_rates(k_p_N14_to_O15_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_N14_to_O15_reaclib);
        rate_eval.dscreened_rates_dT(k_p_N14_to_O15_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 7.0_rt, 14.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_He4_N14_to_F18_reaclib);
    rate_eval.screened_rates(k_He4_N14_to_F18_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_N14_to_F18_reaclib);
        rate_eval.dscreened_rates_dT(k_He4_N14_to_F18_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_N14_to_p_O17_derived);
    rate_eval.screened_rates(k_He4_N14_to_p_O17_derived) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_N14_to_p_O17_derived);
        rate_eval.dscreened_rates_dT(k_He4_N14_to_p_O17_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 7.0_rt, 15.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_p_N15_to_O16_reaclib);
    rate_eval.screened_rates(k_p_N15_to_O16_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_N15_to_O16_reaclib);
        rate_eval.dscreened_rates_dT(k_p_N15_to_O16_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_N15_to_He4_C12_reaclib);
    rate_eval.screened_rates(k_p_N15_to_He4_C12_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_N15_to_He4_C12_reaclib);
        rate_eval.dscreened_rates_dT(k_p_N15_to_He4_C12_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 7.0_rt, 15.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_He4_N15_to_F19_reaclib);
    rate_eval.screened_rates(k_He4_N15_to_F19_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_N15_to_F19_reaclib);
        rate_eval.dscreened_rates_dT(k_He4_N15_to_F19_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_N15_to_p_O18_derived);
    rate_eval.screened_rates(k_He4_N15_to_p_O18_derived) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_N15_to_p_O18_derived);
        rate_eval.dscreened_rates_dT(k_He4_N15_to_p_O18_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 8.0_rt, 14.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_He4_O14_to_Ne18_reaclib);
    rate_eval.screened_rates(k_He4_O14_to_Ne18_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_O14_to_Ne18_reaclib);
        rate_eval.dscreened_rates_dT(k_He4_O14_to_Ne18_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_O14_to_p_F17_reaclib);
    rate_eval.screened_rates(k_He4_O14_to_p_F17_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_O14_to_p_F17_reaclib);
        rate_eval.dscreened_rates_dT(k_He4_O14_to_p_F17_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 8.0_rt, 15.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_He4_O15_to_Ne19_reaclib);
    rate_eval.screened_rates(k_He4_O15_to_Ne19_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_O15_to_Ne19_reaclib);
        rate_eval.dscreened_rates_dT(k_He4_O15_to_Ne19_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_O15_to_p_F18_derived);
    rate_eval.screened_rates(k_He4_O15_to_p_F18_derived) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_O15_to_p_F18_derived);
        rate_eval.dscreened_rates_dT(k_He4_O15_to_p_F18_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 8.0_rt, 16.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_p_O16_to_F17_reaclib);
    rate_eval.screened_rates(k_p_O16_to_F17_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_O16_to_F17_reaclib);
        rate_eval.dscreened_rates_dT(k_p_O16_to_F17_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_O16_to_He4_N13_derived);
    rate_eval.screened_rates(k_p_O16_to_He4_N13_derived) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_O16_to_He4_N13_derived);
        rate_eval.dscreened_rates_dT(k_p_O16_to_He4_N13_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 8.0_rt, 16.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_He4_O16_to_Ne20_reaclib);
    rate_eval.screened_rates(k_He4_O16_to_Ne20_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_O16_to_Ne20_reaclib);
        rate_eval.dscreened_rates_dT(k_He4_O16_to_Ne20_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_O16_to_p_F19_derived);
    rate_eval.screened_rates(k_He4_O16_to_p_F19_derived) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_O16_to_p_F19_derived);
        rate_eval.dscreened_rates_dT(k_He4_O16_to_p_F19_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 8.0_rt, 17.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_p_O17_to_F18_reaclib);
    rate_eval.screened_rates(k_p_O17_to_F18_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_O17_to_F18_reaclib);
        rate_eval.dscreened_rates_dT(k_p_O17_to_F18_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_O17_to_He4_N14_reaclib);
    rate_eval.screened_rates(k_p_O17_to_He4_N14_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_O17_to_He4_N14_reaclib);
        rate_eval.dscreened_rates_dT(k_p_O17_to_He4_N14_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 8.0_rt, 17.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_He4_O17_to_Ne21_reaclib);
    rate_eval.screened_rates(k_He4_O17_to_Ne21_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_O17_to_Ne21_reaclib);
        rate_eval.dscreened_rates_dT(k_He4_O17_to_Ne21_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 8.0_rt, 18.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_p_O18_to_F19_reaclib);
    rate_eval.screened_rates(k_p_O18_to_F19_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_O18_to_F19_reaclib);
        rate_eval.dscreened_rates_dT(k_p_O18_to_F19_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_O18_to_He4_N15_reaclib);
    rate_eval.screened_rates(k_p_O18_to_He4_N15_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_O18_to_He4_N15_reaclib);
        rate_eval.dscreened_rates_dT(k_p_O18_to_He4_N15_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 9.0_rt, 17.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_p_F17_to_Ne18_reaclib);
    rate_eval.screened_rates(k_p_F17_to_Ne18_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_F17_to_Ne18_reaclib);
        rate_eval.dscreened_rates_dT(k_p_F17_to_Ne18_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_F17_to_He4_O14_derived);
    rate_eval.screened_rates(k_p_F17_to_He4_O14_derived) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_F17_to_He4_O14_derived);
        rate_eval.dscreened_rates_dT(k_p_F17_to_He4_O14_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 9.0_rt, 18.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_p_F18_to_Ne19_reaclib);
    rate_eval.screened_rates(k_p_F18_to_Ne19_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_F18_to_Ne19_reaclib);
        rate_eval.dscreened_rates_dT(k_p_F18_to_Ne19_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_F18_to_He4_O15_reaclib);
    rate_eval.screened_rates(k_p_F18_to_He4_O15_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_F18_to_He4_O15_reaclib);
        rate_eval.dscreened_rates_dT(k_p_F18_to_He4_O15_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 9.0_rt, 18.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_He4_F18_to_Na22_reaclib);
    rate_eval.screened_rates(k_He4_F18_to_Na22_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_F18_to_Na22_reaclib);
        rate_eval.dscreened_rates_dT(k_He4_F18_to_Na22_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_F18_to_p_Ne21_reaclib);
    rate_eval.screened_rates(k_He4_F18_to_p_Ne21_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_F18_to_p_Ne21_reaclib);
        rate_eval.dscreened_rates_dT(k_He4_F18_to_p_Ne21_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 9.0_rt, 19.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_p_F19_to_Ne20_reaclib);
    rate_eval.screened_rates(k_p_F19_to_Ne20_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_F19_to_Ne20_reaclib);
        rate_eval.dscreened_rates_dT(k_p_F19_to_Ne20_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_F19_to_He4_O16_reaclib);
    rate_eval.screened_rates(k_p_F19_to_He4_O16_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_F19_to_He4_O16_reaclib);
        rate_eval.dscreened_rates_dT(k_p_F19_to_He4_O16_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 9.0_rt, 19.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_He4_F19_to_Na23_reaclib);
    rate_eval.screened_rates(k_He4_F19_to_Na23_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_F19_to_Na23_reaclib);
        rate_eval.dscreened_rates_dT(k_He4_F19_to_Na23_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 10.0_rt, 18.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_He4_Ne18_to_Mg22_reaclib);
    rate_eval.screened_rates(k_He4_Ne18_to_Mg22_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ne18_to_Mg22_reaclib);
        rate_eval.dscreened_rates_dT(k_He4_Ne18_to_Mg22_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 10.0_rt, 20.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_He4_Ne20_to_Mg24_reaclib);
    rate_eval.screened_rates(k_He4_Ne20_to_Mg24_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ne20_to_Mg24_reaclib);
        rate_eval.dscreened_rates_dT(k_He4_Ne20_to_Mg24_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ne20_to_p_Na23_derived);
    rate_eval.screened_rates(k_He4_Ne20_to_p_Na23_derived) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ne20_to_p_Na23_derived);
        rate_eval.dscreened_rates_dT(k_He4_Ne20_to_p_Na23_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 10.0_rt, 21.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_p_Ne21_to_Na22_reaclib);
    rate_eval.screened_rates(k_p_Ne21_to_Na22_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ne21_to_Na22_reaclib);
        rate_eval.dscreened_rates_dT(k_p_Ne21_to_Na22_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ne21_to_He4_F18_derived);
    rate_eval.screened_rates(k_p_Ne21_to_He4_F18_derived) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ne21_to_He4_F18_derived);
        rate_eval.dscreened_rates_dT(k_p_Ne21_to_He4_F18_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 11.0_rt, 23.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_p_Na23_to_Mg24_reaclib);
    rate_eval.screened_rates(k_p_Na23_to_Mg24_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Na23_to_Mg24_reaclib);
        rate_eval.dscreened_rates_dT(k_p_Na23_to_Mg24_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Na23_to_He4_Ne20_reaclib);
    rate_eval.screened_rates(k_p_Na23_to_He4_Ne20_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Na23_to_He4_Ne20_reaclib);
        rate_eval.dscreened_rates_dT(k_p_Na23_to_He4_Ne20_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 12.0_rt, 24.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_He4_Mg24_to_Si28_reaclib);
    rate_eval.screened_rates(k_He4_Mg24_to_Si28_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mg24_to_Si28_reaclib);
        rate_eval.dscreened_rates_dT(k_He4_Mg24_to_Si28_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Mg24_to_p_Al27_derived);
    rate_eval.screened_rates(k_He4_Mg24_to_p_Al27_derived) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mg24_to_p_Al27_derived);
        rate_eval.dscreened_rates_dT(k_He4_Mg24_to_p_Al27_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 13.0_rt, 27.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_p_Al27_to_Si28_reaclib);
    rate_eval.screened_rates(k_p_Al27_to_Si28_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Al27_to_Si28_reaclib);
        rate_eval.dscreened_rates_dT(k_p_Al27_to_Si28_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Al27_to_He4_Mg24_reaclib);
    rate_eval.screened_rates(k_p_Al27_to_He4_Mg24_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Al27_to_He4_Mg24_reaclib);
        rate_eval.dscreened_rates_dT(k_p_Al27_to_He4_Mg24_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 14.0_rt, 28.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_He4_Si28_to_S32_reaclib);
    rate_eval.screened_rates(k_He4_Si28_to_S32_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Si28_to_S32_reaclib);
        rate_eval.dscreened_rates_dT(k_He4_Si28_to_S32_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Si28_to_p_P31_derived);
    rate_eval.screened_rates(k_He4_Si28_to_p_P31_derived) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Si28_to_p_P31_derived);
        rate_eval.dscreened_rates_dT(k_He4_Si28_to_p_P31_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 15.0_rt, 31.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_p_P31_to_S32_reaclib);
    rate_eval.screened_rates(k_p_P31_to_S32_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_P31_to_S32_reaclib);
        rate_eval.dscreened_rates_dT(k_p_P31_to_S32_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_P31_to_He4_Si28_reaclib);
    rate_eval.screened_rates(k_p_P31_to_He4_Si28_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_P31_to_He4_Si28_reaclib);
        rate_eval.dscreened_rates_dT(k_p_P31_to_He4_Si28_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(6.0_rt, 12.0_rt, 6.0_rt, 12.0_rt);
        static_assert(scn_fac.z1 == 6.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_C12_C12_to_p_Na23_reaclib);
    rate_eval.screened_rates(k_C12_C12_to_p_Na23_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_C12_to_p_Na23_reaclib);
        rate_eval.dscreened_rates_dT(k_C12_C12_to_p_Na23_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_C12_C12_to_He4_Ne20_reaclib);
    rate_eval.screened_rates(k_C12_C12_to_He4_Ne20_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_C12_to_He4_Ne20_reaclib);
        rate_eval.dscreened_rates_dT(k_C12_C12_to_He4_Ne20_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_C12_C12_to_Mg24_modified);
    rate_eval.screened_rates(k_C12_C12_to_Mg24_modified) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_C12_to_Mg24_modified);
        rate_eval.dscreened_rates_dT(k_C12_C12_to_Mg24_modified) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 7.0_rt, 13.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_He4_N13_to_p_O16_reaclib);
    rate_eval.screened_rates(k_He4_N13_to_p_O16_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_N13_to_p_O16_reaclib);
        rate_eval.dscreened_rates_dT(k_He4_N13_to_p_O16_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(6.0_rt, 12.0_rt, 8.0_rt, 16.0_rt);
        static_assert(scn_fac.z1 == 6.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_C12_O16_to_p_Al27_reaclib);
    rate_eval.screened_rates(k_C12_O16_to_p_Al27_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_O16_to_p_Al27_reaclib);
        rate_eval.dscreened_rates_dT(k_C12_O16_to_p_Al27_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_C12_O16_to_He4_Mg24_reaclib);
    rate_eval.screened_rates(k_C12_O16_to_He4_Mg24_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_O16_to_He4_Mg24_reaclib);
        rate_eval.dscreened_rates_dT(k_C12_O16_to_He4_Mg24_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_C12_O16_to_Si28_modified);
    rate_eval.screened_rates(k_C12_O16_to_Si28_modified) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_O16_to_Si28_modified);
        rate_eval.dscreened_rates_dT(k_C12_O16_to_Si28_modified) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(8.0_rt, 16.0_rt, 8.0_rt, 16.0_rt);
        static_assert(scn_fac.z1 == 8.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_O16_O16_to_p_P31_reaclib);
    rate_eval.screened_rates(k_O16_O16_to_p_P31_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_O16_O16_to_p_P31_reaclib);
        rate_eval.dscreened_rates_dT(k_O16_O16_to_p_P31_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_O16_O16_to_He4_Si28_reaclib);
    rate_eval.screened_rates(k_O16_O16_to_He4_Si28_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_O16_O16_to_He4_Si28_reaclib);
        rate_eval.dscreened_rates_dT(k_O16_O16_to_He4_Si28_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_O16_O16_to_S32_modified);
    rate_eval.screened_rates(k_O16_O16_to_S32_modified) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_O16_O16_to_S32_modified);
        rate_eval.dscreened_rates_dT(k_O16_O16_to_S32_modified) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 10.0_rt, 19.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_He4_Ne19_to_p_Na22_reaclib);
    rate_eval.screened_rates(k_He4_Ne19_to_p_Na22_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ne19_to_p_Na22_reaclib);
        rate_eval.dscreened_rates_dT(k_He4_Ne19_to_p_Na22_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 10.0_rt, 20.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_p_Ne20_to_He4_F17_reaclib);
    rate_eval.screened_rates(k_p_Ne20_to_He4_F17_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ne20_to_He4_F17_reaclib);
        rate_eval.dscreened_rates_dT(k_p_Ne20_to_He4_F17_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 2.0_rt, 4.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    {
        constexpr auto scn_fac2 = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 4.0_rt, 8.0_rt);
        static_assert(scn_fac2.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac2, scor2, dscor2_dt);
    }

    ratraw = rate_eval.screened_rates(k_He4_He4_He4_to_C12_reaclib);
    rate_eval.screened_rates(k_He4_He4_He4_to_C12_reaclib) *= scor * scor2;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_He4_He4_to_C12_reaclib);
        rate_eval.dscreened_rates_dT(k_He4_He4_He4_to_C12_reaclib) = ratraw * (scor * dscor2_dt + dscor_dt * scor2) + dratraw_dT * scor * scor2;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 9.0_rt, 17.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_He4_F17_to_p_Ne20_derived);
    rate_eval.screened_rates(k_He4_F17_to_p_Ne20_derived) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_F17_to_p_Ne20_derived);
        rate_eval.dscreened_rates_dT(k_He4_F17_to_p_Ne20_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 11.0_rt, 22.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_p_Na22_to_He4_Ne19_derived);
    rate_eval.screened_rates(k_p_Na22_to_He4_Ne19_derived) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Na22_to_He4_Ne19_derived);
        rate_eval.dscreened_rates_dT(k_p_Na22_to_He4_Ne19_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 16.0_rt, 32.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_He4_S32_to_Ar36_removed);
    rate_eval.screened_rates(k_He4_S32_to_Ar36_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_S32_to_Ar36_removed);
        rate_eval.dscreened_rates_dT(k_He4_S32_to_Ar36_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_S32_to_p_Cl35_removed);
    rate_eval.screened_rates(k_He4_S32_to_p_Cl35_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_S32_to_p_Cl35_removed);
        rate_eval.dscreened_rates_dT(k_He4_S32_to_p_Cl35_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 17.0_rt, 35.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_p_Cl35_to_Ar36_removed);
    rate_eval.screened_rates(k_p_Cl35_to_Ar36_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cl35_to_Ar36_removed);
        rate_eval.dscreened_rates_dT(k_p_Cl35_to_Ar36_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cl35_to_He4_S32_removed);
    rate_eval.screened_rates(k_p_Cl35_to_He4_S32_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cl35_to_He4_S32_removed);
        rate_eval.dscreened_rates_dT(k_p_Cl35_to_He4_S32_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 18.0_rt, 36.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_He4_Ar36_to_Ca40_removed);
    rate_eval.screened_rates(k_He4_Ar36_to_Ca40_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ar36_to_Ca40_removed);
        rate_eval.dscreened_rates_dT(k_He4_Ar36_to_Ca40_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ar36_to_p_K39_removed);
    rate_eval.screened_rates(k_He4_Ar36_to_p_K39_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ar36_to_p_K39_removed);
        rate_eval.dscreened_rates_dT(k_He4_Ar36_to_p_K39_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 19.0_rt, 39.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_p_K39_to_Ca40_removed);
    rate_eval.screened_rates(k_p_K39_to_Ca40_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_K39_to_Ca40_removed);
        rate_eval.dscreened_rates_dT(k_p_K39_to_Ca40_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_K39_to_He4_Ar36_removed);
    rate_eval.screened_rates(k_p_K39_to_He4_Ar36_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_K39_to_He4_Ar36_removed);
        rate_eval.dscreened_rates_dT(k_p_K39_to_He4_Ar36_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 20.0_rt, 40.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_He4_Ca40_to_Ti44_removed);
    rate_eval.screened_rates(k_He4_Ca40_to_Ti44_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca40_to_Ti44_removed);
        rate_eval.dscreened_rates_dT(k_He4_Ca40_to_Ti44_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ca40_to_p_Sc43_removed);
    rate_eval.screened_rates(k_He4_Ca40_to_p_Sc43_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca40_to_p_Sc43_removed);
        rate_eval.dscreened_rates_dT(k_He4_Ca40_to_p_Sc43_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 21.0_rt, 43.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_p_Sc43_to_Ti44_removed);
    rate_eval.screened_rates(k_p_Sc43_to_Ti44_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Sc43_to_Ti44_removed);
        rate_eval.dscreened_rates_dT(k_p_Sc43_to_Ti44_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Sc43_to_He4_Ca40_removed);
    rate_eval.screened_rates(k_p_Sc43_to_He4_Ca40_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Sc43_to_He4_Ca40_removed);
        rate_eval.dscreened_rates_dT(k_p_Sc43_to_He4_Ca40_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 22.0_rt, 44.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_He4_Ti44_to_Cr48_removed);
    rate_eval.screened_rates(k_He4_Ti44_to_Cr48_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ti44_to_Cr48_removed);
        rate_eval.dscreened_rates_dT(k_He4_Ti44_to_Cr48_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ti44_to_p_V47_removed);
    rate_eval.screened_rates(k_He4_Ti44_to_p_V47_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ti44_to_p_V47_removed);
        rate_eval.dscreened_rates_dT(k_He4_Ti44_to_p_V47_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 23.0_rt, 47.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_p_V47_to_Cr48_removed);
    rate_eval.screened_rates(k_p_V47_to_Cr48_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_V47_to_Cr48_removed);
        rate_eval.dscreened_rates_dT(k_p_V47_to_Cr48_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_V47_to_He4_Ti44_removed);
    rate_eval.screened_rates(k_p_V47_to_He4_Ti44_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_V47_to_He4_Ti44_removed);
        rate_eval.dscreened_rates_dT(k_p_V47_to_He4_Ti44_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 24.0_rt, 48.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_He4_Cr48_to_Fe52_removed);
    rate_eval.screened_rates(k_He4_Cr48_to_Fe52_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cr48_to_Fe52_removed);
        rate_eval.dscreened_rates_dT(k_He4_Cr48_to_Fe52_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cr48_to_p_Mn51_removed);
    rate_eval.screened_rates(k_He4_Cr48_to_p_Mn51_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cr48_to_p_Mn51_removed);
        rate_eval.dscreened_rates_dT(k_He4_Cr48_to_p_Mn51_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 25.0_rt, 51.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_p_Mn51_to_Fe52_removed);
    rate_eval.screened_rates(k_p_Mn51_to_Fe52_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mn51_to_Fe52_removed);
        rate_eval.dscreened_rates_dT(k_p_Mn51_to_Fe52_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Mn51_to_He4_Cr48_removed);
    rate_eval.screened_rates(k_p_Mn51_to_He4_Cr48_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mn51_to_He4_Cr48_removed);
        rate_eval.dscreened_rates_dT(k_p_Mn51_to_He4_Cr48_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 26.0_rt, 52.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_He4_Fe52_to_Ni56_removed);
    rate_eval.screened_rates(k_He4_Fe52_to_Ni56_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe52_to_Ni56_removed);
        rate_eval.dscreened_rates_dT(k_He4_Fe52_to_Ni56_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Fe52_to_p_Co55_removed);
    rate_eval.screened_rates(k_He4_Fe52_to_p_Co55_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe52_to_p_Co55_removed);
        rate_eval.dscreened_rates_dT(k_He4_Fe52_to_p_Co55_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 27.0_rt, 55.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_p_Co55_to_Ni56_removed);
    rate_eval.screened_rates(k_p_Co55_to_Ni56_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Co55_to_Ni56_removed);
        rate_eval.dscreened_rates_dT(k_p_Co55_to_Ni56_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Co55_to_He4_Fe52_removed);
    rate_eval.screened_rates(k_p_Co55_to_He4_Fe52_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Co55_to_He4_Fe52_removed);
        rate_eval.dscreened_rates_dT(k_p_Co55_to_He4_Fe52_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }

#endif

    // Fill approximate rates

    fill_approx_rates<do_T_derivatives, T>(tfactors, state.rho, Y, rate_eval);

    // Calculate tabular rates

    [[maybe_unused]] amrex::Real rate, drate_dt, edot_nu, edot_gamma;

    rate_eval.enuc_weak = 0.0_rt;

    amrex::Real log_temp = std::log10(state.T);
    amrex::Real log_rhoy = std::log10(rhoy);

    tabular_evaluate(j_F17_O17_meta, j_F17_O17_rhoy, j_F17_O17_temp, j_F17_O17_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_F17_to_O17_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_F17_to_O17_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(F17) * (edot_nu + edot_gamma);

    tabular_evaluate(j_O17_F17_meta, j_O17_F17_rhoy, j_O17_F17_temp, j_O17_F17_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_O17_to_F17_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_O17_to_F17_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(O17) * (edot_nu + edot_gamma);

    tabular_evaluate(j_F18_Ne18_meta, j_F18_Ne18_rhoy, j_F18_Ne18_temp, j_F18_Ne18_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_F18_to_Ne18_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_F18_to_Ne18_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(F18) * (edot_nu + edot_gamma);

    tabular_evaluate(j_F18_O18_meta, j_F18_O18_rhoy, j_F18_O18_temp, j_F18_O18_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_F18_to_O18_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_F18_to_O18_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(F18) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Ne18_F18_meta, j_Ne18_F18_rhoy, j_Ne18_F18_temp, j_Ne18_F18_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Ne18_to_F18_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ne18_to_F18_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(Ne18) * (edot_nu + edot_gamma);

    tabular_evaluate(j_O18_F18_meta, j_O18_F18_rhoy, j_O18_F18_temp, j_O18_F18_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_O18_to_F18_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_O18_to_F18_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(O18) * (edot_nu + edot_gamma);

    tabular_evaluate(j_F19_Ne19_meta, j_F19_Ne19_rhoy, j_F19_Ne19_temp, j_F19_Ne19_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_F19_to_Ne19_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_F19_to_Ne19_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(F19) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Ne19_F19_meta, j_Ne19_F19_rhoy, j_Ne19_F19_temp, j_Ne19_F19_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Ne19_to_F19_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ne19_to_F19_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(Ne19) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Mg22_Na22_meta, j_Mg22_Na22_rhoy, j_Mg22_Na22_temp, j_Mg22_Na22_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Mg22_to_Na22_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Mg22_to_Na22_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(Mg22) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Na22_Mg22_meta, j_Na22_Mg22_rhoy, j_Na22_Mg22_temp, j_Na22_Mg22_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Na22_to_Mg22_weaktab) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Na22_to_Mg22_weaktab) = drate_dt;
    }
    rate_eval.enuc_weak += C::n_A * Y(Na22) * (edot_nu + edot_gamma);


}


// backwards-compatible wrapper
template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
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
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void get_ydot_weak(const burn_t& state,
             amrex::Array1D<amrex::Real, 1, neqs>& ydot_nuc,
             amrex::Real& enuc_weak,
             [[maybe_unused]] const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {
    ///
    /// Calculate Ydots contribute only from weak reactions.
    /// This is used to calculate dyedt and energy generation from
    /// weak reactions for self-consistent NSE
    ///


    // initialize ydot_nuc to 0

    for (int i = 1; i <= neqs; ++i) {
        ydot_nuc(i) = 0.0_rt;
    }

    rate_t rate_eval;

    [[maybe_unused]] amrex::Real rate, drate_dt, edot_nu, edot_gamma;
    [[maybe_unused]] amrex::Real rhoy = state.rho * state.y_e;

    rate_eval.enuc_weak = 0.0_rt;

    // Calculate tabular rates and get ydot_weak

    amrex::Real log_temp = std::log10(state.T);
    amrex::Real log_rhoy = std::log10(rhoy);

    tabular_evaluate(j_F17_O17_meta, j_F17_O17_rhoy, j_F17_O17_temp, j_F17_O17_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_F17_to_O17_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(F17) * (edot_nu + edot_gamma);

    tabular_evaluate(j_O17_F17_meta, j_O17_F17_rhoy, j_O17_F17_temp, j_O17_F17_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_O17_to_F17_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(O17) * (edot_nu + edot_gamma);

    tabular_evaluate(j_F18_Ne18_meta, j_F18_Ne18_rhoy, j_F18_Ne18_temp, j_F18_Ne18_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_F18_to_Ne18_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(F18) * (edot_nu + edot_gamma);

    tabular_evaluate(j_F18_O18_meta, j_F18_O18_rhoy, j_F18_O18_temp, j_F18_O18_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_F18_to_O18_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(F18) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Ne18_F18_meta, j_Ne18_F18_rhoy, j_Ne18_F18_temp, j_Ne18_F18_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Ne18_to_F18_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(Ne18) * (edot_nu + edot_gamma);

    tabular_evaluate(j_O18_F18_meta, j_O18_F18_rhoy, j_O18_F18_temp, j_O18_F18_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_O18_to_F18_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(O18) * (edot_nu + edot_gamma);

    tabular_evaluate(j_F19_Ne19_meta, j_F19_Ne19_rhoy, j_F19_Ne19_temp, j_F19_Ne19_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_F19_to_Ne19_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(F19) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Ne19_F19_meta, j_Ne19_F19_rhoy, j_Ne19_F19_temp, j_Ne19_F19_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Ne19_to_F19_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(Ne19) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Mg22_Na22_meta, j_Mg22_Na22_rhoy, j_Mg22_Na22_temp, j_Mg22_Na22_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Mg22_to_Na22_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(Mg22) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Na22_Mg22_meta, j_Na22_Mg22_rhoy, j_Na22_Mg22_temp, j_Na22_Mg22_data,
                     log_rhoy, log_temp, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Na22_to_Mg22_weaktab) = rate;
    rate_eval.enuc_weak += C::n_A * Y(Na22) * (edot_nu + edot_gamma);

    auto screened_rates = rate_eval.screened_rates;

    ydot_nuc(H1) = 0.0_rt;

    ydot_nuc(He4) = 0.0_rt;

    ydot_nuc(C12) = 0.0_rt;

    ydot_nuc(C13) =
        screened_rates(k_N13_to_C13_reaclib)*Y(N13);

    ydot_nuc(N13) =
        -screened_rates(k_N13_to_C13_reaclib)*Y(N13);

    ydot_nuc(N14) =
        screened_rates(k_O14_to_N14_reaclib)*Y(O14);

    ydot_nuc(N15) =
        screened_rates(k_O15_to_N15_reaclib)*Y(O15);

    ydot_nuc(O14) =
        -screened_rates(k_O14_to_N14_reaclib)*Y(O14);

    ydot_nuc(O15) =
        -screened_rates(k_O15_to_N15_reaclib)*Y(O15);

    ydot_nuc(O16) = 0.0_rt;

    ydot_nuc(O17) =
        (screened_rates(k_F17_to_O17_weaktab)*Y(F17) + -screened_rates(k_O17_to_F17_weaktab)*Y(O17));

    ydot_nuc(O18) =
        (screened_rates(k_F18_to_O18_weaktab)*Y(F18) + -screened_rates(k_O18_to_F18_weaktab)*Y(O18));

    ydot_nuc(F17) =
        (-screened_rates(k_F17_to_O17_weaktab)*Y(F17) + screened_rates(k_O17_to_F17_weaktab)*Y(O17));

    ydot_nuc(F18) =
        (-screened_rates(k_F18_to_O18_weaktab)*Y(F18) + screened_rates(k_O18_to_F18_weaktab)*Y(O18)) +
        (screened_rates(k_Ne18_to_F18_weaktab)*Y(Ne18) + -screened_rates(k_F18_to_Ne18_weaktab)*Y(F18));

    ydot_nuc(F19) =
        (screened_rates(k_Ne19_to_F19_weaktab)*Y(Ne19) + -screened_rates(k_F19_to_Ne19_weaktab)*Y(F19));

    ydot_nuc(Ne18) =
        (-screened_rates(k_Ne18_to_F18_weaktab)*Y(Ne18) + screened_rates(k_F18_to_Ne18_weaktab)*Y(F18));

    ydot_nuc(Ne19) =
        (-screened_rates(k_Ne19_to_F19_weaktab)*Y(Ne19) + screened_rates(k_F19_to_Ne19_weaktab)*Y(F19));

    ydot_nuc(Ne20) = 0.0_rt;

    ydot_nuc(Ne21) = 0.0_rt;

    ydot_nuc(Na22) =
        (screened_rates(k_Mg22_to_Na22_weaktab)*Y(Mg22) + -screened_rates(k_Na22_to_Mg22_weaktab)*Y(Na22));

    ydot_nuc(Na23) = 0.0_rt;

    ydot_nuc(Mg22) =
        (-screened_rates(k_Mg22_to_Na22_weaktab)*Y(Mg22) + screened_rates(k_Na22_to_Mg22_weaktab)*Y(Na22));

    ydot_nuc(Mg24) = 0.0_rt;

    ydot_nuc(Al27) = 0.0_rt;

    ydot_nuc(Si28) = 0.0_rt;

    ydot_nuc(P31) = 0.0_rt;

    ydot_nuc(S32) = 0.0_rt;

    ydot_nuc(Ar36) = 0.0_rt;

    ydot_nuc(Ca40) = 0.0_rt;

    ydot_nuc(Ti44) = 0.0_rt;

    ydot_nuc(Cr48) = 0.0_rt;

    ydot_nuc(Fe52) = 0.0_rt;

    ydot_nuc(Ni56) = 0.0_rt;

    enuc_weak = rate_eval.enuc_weak;
}
#endif


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rhs_nuc(const burn_t& state,
             amrex::Array1D<amrex::Real, 1, neqs>& ydot_nuc,
             const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
             const amrex::Array1D<amrex::Real, 1, NumRates>& screened_rates) {

    using namespace Rates;

    ydot_nuc(H1) =
        (-screened_rates(k_p_C12_to_N13_reaclib)*Y(C12)*Y(H1)*state.rho + screened_rates(k_N13_to_p_C12_derived)*Y(N13)) +
        (-screened_rates(k_p_C13_to_N14_reaclib)*Y(C13)*Y(H1)*state.rho + screened_rates(k_N14_to_p_C13_derived)*Y(N14)) +
        (-screened_rates(k_p_N13_to_O14_reaclib)*Y(N13)*Y(H1)*state.rho + screened_rates(k_O14_to_p_N13_derived)*Y(O14)) +
        (-screened_rates(k_p_N14_to_O15_reaclib)*Y(N14)*Y(H1)*state.rho + screened_rates(k_O15_to_p_N14_derived)*Y(O15)) +
        (-screened_rates(k_p_N15_to_O16_reaclib)*Y(N15)*Y(H1)*state.rho + screened_rates(k_O16_to_p_N15_derived)*Y(O16)) +
        (-screened_rates(k_p_O16_to_F17_reaclib)*Y(O16)*Y(H1)*state.rho + screened_rates(k_F17_to_p_O16_derived)*Y(F17)) +
        (-screened_rates(k_p_O17_to_F18_reaclib)*Y(O17)*Y(H1)*state.rho + screened_rates(k_F18_to_p_O17_derived)*Y(F18)) +
        (-screened_rates(k_p_O18_to_F19_reaclib)*Y(O18)*Y(H1)*state.rho + screened_rates(k_F19_to_p_O18_derived)*Y(F19)) +
        (-screened_rates(k_p_F17_to_Ne18_reaclib)*Y(F17)*Y(H1)*state.rho + screened_rates(k_Ne18_to_p_F17_derived)*Y(Ne18)) +
        (-screened_rates(k_p_F18_to_Ne19_reaclib)*Y(F18)*Y(H1)*state.rho + screened_rates(k_Ne19_to_p_F18_derived)*Y(Ne19)) +
        (-screened_rates(k_p_F19_to_Ne20_reaclib)*Y(F19)*Y(H1)*state.rho + screened_rates(k_Ne20_to_p_F19_derived)*Y(Ne20)) +
        (-screened_rates(k_p_Ne21_to_Na22_reaclib)*Y(Ne21)*Y(H1)*state.rho + screened_rates(k_Na22_to_p_Ne21_derived)*Y(Na22)) +
        (-screened_rates(k_p_Na23_to_Mg24_reaclib)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_Mg24_to_p_Na23_derived)*Y(Mg24)) +
        (-screened_rates(k_p_Al27_to_Si28_reaclib)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_Si28_to_p_Al27_derived)*Y(Si28)) +
        (-screened_rates(k_p_P31_to_S32_reaclib)*Y(P31)*Y(H1)*state.rho + screened_rates(k_S32_to_p_P31_derived)*Y(S32)) +
        0.5*screened_rates(k_C12_C12_to_p_Na23_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho +
        (screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(He4)*Y(N13)*state.rho + -screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*Y(H1)*state.rho) +
        (-screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(N15)*Y(H1)*state.rho + screened_rates(k_He4_C12_to_p_N15_derived)*Y(C12)*Y(He4)*state.rho) +
        (screened_rates(k_He4_O14_to_p_F17_reaclib)*Y(He4)*Y(O14)*state.rho + -screened_rates(k_p_F17_to_He4_O14_derived)*Y(F17)*Y(H1)*state.rho) +
        screened_rates(k_C12_O16_to_p_Al27_reaclib)*Y(C12)*Y(O16)*state.rho +
        0.5*screened_rates(k_O16_O16_to_p_P31_reaclib)*amrex::Math::powi<2>(Y(O16))*state.rho +
        (-screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(O17)*Y(H1)*state.rho + screened_rates(k_He4_N14_to_p_O17_derived)*Y(He4)*Y(N14)*state.rho) +
        (-screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(O18)*Y(H1)*state.rho + screened_rates(k_He4_N15_to_p_O18_derived)*Y(He4)*Y(N15)*state.rho) +
        (-screened_rates(k_p_F18_to_He4_O15_reaclib)*Y(F18)*Y(H1)*state.rho + screened_rates(k_He4_O15_to_p_F18_derived)*Y(He4)*Y(O15)*state.rho) +
        (screened_rates(k_He4_F18_to_p_Ne21_reaclib)*Y(F18)*Y(He4)*state.rho + -screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(Ne21)*Y(H1)*state.rho) +
        (-screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(F19)*Y(H1)*state.rho + screened_rates(k_He4_O16_to_p_F19_derived)*Y(He4)*Y(O16)*state.rho) +
        (screened_rates(k_He4_Ne19_to_p_Na22_reaclib)*Y(He4)*Y(Ne19)*state.rho + -screened_rates(k_p_Na22_to_He4_Ne19_derived)*Y(Na22)*Y(H1)*state.rho) +
        (-screened_rates(k_p_Na23_to_He4_Ne20_reaclib)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(He4)*Y(Ne20)*state.rho) +
        (-screened_rates(k_p_Al27_to_He4_Mg24_reaclib)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(He4)*Y(Mg24)*state.rho) +
        (-screened_rates(k_p_P31_to_He4_Si28_reaclib)*Y(P31)*Y(H1)*state.rho + screened_rates(k_He4_Si28_to_p_P31_derived)*Y(He4)*Y(Si28)*state.rho) +
        (screened_rates(k_He4_F17_to_p_Ne20_derived)*Y(F17)*Y(He4)*state.rho + -screened_rates(k_p_Ne20_to_He4_F17_reaclib)*Y(Ne20)*Y(H1)*state.rho);

    ydot_nuc(He4) =
        (-screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*Y(He4)*state.rho + screened_rates(k_O16_to_He4_C12_derived)*Y(O16)) +
        (-screened_rates(k_He4_N14_to_F18_reaclib)*Y(He4)*Y(N14)*state.rho + screened_rates(k_F18_to_He4_N14_derived)*Y(F18)) +
        (-screened_rates(k_He4_N15_to_F19_reaclib)*Y(He4)*Y(N15)*state.rho + screened_rates(k_F19_to_He4_N15_derived)*Y(F19)) +
        (-screened_rates(k_He4_O14_to_Ne18_reaclib)*Y(He4)*Y(O14)*state.rho + screened_rates(k_Ne18_to_He4_O14_derived)*Y(Ne18)) +
        (-screened_rates(k_He4_O15_to_Ne19_reaclib)*Y(He4)*Y(O15)*state.rho + screened_rates(k_Ne19_to_He4_O15_derived)*Y(Ne19)) +
        (-screened_rates(k_He4_O16_to_Ne20_reaclib)*Y(He4)*Y(O16)*state.rho + screened_rates(k_Ne20_to_He4_O16_derived)*Y(Ne20)) +
        (-screened_rates(k_He4_O17_to_Ne21_reaclib)*Y(He4)*Y(O17)*state.rho + screened_rates(k_Ne21_to_He4_O17_derived)*Y(Ne21)) +
        (-screened_rates(k_He4_F18_to_Na22_reaclib)*Y(F18)*Y(He4)*state.rho + screened_rates(k_Na22_to_He4_F18_derived)*Y(Na22)) +
        (-screened_rates(k_He4_F19_to_Na23_reaclib)*Y(F19)*Y(He4)*state.rho + screened_rates(k_Na23_to_He4_F19_derived)*Y(Na23)) +
        (-screened_rates(k_He4_Ne18_to_Mg22_reaclib)*Y(He4)*Y(Ne18)*state.rho + screened_rates(k_Mg22_to_He4_Ne18_derived)*Y(Mg22)) +
        (-screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(He4)*Y(Ne20)*state.rho + screened_rates(k_Mg24_to_He4_Ne20_derived)*Y(Mg24)) +
        (-screened_rates(k_He4_Mg24_to_Si28_reaclib)*Y(He4)*Y(Mg24)*state.rho + screened_rates(k_Si28_to_He4_Mg24_derived)*Y(Si28)) +
        (-screened_rates(k_He4_Si28_to_S32_reaclib)*Y(He4)*Y(Si28)*state.rho + screened_rates(k_S32_to_He4_Si28_derived)*Y(S32)) +
        0.5*screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho +
        (-screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(He4)*Y(N13)*state.rho + screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*Y(H1)*state.rho) +
        (screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(N15)*Y(H1)*state.rho + -screened_rates(k_He4_C12_to_p_N15_derived)*Y(C12)*Y(He4)*state.rho) +
        (-screened_rates(k_He4_O14_to_p_F17_reaclib)*Y(He4)*Y(O14)*state.rho + screened_rates(k_p_F17_to_He4_O14_derived)*Y(F17)*Y(H1)*state.rho) +
        screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(C12)*Y(O16)*state.rho +
        0.5*screened_rates(k_O16_O16_to_He4_Si28_reaclib)*amrex::Math::powi<2>(Y(O16))*state.rho +
        (screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(O17)*Y(H1)*state.rho + -screened_rates(k_He4_N14_to_p_O17_derived)*Y(He4)*Y(N14)*state.rho) +
        (screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(O18)*Y(H1)*state.rho + -screened_rates(k_He4_N15_to_p_O18_derived)*Y(He4)*Y(N15)*state.rho) +
        (screened_rates(k_p_F18_to_He4_O15_reaclib)*Y(F18)*Y(H1)*state.rho + -screened_rates(k_He4_O15_to_p_F18_derived)*Y(He4)*Y(O15)*state.rho) +
        (-screened_rates(k_He4_F18_to_p_Ne21_reaclib)*Y(F18)*Y(He4)*state.rho + screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(Ne21)*Y(H1)*state.rho) +
        (screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(F19)*Y(H1)*state.rho + -screened_rates(k_He4_O16_to_p_F19_derived)*Y(He4)*Y(O16)*state.rho) +
        (-screened_rates(k_He4_Ne19_to_p_Na22_reaclib)*Y(He4)*Y(Ne19)*state.rho + screened_rates(k_p_Na22_to_He4_Ne19_derived)*Y(Na22)*Y(H1)*state.rho) +
        (screened_rates(k_p_Na23_to_He4_Ne20_reaclib)*Y(Na23)*Y(H1)*state.rho + -screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(He4)*Y(Ne20)*state.rho) +
        (screened_rates(k_p_Al27_to_He4_Mg24_reaclib)*Y(Al27)*Y(H1)*state.rho + -screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(He4)*Y(Mg24)*state.rho) +
        (screened_rates(k_p_P31_to_He4_Si28_reaclib)*Y(P31)*Y(H1)*state.rho + -screened_rates(k_He4_Si28_to_p_P31_derived)*Y(He4)*Y(Si28)*state.rho) +
        (-0.5*screened_rates(k_He4_He4_He4_to_C12_reaclib)*amrex::Math::powi<3>(Y(He4))*amrex::Math::powi<2>(state.rho) + 3.0*screened_rates(k_C12_to_He4_He4_He4_derived)*Y(C12)) +
        (-screened_rates(k_He4_F17_to_p_Ne20_derived)*Y(F17)*Y(He4)*state.rho + screened_rates(k_p_Ne20_to_He4_F17_reaclib)*Y(Ne20)*Y(H1)*state.rho) +
        (-screened_rates(k_S32_He4_to_Ar36_approx)*Y(He4)*Y(S32)*state.rho + screened_rates(k_Ar36_to_S32_He4_approx)*Y(Ar36)) +
        (-screened_rates(k_Ar36_He4_to_Ca40_approx)*Y(Ar36)*Y(He4)*state.rho + screened_rates(k_Ca40_to_Ar36_He4_approx)*Y(Ca40)) +
        (-screened_rates(k_Ca40_He4_to_Ti44_approx)*Y(Ca40)*Y(He4)*state.rho + screened_rates(k_Ti44_to_Ca40_He4_approx)*Y(Ti44)) +
        (-screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(He4)*Y(Ti44)*state.rho + screened_rates(k_Cr48_to_Ti44_He4_approx)*Y(Cr48)) +
        (-screened_rates(k_Cr48_He4_to_Fe52_approx)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_Fe52_to_Cr48_He4_approx)*Y(Fe52)) +
        (-screened_rates(k_Fe52_He4_to_Ni56_approx)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_Ni56_to_Fe52_He4_approx)*Y(Ni56));

    ydot_nuc(C12) =
        (-screened_rates(k_p_C12_to_N13_reaclib)*Y(C12)*Y(H1)*state.rho + screened_rates(k_N13_to_p_C12_derived)*Y(N13)) +
        (-screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*Y(He4)*state.rho + screened_rates(k_O16_to_He4_C12_derived)*Y(O16)) +
        -screened_rates(k_C12_C12_to_p_Na23_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho +
        -screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho +
        (screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(N15)*Y(H1)*state.rho + -screened_rates(k_He4_C12_to_p_N15_derived)*Y(C12)*Y(He4)*state.rho) +
        -screened_rates(k_C12_O16_to_p_Al27_reaclib)*Y(C12)*Y(O16)*state.rho +
        -screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(C12)*Y(O16)*state.rho +
        (0.16666666666666667*screened_rates(k_He4_He4_He4_to_C12_reaclib)*amrex::Math::powi<3>(Y(He4))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_C12_to_He4_He4_He4_derived)*Y(C12)) +
        -screened_rates(k_C12_C12_to_Mg24_modified)*amrex::Math::powi<2>(Y(C12))*state.rho +
        -screened_rates(k_C12_O16_to_Si28_modified)*Y(C12)*Y(O16)*state.rho;

    ydot_nuc(C13) =
        screened_rates(k_N13_to_C13_reaclib)*Y(N13) +
        (-screened_rates(k_p_C13_to_N14_reaclib)*Y(C13)*Y(H1)*state.rho + screened_rates(k_N14_to_p_C13_derived)*Y(N14));

    ydot_nuc(N13) =
        -screened_rates(k_N13_to_C13_reaclib)*Y(N13) +
        (screened_rates(k_p_C12_to_N13_reaclib)*Y(C12)*Y(H1)*state.rho + -screened_rates(k_N13_to_p_C12_derived)*Y(N13)) +
        (-screened_rates(k_p_N13_to_O14_reaclib)*Y(N13)*Y(H1)*state.rho + screened_rates(k_O14_to_p_N13_derived)*Y(O14)) +
        (-screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(He4)*Y(N13)*state.rho + screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*Y(H1)*state.rho);

    ydot_nuc(N14) =
        screened_rates(k_O14_to_N14_reaclib)*Y(O14) +
        (screened_rates(k_p_C13_to_N14_reaclib)*Y(C13)*Y(H1)*state.rho + -screened_rates(k_N14_to_p_C13_derived)*Y(N14)) +
        (-screened_rates(k_p_N14_to_O15_reaclib)*Y(N14)*Y(H1)*state.rho + screened_rates(k_O15_to_p_N14_derived)*Y(O15)) +
        (-screened_rates(k_He4_N14_to_F18_reaclib)*Y(He4)*Y(N14)*state.rho + screened_rates(k_F18_to_He4_N14_derived)*Y(F18)) +
        (screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(O17)*Y(H1)*state.rho + -screened_rates(k_He4_N14_to_p_O17_derived)*Y(He4)*Y(N14)*state.rho);

    ydot_nuc(N15) =
        screened_rates(k_O15_to_N15_reaclib)*Y(O15) +
        (-screened_rates(k_p_N15_to_O16_reaclib)*Y(N15)*Y(H1)*state.rho + screened_rates(k_O16_to_p_N15_derived)*Y(O16)) +
        (-screened_rates(k_He4_N15_to_F19_reaclib)*Y(He4)*Y(N15)*state.rho + screened_rates(k_F19_to_He4_N15_derived)*Y(F19)) +
        (-screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(N15)*Y(H1)*state.rho + screened_rates(k_He4_C12_to_p_N15_derived)*Y(C12)*Y(He4)*state.rho) +
        (screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(O18)*Y(H1)*state.rho + -screened_rates(k_He4_N15_to_p_O18_derived)*Y(He4)*Y(N15)*state.rho);

    ydot_nuc(O14) =
        -screened_rates(k_O14_to_N14_reaclib)*Y(O14) +
        (screened_rates(k_p_N13_to_O14_reaclib)*Y(N13)*Y(H1)*state.rho + -screened_rates(k_O14_to_p_N13_derived)*Y(O14)) +
        (-screened_rates(k_He4_O14_to_Ne18_reaclib)*Y(He4)*Y(O14)*state.rho + screened_rates(k_Ne18_to_He4_O14_derived)*Y(Ne18)) +
        (-screened_rates(k_He4_O14_to_p_F17_reaclib)*Y(He4)*Y(O14)*state.rho + screened_rates(k_p_F17_to_He4_O14_derived)*Y(F17)*Y(H1)*state.rho);

    ydot_nuc(O15) =
        -screened_rates(k_O15_to_N15_reaclib)*Y(O15) +
        (screened_rates(k_p_N14_to_O15_reaclib)*Y(N14)*Y(H1)*state.rho + -screened_rates(k_O15_to_p_N14_derived)*Y(O15)) +
        (-screened_rates(k_He4_O15_to_Ne19_reaclib)*Y(He4)*Y(O15)*state.rho + screened_rates(k_Ne19_to_He4_O15_derived)*Y(Ne19)) +
        (screened_rates(k_p_F18_to_He4_O15_reaclib)*Y(F18)*Y(H1)*state.rho + -screened_rates(k_He4_O15_to_p_F18_derived)*Y(He4)*Y(O15)*state.rho);

    ydot_nuc(O16) =
        (screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*Y(He4)*state.rho + -screened_rates(k_O16_to_He4_C12_derived)*Y(O16)) +
        (screened_rates(k_p_N15_to_O16_reaclib)*Y(N15)*Y(H1)*state.rho + -screened_rates(k_O16_to_p_N15_derived)*Y(O16)) +
        (-screened_rates(k_p_O16_to_F17_reaclib)*Y(O16)*Y(H1)*state.rho + screened_rates(k_F17_to_p_O16_derived)*Y(F17)) +
        (-screened_rates(k_He4_O16_to_Ne20_reaclib)*Y(He4)*Y(O16)*state.rho + screened_rates(k_Ne20_to_He4_O16_derived)*Y(Ne20)) +
        (screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(He4)*Y(N13)*state.rho + -screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*Y(H1)*state.rho) +
        -screened_rates(k_C12_O16_to_p_Al27_reaclib)*Y(C12)*Y(O16)*state.rho +
        -screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(C12)*Y(O16)*state.rho +
        -screened_rates(k_O16_O16_to_p_P31_reaclib)*amrex::Math::powi<2>(Y(O16))*state.rho +
        -screened_rates(k_O16_O16_to_He4_Si28_reaclib)*amrex::Math::powi<2>(Y(O16))*state.rho +
        (screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(F19)*Y(H1)*state.rho + -screened_rates(k_He4_O16_to_p_F19_derived)*Y(He4)*Y(O16)*state.rho) +
        -screened_rates(k_O16_O16_to_S32_modified)*amrex::Math::powi<2>(Y(O16))*state.rho +
        -screened_rates(k_C12_O16_to_Si28_modified)*Y(C12)*Y(O16)*state.rho;

    ydot_nuc(O17) =
        (-screened_rates(k_p_O17_to_F18_reaclib)*Y(O17)*Y(H1)*state.rho + screened_rates(k_F18_to_p_O17_derived)*Y(F18)) +
        (-screened_rates(k_He4_O17_to_Ne21_reaclib)*Y(He4)*Y(O17)*state.rho + screened_rates(k_Ne21_to_He4_O17_derived)*Y(Ne21)) +
        (-screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(O17)*Y(H1)*state.rho + screened_rates(k_He4_N14_to_p_O17_derived)*Y(He4)*Y(N14)*state.rho) +
        (screened_rates(k_F17_to_O17_weaktab)*Y(F17) + -screened_rates(k_O17_to_F17_weaktab)*Y(O17));

    ydot_nuc(O18) =
        (-screened_rates(k_p_O18_to_F19_reaclib)*Y(O18)*Y(H1)*state.rho + screened_rates(k_F19_to_p_O18_derived)*Y(F19)) +
        (-screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(O18)*Y(H1)*state.rho + screened_rates(k_He4_N15_to_p_O18_derived)*Y(He4)*Y(N15)*state.rho) +
        (screened_rates(k_F18_to_O18_weaktab)*Y(F18) + -screened_rates(k_O18_to_F18_weaktab)*Y(O18));

    ydot_nuc(F17) =
        (screened_rates(k_p_O16_to_F17_reaclib)*Y(O16)*Y(H1)*state.rho + -screened_rates(k_F17_to_p_O16_derived)*Y(F17)) +
        (-screened_rates(k_p_F17_to_Ne18_reaclib)*Y(F17)*Y(H1)*state.rho + screened_rates(k_Ne18_to_p_F17_derived)*Y(Ne18)) +
        (screened_rates(k_He4_O14_to_p_F17_reaclib)*Y(He4)*Y(O14)*state.rho + -screened_rates(k_p_F17_to_He4_O14_derived)*Y(F17)*Y(H1)*state.rho) +
        (-screened_rates(k_He4_F17_to_p_Ne20_derived)*Y(F17)*Y(He4)*state.rho + screened_rates(k_p_Ne20_to_He4_F17_reaclib)*Y(Ne20)*Y(H1)*state.rho) +
        (-screened_rates(k_F17_to_O17_weaktab)*Y(F17) + screened_rates(k_O17_to_F17_weaktab)*Y(O17));

    ydot_nuc(F18) =
        (screened_rates(k_He4_N14_to_F18_reaclib)*Y(He4)*Y(N14)*state.rho + -screened_rates(k_F18_to_He4_N14_derived)*Y(F18)) +
        (screened_rates(k_p_O17_to_F18_reaclib)*Y(O17)*Y(H1)*state.rho + -screened_rates(k_F18_to_p_O17_derived)*Y(F18)) +
        (-screened_rates(k_p_F18_to_Ne19_reaclib)*Y(F18)*Y(H1)*state.rho + screened_rates(k_Ne19_to_p_F18_derived)*Y(Ne19)) +
        (-screened_rates(k_He4_F18_to_Na22_reaclib)*Y(F18)*Y(He4)*state.rho + screened_rates(k_Na22_to_He4_F18_derived)*Y(Na22)) +
        (-screened_rates(k_p_F18_to_He4_O15_reaclib)*Y(F18)*Y(H1)*state.rho + screened_rates(k_He4_O15_to_p_F18_derived)*Y(He4)*Y(O15)*state.rho) +
        (-screened_rates(k_He4_F18_to_p_Ne21_reaclib)*Y(F18)*Y(He4)*state.rho + screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(Ne21)*Y(H1)*state.rho) +
        (-screened_rates(k_F18_to_O18_weaktab)*Y(F18) + screened_rates(k_O18_to_F18_weaktab)*Y(O18)) +
        (screened_rates(k_Ne18_to_F18_weaktab)*Y(Ne18) + -screened_rates(k_F18_to_Ne18_weaktab)*Y(F18));

    ydot_nuc(F19) =
        (screened_rates(k_He4_N15_to_F19_reaclib)*Y(He4)*Y(N15)*state.rho + -screened_rates(k_F19_to_He4_N15_derived)*Y(F19)) +
        (screened_rates(k_p_O18_to_F19_reaclib)*Y(O18)*Y(H1)*state.rho + -screened_rates(k_F19_to_p_O18_derived)*Y(F19)) +
        (-screened_rates(k_p_F19_to_Ne20_reaclib)*Y(F19)*Y(H1)*state.rho + screened_rates(k_Ne20_to_p_F19_derived)*Y(Ne20)) +
        (-screened_rates(k_He4_F19_to_Na23_reaclib)*Y(F19)*Y(He4)*state.rho + screened_rates(k_Na23_to_He4_F19_derived)*Y(Na23)) +
        (-screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(F19)*Y(H1)*state.rho + screened_rates(k_He4_O16_to_p_F19_derived)*Y(He4)*Y(O16)*state.rho) +
        (screened_rates(k_Ne19_to_F19_weaktab)*Y(Ne19) + -screened_rates(k_F19_to_Ne19_weaktab)*Y(F19));

    ydot_nuc(Ne18) =
        (screened_rates(k_He4_O14_to_Ne18_reaclib)*Y(He4)*Y(O14)*state.rho + -screened_rates(k_Ne18_to_He4_O14_derived)*Y(Ne18)) +
        (screened_rates(k_p_F17_to_Ne18_reaclib)*Y(F17)*Y(H1)*state.rho + -screened_rates(k_Ne18_to_p_F17_derived)*Y(Ne18)) +
        (-screened_rates(k_He4_Ne18_to_Mg22_reaclib)*Y(He4)*Y(Ne18)*state.rho + screened_rates(k_Mg22_to_He4_Ne18_derived)*Y(Mg22)) +
        (-screened_rates(k_Ne18_to_F18_weaktab)*Y(Ne18) + screened_rates(k_F18_to_Ne18_weaktab)*Y(F18));

    ydot_nuc(Ne19) =
        (screened_rates(k_He4_O15_to_Ne19_reaclib)*Y(He4)*Y(O15)*state.rho + -screened_rates(k_Ne19_to_He4_O15_derived)*Y(Ne19)) +
        (screened_rates(k_p_F18_to_Ne19_reaclib)*Y(F18)*Y(H1)*state.rho + -screened_rates(k_Ne19_to_p_F18_derived)*Y(Ne19)) +
        (-screened_rates(k_He4_Ne19_to_p_Na22_reaclib)*Y(He4)*Y(Ne19)*state.rho + screened_rates(k_p_Na22_to_He4_Ne19_derived)*Y(Na22)*Y(H1)*state.rho) +
        (-screened_rates(k_Ne19_to_F19_weaktab)*Y(Ne19) + screened_rates(k_F19_to_Ne19_weaktab)*Y(F19));

    ydot_nuc(Ne20) =
        (screened_rates(k_He4_O16_to_Ne20_reaclib)*Y(He4)*Y(O16)*state.rho + -screened_rates(k_Ne20_to_He4_O16_derived)*Y(Ne20)) +
        (screened_rates(k_p_F19_to_Ne20_reaclib)*Y(F19)*Y(H1)*state.rho + -screened_rates(k_Ne20_to_p_F19_derived)*Y(Ne20)) +
        (-screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(He4)*Y(Ne20)*state.rho + screened_rates(k_Mg24_to_He4_Ne20_derived)*Y(Mg24)) +
        0.5*screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho +
        (screened_rates(k_p_Na23_to_He4_Ne20_reaclib)*Y(Na23)*Y(H1)*state.rho + -screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(He4)*Y(Ne20)*state.rho) +
        (screened_rates(k_He4_F17_to_p_Ne20_derived)*Y(F17)*Y(He4)*state.rho + -screened_rates(k_p_Ne20_to_He4_F17_reaclib)*Y(Ne20)*Y(H1)*state.rho);

    ydot_nuc(Ne21) =
        (screened_rates(k_He4_O17_to_Ne21_reaclib)*Y(He4)*Y(O17)*state.rho + -screened_rates(k_Ne21_to_He4_O17_derived)*Y(Ne21)) +
        (-screened_rates(k_p_Ne21_to_Na22_reaclib)*Y(Ne21)*Y(H1)*state.rho + screened_rates(k_Na22_to_p_Ne21_derived)*Y(Na22)) +
        (screened_rates(k_He4_F18_to_p_Ne21_reaclib)*Y(F18)*Y(He4)*state.rho + -screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(Ne21)*Y(H1)*state.rho);

    ydot_nuc(Na22) =
        (screened_rates(k_He4_F18_to_Na22_reaclib)*Y(F18)*Y(He4)*state.rho + -screened_rates(k_Na22_to_He4_F18_derived)*Y(Na22)) +
        (screened_rates(k_p_Ne21_to_Na22_reaclib)*Y(Ne21)*Y(H1)*state.rho + -screened_rates(k_Na22_to_p_Ne21_derived)*Y(Na22)) +
        (screened_rates(k_He4_Ne19_to_p_Na22_reaclib)*Y(He4)*Y(Ne19)*state.rho + -screened_rates(k_p_Na22_to_He4_Ne19_derived)*Y(Na22)*Y(H1)*state.rho) +
        (screened_rates(k_Mg22_to_Na22_weaktab)*Y(Mg22) + -screened_rates(k_Na22_to_Mg22_weaktab)*Y(Na22));

    ydot_nuc(Na23) =
        (screened_rates(k_He4_F19_to_Na23_reaclib)*Y(F19)*Y(He4)*state.rho + -screened_rates(k_Na23_to_He4_F19_derived)*Y(Na23)) +
        (-screened_rates(k_p_Na23_to_Mg24_reaclib)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_Mg24_to_p_Na23_derived)*Y(Mg24)) +
        0.5*screened_rates(k_C12_C12_to_p_Na23_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho +
        (-screened_rates(k_p_Na23_to_He4_Ne20_reaclib)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(He4)*Y(Ne20)*state.rho);

    ydot_nuc(Mg22) =
        (screened_rates(k_He4_Ne18_to_Mg22_reaclib)*Y(He4)*Y(Ne18)*state.rho + -screened_rates(k_Mg22_to_He4_Ne18_derived)*Y(Mg22)) +
        (-screened_rates(k_Mg22_to_Na22_weaktab)*Y(Mg22) + screened_rates(k_Na22_to_Mg22_weaktab)*Y(Na22));

    ydot_nuc(Mg24) =
        (screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(He4)*Y(Ne20)*state.rho + -screened_rates(k_Mg24_to_He4_Ne20_derived)*Y(Mg24)) +
        (screened_rates(k_p_Na23_to_Mg24_reaclib)*Y(Na23)*Y(H1)*state.rho + -screened_rates(k_Mg24_to_p_Na23_derived)*Y(Mg24)) +
        (-screened_rates(k_He4_Mg24_to_Si28_reaclib)*Y(He4)*Y(Mg24)*state.rho + screened_rates(k_Si28_to_He4_Mg24_derived)*Y(Si28)) +
        screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(C12)*Y(O16)*state.rho +
        (screened_rates(k_p_Al27_to_He4_Mg24_reaclib)*Y(Al27)*Y(H1)*state.rho + -screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(He4)*Y(Mg24)*state.rho) +
        0.5*screened_rates(k_C12_C12_to_Mg24_modified)*amrex::Math::powi<2>(Y(C12))*state.rho;

    ydot_nuc(Al27) =
        (-screened_rates(k_p_Al27_to_Si28_reaclib)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_Si28_to_p_Al27_derived)*Y(Si28)) +
        screened_rates(k_C12_O16_to_p_Al27_reaclib)*Y(C12)*Y(O16)*state.rho +
        (-screened_rates(k_p_Al27_to_He4_Mg24_reaclib)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(He4)*Y(Mg24)*state.rho);

    ydot_nuc(Si28) =
        (screened_rates(k_He4_Mg24_to_Si28_reaclib)*Y(He4)*Y(Mg24)*state.rho + -screened_rates(k_Si28_to_He4_Mg24_derived)*Y(Si28)) +
        (screened_rates(k_p_Al27_to_Si28_reaclib)*Y(Al27)*Y(H1)*state.rho + -screened_rates(k_Si28_to_p_Al27_derived)*Y(Si28)) +
        (-screened_rates(k_He4_Si28_to_S32_reaclib)*Y(He4)*Y(Si28)*state.rho + screened_rates(k_S32_to_He4_Si28_derived)*Y(S32)) +
        0.5*screened_rates(k_O16_O16_to_He4_Si28_reaclib)*amrex::Math::powi<2>(Y(O16))*state.rho +
        (screened_rates(k_p_P31_to_He4_Si28_reaclib)*Y(P31)*Y(H1)*state.rho + -screened_rates(k_He4_Si28_to_p_P31_derived)*Y(He4)*Y(Si28)*state.rho) +
        screened_rates(k_C12_O16_to_Si28_modified)*Y(C12)*Y(O16)*state.rho;

    ydot_nuc(P31) =
        (-screened_rates(k_p_P31_to_S32_reaclib)*Y(P31)*Y(H1)*state.rho + screened_rates(k_S32_to_p_P31_derived)*Y(S32)) +
        0.5*screened_rates(k_O16_O16_to_p_P31_reaclib)*amrex::Math::powi<2>(Y(O16))*state.rho +
        (-screened_rates(k_p_P31_to_He4_Si28_reaclib)*Y(P31)*Y(H1)*state.rho + screened_rates(k_He4_Si28_to_p_P31_derived)*Y(He4)*Y(Si28)*state.rho);

    ydot_nuc(S32) =
        (screened_rates(k_He4_Si28_to_S32_reaclib)*Y(He4)*Y(Si28)*state.rho + -screened_rates(k_S32_to_He4_Si28_derived)*Y(S32)) +
        (screened_rates(k_p_P31_to_S32_reaclib)*Y(P31)*Y(H1)*state.rho + -screened_rates(k_S32_to_p_P31_derived)*Y(S32)) +
        0.5*screened_rates(k_O16_O16_to_S32_modified)*amrex::Math::powi<2>(Y(O16))*state.rho +
        (-screened_rates(k_S32_He4_to_Ar36_approx)*Y(He4)*Y(S32)*state.rho + screened_rates(k_Ar36_to_S32_He4_approx)*Y(Ar36));

    ydot_nuc(Ar36) =
        (screened_rates(k_S32_He4_to_Ar36_approx)*Y(He4)*Y(S32)*state.rho + -screened_rates(k_Ar36_to_S32_He4_approx)*Y(Ar36)) +
        (-screened_rates(k_Ar36_He4_to_Ca40_approx)*Y(Ar36)*Y(He4)*state.rho + screened_rates(k_Ca40_to_Ar36_He4_approx)*Y(Ca40));

    ydot_nuc(Ca40) =
        (screened_rates(k_Ar36_He4_to_Ca40_approx)*Y(Ar36)*Y(He4)*state.rho + -screened_rates(k_Ca40_to_Ar36_He4_approx)*Y(Ca40)) +
        (-screened_rates(k_Ca40_He4_to_Ti44_approx)*Y(Ca40)*Y(He4)*state.rho + screened_rates(k_Ti44_to_Ca40_He4_approx)*Y(Ti44));

    ydot_nuc(Ti44) =
        (screened_rates(k_Ca40_He4_to_Ti44_approx)*Y(Ca40)*Y(He4)*state.rho + -screened_rates(k_Ti44_to_Ca40_He4_approx)*Y(Ti44)) +
        (-screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(He4)*Y(Ti44)*state.rho + screened_rates(k_Cr48_to_Ti44_He4_approx)*Y(Cr48));

    ydot_nuc(Cr48) =
        (screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(He4)*Y(Ti44)*state.rho + -screened_rates(k_Cr48_to_Ti44_He4_approx)*Y(Cr48)) +
        (-screened_rates(k_Cr48_He4_to_Fe52_approx)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_Fe52_to_Cr48_He4_approx)*Y(Fe52));

    ydot_nuc(Fe52) =
        (screened_rates(k_Cr48_He4_to_Fe52_approx)*Y(Cr48)*Y(He4)*state.rho + -screened_rates(k_Fe52_to_Cr48_He4_approx)*Y(Fe52)) +
        (-screened_rates(k_Fe52_He4_to_Ni56_approx)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_Ni56_to_Fe52_He4_approx)*Y(Ni56));

    ydot_nuc(Ni56) =
        (screened_rates(k_Fe52_He4_to_Ni56_approx)*Y(Fe52)*Y(He4)*state.rho + -screened_rates(k_Ni56_to_Fe52_He4_approx)*Y(Ni56));

}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
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
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void jac_nuc(const burn_t& state,
             MatrixType& jac,
             const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
             const amrex::Array1D<amrex::Real, 1, NumRates>& screened_rates)
{

    amrex::Real scratch;

    scratch = -screened_rates(k_p_Al27_to_He4_Mg24_reaclib)*Y(Al27)*state.rho - screened_rates(k_p_Al27_to_Si28_reaclib)*Y(Al27)*state.rho - screened_rates(k_p_C12_to_N13_reaclib)*Y(C12)*state.rho - screened_rates(k_p_C13_to_N14_reaclib)*Y(C13)*state.rho - screened_rates(k_p_F17_to_He4_O14_derived)*Y(F17)*state.rho - screened_rates(k_p_F17_to_Ne18_reaclib)*Y(F17)*state.rho - screened_rates(k_p_F18_to_He4_O15_reaclib)*Y(F18)*state.rho - screened_rates(k_p_F18_to_Ne19_reaclib)*Y(F18)*state.rho - screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(F19)*state.rho - screened_rates(k_p_F19_to_Ne20_reaclib)*Y(F19)*state.rho - screened_rates(k_p_N13_to_O14_reaclib)*Y(N13)*state.rho - screened_rates(k_p_N14_to_O15_reaclib)*Y(N14)*state.rho - screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(N15)*state.rho - screened_rates(k_p_N15_to_O16_reaclib)*Y(N15)*state.rho - screened_rates(k_p_Na22_to_He4_Ne19_derived)*Y(Na22)*state.rho - screened_rates(k_p_Na23_to_He4_Ne20_reaclib)*Y(Na23)*state.rho - screened_rates(k_p_Na23_to_Mg24_reaclib)*Y(Na23)*state.rho - screened_rates(k_p_Ne20_to_He4_F17_reaclib)*Y(Ne20)*state.rho - screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(Ne21)*state.rho - screened_rates(k_p_Ne21_to_Na22_reaclib)*Y(Ne21)*state.rho - screened_rates(k_p_O16_to_F17_reaclib)*Y(O16)*state.rho - screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*state.rho - screened_rates(k_p_O17_to_F18_reaclib)*Y(O17)*state.rho - screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(O17)*state.rho - screened_rates(k_p_O18_to_F19_reaclib)*Y(O18)*state.rho - screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(O18)*state.rho - screened_rates(k_p_P31_to_He4_Si28_reaclib)*Y(P31)*state.rho - screened_rates(k_p_P31_to_S32_reaclib)*Y(P31)*state.rho;
    jac.set(H1, H1, scratch);

    scratch = screened_rates(k_He4_C12_to_p_N15_derived)*Y(C12)*state.rho + screened_rates(k_He4_F17_to_p_Ne20_derived)*Y(F17)*state.rho + screened_rates(k_He4_F18_to_p_Ne21_reaclib)*Y(F18)*state.rho + screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(Mg24)*state.rho + screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(N13)*state.rho + screened_rates(k_He4_N14_to_p_O17_derived)*Y(N14)*state.rho + screened_rates(k_He4_N15_to_p_O18_derived)*Y(N15)*state.rho + screened_rates(k_He4_Ne19_to_p_Na22_reaclib)*Y(Ne19)*state.rho + screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(Ne20)*state.rho + screened_rates(k_He4_O14_to_p_F17_reaclib)*Y(O14)*state.rho + screened_rates(k_He4_O15_to_p_F18_derived)*Y(O15)*state.rho + screened_rates(k_He4_O16_to_p_F19_derived)*Y(O16)*state.rho + screened_rates(k_He4_Si28_to_p_P31_derived)*Y(Si28)*state.rho;
    jac.set(H1, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_p_Na23_reaclib)*Y(C12)*state.rho + screened_rates(k_C12_O16_to_p_Al27_reaclib)*Y(O16)*state.rho + screened_rates(k_He4_C12_to_p_N15_derived)*Y(He4)*state.rho - screened_rates(k_p_C12_to_N13_reaclib)*Y(H1)*state.rho;
    jac.set(H1, C12, scratch);

    scratch = -screened_rates(k_p_C13_to_N14_reaclib)*Y(H1)*state.rho;
    jac.set(H1, C13, scratch);

    scratch = screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(He4)*state.rho + screened_rates(k_N13_to_p_C12_derived) - screened_rates(k_p_N13_to_O14_reaclib)*Y(H1)*state.rho;
    jac.set(H1, N13, scratch);

    scratch = screened_rates(k_He4_N14_to_p_O17_derived)*Y(He4)*state.rho + screened_rates(k_N14_to_p_C13_derived) - screened_rates(k_p_N14_to_O15_reaclib)*Y(H1)*state.rho;
    jac.set(H1, N14, scratch);

    scratch = screened_rates(k_He4_N15_to_p_O18_derived)*Y(He4)*state.rho - screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(H1)*state.rho - screened_rates(k_p_N15_to_O16_reaclib)*Y(H1)*state.rho;
    jac.set(H1, N15, scratch);

    scratch = screened_rates(k_He4_O14_to_p_F17_reaclib)*Y(He4)*state.rho + screened_rates(k_O14_to_p_N13_derived);
    jac.set(H1, O14, scratch);

    scratch = screened_rates(k_He4_O15_to_p_F18_derived)*Y(He4)*state.rho + screened_rates(k_O15_to_p_N14_derived);
    jac.set(H1, O15, scratch);

    scratch = screened_rates(k_C12_O16_to_p_Al27_reaclib)*Y(C12)*state.rho + screened_rates(k_He4_O16_to_p_F19_derived)*Y(He4)*state.rho + 1.0*screened_rates(k_O16_O16_to_p_P31_reaclib)*Y(O16)*state.rho + screened_rates(k_O16_to_p_N15_derived) - screened_rates(k_p_O16_to_F17_reaclib)*Y(H1)*state.rho - screened_rates(k_p_O16_to_He4_N13_derived)*Y(H1)*state.rho;
    jac.set(H1, O16, scratch);

    scratch = -screened_rates(k_p_O17_to_F18_reaclib)*Y(H1)*state.rho - screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(H1)*state.rho;
    jac.set(H1, O17, scratch);

    scratch = -screened_rates(k_p_O18_to_F19_reaclib)*Y(H1)*state.rho - screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(H1)*state.rho;
    jac.set(H1, O18, scratch);

    scratch = screened_rates(k_F17_to_p_O16_derived) + screened_rates(k_He4_F17_to_p_Ne20_derived)*Y(He4)*state.rho - screened_rates(k_p_F17_to_He4_O14_derived)*Y(H1)*state.rho - screened_rates(k_p_F17_to_Ne18_reaclib)*Y(H1)*state.rho;
    jac.set(H1, F17, scratch);

    scratch = screened_rates(k_F18_to_p_O17_derived) + screened_rates(k_He4_F18_to_p_Ne21_reaclib)*Y(He4)*state.rho - screened_rates(k_p_F18_to_He4_O15_reaclib)*Y(H1)*state.rho - screened_rates(k_p_F18_to_Ne19_reaclib)*Y(H1)*state.rho;
    jac.set(H1, F18, scratch);

    scratch = screened_rates(k_F19_to_p_O18_derived) - screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(H1)*state.rho - screened_rates(k_p_F19_to_Ne20_reaclib)*Y(H1)*state.rho;
    jac.set(H1, F19, scratch);

    scratch = screened_rates(k_Ne18_to_p_F17_derived);
    jac.set(H1, Ne18, scratch);

    scratch = screened_rates(k_He4_Ne19_to_p_Na22_reaclib)*Y(He4)*state.rho + screened_rates(k_Ne19_to_p_F18_derived);
    jac.set(H1, Ne19, scratch);

    scratch = screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(He4)*state.rho + screened_rates(k_Ne20_to_p_F19_derived) - screened_rates(k_p_Ne20_to_He4_F17_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Ne20, scratch);

    scratch = -screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(H1)*state.rho - screened_rates(k_p_Ne21_to_Na22_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Ne21, scratch);

    scratch = screened_rates(k_Na22_to_p_Ne21_derived) - screened_rates(k_p_Na22_to_He4_Ne19_derived)*Y(H1)*state.rho;
    jac.set(H1, Na22, scratch);

    scratch = -screened_rates(k_p_Na23_to_He4_Ne20_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Na23_to_Mg24_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Na23, scratch);

    jac.set(H1, Mg22, 0.0);

    scratch = screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(He4)*state.rho + screened_rates(k_Mg24_to_p_Na23_derived);
    jac.set(H1, Mg24, scratch);

    scratch = -screened_rates(k_p_Al27_to_He4_Mg24_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Al27_to_Si28_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Al27, scratch);

    scratch = screened_rates(k_He4_Si28_to_p_P31_derived)*Y(He4)*state.rho + screened_rates(k_Si28_to_p_Al27_derived);
    jac.set(H1, Si28, scratch);

    scratch = -screened_rates(k_p_P31_to_He4_Si28_reaclib)*Y(H1)*state.rho - screened_rates(k_p_P31_to_S32_reaclib)*Y(H1)*state.rho;
    jac.set(H1, P31, scratch);

    scratch = screened_rates(k_S32_to_p_P31_derived);
    jac.set(H1, S32, scratch);

    jac.set(H1, Ar36, 0.0);

    jac.set(H1, Ca40, 0.0);

    jac.set(H1, Ti44, 0.0);

    jac.set(H1, Cr48, 0.0);

    jac.set(H1, Fe52, 0.0);

    jac.set(H1, Ni56, 0.0);

    scratch = screened_rates(k_p_Al27_to_He4_Mg24_reaclib)*Y(Al27)*state.rho + screened_rates(k_p_F17_to_He4_O14_derived)*Y(F17)*state.rho + screened_rates(k_p_F18_to_He4_O15_reaclib)*Y(F18)*state.rho + screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(F19)*state.rho + screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(N15)*state.rho + screened_rates(k_p_Na22_to_He4_Ne19_derived)*Y(Na22)*state.rho + screened_rates(k_p_Na23_to_He4_Ne20_reaclib)*Y(Na23)*state.rho + screened_rates(k_p_Ne20_to_He4_F17_reaclib)*Y(Ne20)*state.rho + screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(Ne21)*state.rho + screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*state.rho + screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(O17)*state.rho + screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(O18)*state.rho + screened_rates(k_p_P31_to_He4_Si28_reaclib)*Y(P31)*state.rho;
    jac.set(He4, H1, scratch);

    scratch = -screened_rates(k_Ar36_He4_to_Ca40_approx)*Y(Ar36)*state.rho - screened_rates(k_Ca40_He4_to_Ti44_approx)*Y(Ca40)*state.rho - screened_rates(k_Cr48_He4_to_Fe52_approx)*Y(Cr48)*state.rho - screened_rates(k_Fe52_He4_to_Ni56_approx)*Y(Fe52)*state.rho - screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*state.rho - screened_rates(k_He4_C12_to_p_N15_derived)*Y(C12)*state.rho - screened_rates(k_He4_F17_to_p_Ne20_derived)*Y(F17)*state.rho - screened_rates(k_He4_F18_to_Na22_reaclib)*Y(F18)*state.rho - screened_rates(k_He4_F18_to_p_Ne21_reaclib)*Y(F18)*state.rho - screened_rates(k_He4_F19_to_Na23_reaclib)*Y(F19)*state.rho - 1.5*screened_rates(k_He4_He4_He4_to_C12_reaclib)*amrex::Math::powi<2>(Y(He4))*amrex::Math::powi<2>(state.rho) - screened_rates(k_He4_Mg24_to_Si28_reaclib)*Y(Mg24)*state.rho - screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(Mg24)*state.rho - screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(N13)*state.rho - screened_rates(k_He4_N14_to_F18_reaclib)*Y(N14)*state.rho - screened_rates(k_He4_N14_to_p_O17_derived)*Y(N14)*state.rho - screened_rates(k_He4_N15_to_F19_reaclib)*Y(N15)*state.rho - screened_rates(k_He4_N15_to_p_O18_derived)*Y(N15)*state.rho - screened_rates(k_He4_Ne18_to_Mg22_reaclib)*Y(Ne18)*state.rho - screened_rates(k_He4_Ne19_to_p_Na22_reaclib)*Y(Ne19)*state.rho - screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(Ne20)*state.rho - screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(Ne20)*state.rho - screened_rates(k_He4_O14_to_Ne18_reaclib)*Y(O14)*state.rho - screened_rates(k_He4_O14_to_p_F17_reaclib)*Y(O14)*state.rho - screened_rates(k_He4_O15_to_Ne19_reaclib)*Y(O15)*state.rho - screened_rates(k_He4_O15_to_p_F18_derived)*Y(O15)*state.rho - screened_rates(k_He4_O16_to_Ne20_reaclib)*Y(O16)*state.rho - screened_rates(k_He4_O16_to_p_F19_derived)*Y(O16)*state.rho - screened_rates(k_He4_O17_to_Ne21_reaclib)*Y(O17)*state.rho - screened_rates(k_He4_Si28_to_S32_reaclib)*Y(Si28)*state.rho - screened_rates(k_He4_Si28_to_p_P31_derived)*Y(Si28)*state.rho - screened_rates(k_S32_He4_to_Ar36_approx)*Y(S32)*state.rho - screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(Ti44)*state.rho;
    jac.set(He4, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*Y(C12)*state.rho + screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(O16)*state.rho + 3.0*screened_rates(k_C12_to_He4_He4_He4_derived) - screened_rates(k_He4_C12_to_O16_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_C12_to_p_N15_derived)*Y(He4)*state.rho;
    jac.set(He4, C12, scratch);

    jac.set(He4, C13, 0.0);

    scratch = -screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(He4)*state.rho;
    jac.set(He4, N13, scratch);

    scratch = -screened_rates(k_He4_N14_to_F18_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_N14_to_p_O17_derived)*Y(He4)*state.rho;
    jac.set(He4, N14, scratch);

    scratch = -screened_rates(k_He4_N15_to_F19_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_N15_to_p_O18_derived)*Y(He4)*state.rho + screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(H1)*state.rho;
    jac.set(He4, N15, scratch);

    scratch = -screened_rates(k_He4_O14_to_Ne18_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_O14_to_p_F17_reaclib)*Y(He4)*state.rho;
    jac.set(He4, O14, scratch);

    scratch = -screened_rates(k_He4_O15_to_Ne19_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_O15_to_p_F18_derived)*Y(He4)*state.rho;
    jac.set(He4, O15, scratch);

    scratch = screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(C12)*state.rho - screened_rates(k_He4_O16_to_Ne20_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_O16_to_p_F19_derived)*Y(He4)*state.rho + 1.0*screened_rates(k_O16_O16_to_He4_Si28_reaclib)*Y(O16)*state.rho + screened_rates(k_O16_to_He4_C12_derived) + screened_rates(k_p_O16_to_He4_N13_derived)*Y(H1)*state.rho;
    jac.set(He4, O16, scratch);

    scratch = -screened_rates(k_He4_O17_to_Ne21_reaclib)*Y(He4)*state.rho + screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(H1)*state.rho;
    jac.set(He4, O17, scratch);

    scratch = screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(H1)*state.rho;
    jac.set(He4, O18, scratch);

    scratch = -screened_rates(k_He4_F17_to_p_Ne20_derived)*Y(He4)*state.rho + screened_rates(k_p_F17_to_He4_O14_derived)*Y(H1)*state.rho;
    jac.set(He4, F17, scratch);

    scratch = screened_rates(k_F18_to_He4_N14_derived) - screened_rates(k_He4_F18_to_Na22_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_F18_to_p_Ne21_reaclib)*Y(He4)*state.rho + screened_rates(k_p_F18_to_He4_O15_reaclib)*Y(H1)*state.rho;
    jac.set(He4, F18, scratch);

    scratch = screened_rates(k_F19_to_He4_N15_derived) - screened_rates(k_He4_F19_to_Na23_reaclib)*Y(He4)*state.rho + screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(H1)*state.rho;
    jac.set(He4, F19, scratch);

    scratch = -screened_rates(k_He4_Ne18_to_Mg22_reaclib)*Y(He4)*state.rho + screened_rates(k_Ne18_to_He4_O14_derived);
    jac.set(He4, Ne18, scratch);

    scratch = -screened_rates(k_He4_Ne19_to_p_Na22_reaclib)*Y(He4)*state.rho + screened_rates(k_Ne19_to_He4_O15_derived);
    jac.set(He4, Ne19, scratch);

    scratch = -screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(He4)*state.rho + screened_rates(k_Ne20_to_He4_O16_derived) + screened_rates(k_p_Ne20_to_He4_F17_reaclib)*Y(H1)*state.rho;
    jac.set(He4, Ne20, scratch);

    scratch = screened_rates(k_Ne21_to_He4_O17_derived) + screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(H1)*state.rho;
    jac.set(He4, Ne21, scratch);

    scratch = screened_rates(k_Na22_to_He4_F18_derived) + screened_rates(k_p_Na22_to_He4_Ne19_derived)*Y(H1)*state.rho;
    jac.set(He4, Na22, scratch);

    scratch = screened_rates(k_Na23_to_He4_F19_derived) + screened_rates(k_p_Na23_to_He4_Ne20_reaclib)*Y(H1)*state.rho;
    jac.set(He4, Na23, scratch);

    scratch = screened_rates(k_Mg22_to_He4_Ne18_derived);
    jac.set(He4, Mg22, scratch);

    scratch = -screened_rates(k_He4_Mg24_to_Si28_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(He4)*state.rho + screened_rates(k_Mg24_to_He4_Ne20_derived);
    jac.set(He4, Mg24, scratch);

    scratch = screened_rates(k_p_Al27_to_He4_Mg24_reaclib)*Y(H1)*state.rho;
    jac.set(He4, Al27, scratch);

    scratch = -screened_rates(k_He4_Si28_to_S32_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Si28_to_p_P31_derived)*Y(He4)*state.rho + screened_rates(k_Si28_to_He4_Mg24_derived);
    jac.set(He4, Si28, scratch);

    scratch = screened_rates(k_p_P31_to_He4_Si28_reaclib)*Y(H1)*state.rho;
    jac.set(He4, P31, scratch);

    scratch = -screened_rates(k_S32_He4_to_Ar36_approx)*Y(He4)*state.rho + screened_rates(k_S32_to_He4_Si28_derived);
    jac.set(He4, S32, scratch);

    scratch = -screened_rates(k_Ar36_He4_to_Ca40_approx)*Y(He4)*state.rho + screened_rates(k_Ar36_to_S32_He4_approx);
    jac.set(He4, Ar36, scratch);

    scratch = -screened_rates(k_Ca40_He4_to_Ti44_approx)*Y(He4)*state.rho + screened_rates(k_Ca40_to_Ar36_He4_approx);
    jac.set(He4, Ca40, scratch);

    scratch = -screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(He4)*state.rho + screened_rates(k_Ti44_to_Ca40_He4_approx);
    jac.set(He4, Ti44, scratch);

    scratch = -screened_rates(k_Cr48_He4_to_Fe52_approx)*Y(He4)*state.rho + screened_rates(k_Cr48_to_Ti44_He4_approx);
    jac.set(He4, Cr48, scratch);

    scratch = -screened_rates(k_Fe52_He4_to_Ni56_approx)*Y(He4)*state.rho + screened_rates(k_Fe52_to_Cr48_He4_approx);
    jac.set(He4, Fe52, scratch);

    scratch = screened_rates(k_Ni56_to_Fe52_He4_approx);
    jac.set(He4, Ni56, scratch);

    scratch = -screened_rates(k_p_C12_to_N13_reaclib)*Y(C12)*state.rho + screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(N15)*state.rho;
    jac.set(C12, H1, scratch);

    scratch = -screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*state.rho - screened_rates(k_He4_C12_to_p_N15_derived)*Y(C12)*state.rho + 0.5*screened_rates(k_He4_He4_He4_to_C12_reaclib)*amrex::Math::powi<2>(Y(He4))*amrex::Math::powi<2>(state.rho);
    jac.set(C12, He4, scratch);

    scratch = -2.0*screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*Y(C12)*state.rho - 2.0*screened_rates(k_C12_C12_to_Mg24_modified)*Y(C12)*state.rho - 2.0*screened_rates(k_C12_C12_to_p_Na23_reaclib)*Y(C12)*state.rho - screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(O16)*state.rho - screened_rates(k_C12_O16_to_Si28_modified)*Y(O16)*state.rho - screened_rates(k_C12_O16_to_p_Al27_reaclib)*Y(O16)*state.rho - screened_rates(k_C12_to_He4_He4_He4_derived) - screened_rates(k_He4_C12_to_O16_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_C12_to_p_N15_derived)*Y(He4)*state.rho - screened_rates(k_p_C12_to_N13_reaclib)*Y(H1)*state.rho;
    jac.set(C12, C12, scratch);

    jac.set(C12, C13, 0.0);

    scratch = screened_rates(k_N13_to_p_C12_derived);
    jac.set(C12, N13, scratch);

    jac.set(C12, N14, 0.0);

    scratch = screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(H1)*state.rho;
    jac.set(C12, N15, scratch);

    jac.set(C12, O14, 0.0);

    jac.set(C12, O15, 0.0);

    scratch = -screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(C12)*state.rho - screened_rates(k_C12_O16_to_Si28_modified)*Y(C12)*state.rho - screened_rates(k_C12_O16_to_p_Al27_reaclib)*Y(C12)*state.rho + screened_rates(k_O16_to_He4_C12_derived);
    jac.set(C12, O16, scratch);

    jac.set(C12, O17, 0.0);

    jac.set(C12, O18, 0.0);

    jac.set(C12, F17, 0.0);

    jac.set(C12, F18, 0.0);

    jac.set(C12, F19, 0.0);

    jac.set(C12, Ne18, 0.0);

    jac.set(C12, Ne19, 0.0);

    jac.set(C12, Ne20, 0.0);

    jac.set(C12, Ne21, 0.0);

    jac.set(C12, Na22, 0.0);

    jac.set(C12, Na23, 0.0);

    jac.set(C12, Mg22, 0.0);

    jac.set(C12, Mg24, 0.0);

    jac.set(C12, Al27, 0.0);

    jac.set(C12, Si28, 0.0);

    jac.set(C12, P31, 0.0);

    jac.set(C12, S32, 0.0);

    jac.set(C12, Ar36, 0.0);

    jac.set(C12, Ca40, 0.0);

    jac.set(C12, Ti44, 0.0);

    jac.set(C12, Cr48, 0.0);

    jac.set(C12, Fe52, 0.0);

    jac.set(C12, Ni56, 0.0);

    scratch = -screened_rates(k_p_C13_to_N14_reaclib)*Y(C13)*state.rho;
    jac.set(C13, H1, scratch);

    jac.set(C13, He4, 0.0);

    jac.set(C13, C12, 0.0);

    scratch = -screened_rates(k_p_C13_to_N14_reaclib)*Y(H1)*state.rho;
    jac.set(C13, C13, scratch);

    scratch = screened_rates(k_N13_to_C13_reaclib);
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

    jac.set(C13, Ne21, 0.0);

    jac.set(C13, Na22, 0.0);

    jac.set(C13, Na23, 0.0);

    jac.set(C13, Mg22, 0.0);

    jac.set(C13, Mg24, 0.0);

    jac.set(C13, Al27, 0.0);

    jac.set(C13, Si28, 0.0);

    jac.set(C13, P31, 0.0);

    jac.set(C13, S32, 0.0);

    jac.set(C13, Ar36, 0.0);

    jac.set(C13, Ca40, 0.0);

    jac.set(C13, Ti44, 0.0);

    jac.set(C13, Cr48, 0.0);

    jac.set(C13, Fe52, 0.0);

    jac.set(C13, Ni56, 0.0);

    scratch = screened_rates(k_p_C12_to_N13_reaclib)*Y(C12)*state.rho - screened_rates(k_p_N13_to_O14_reaclib)*Y(N13)*state.rho + screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*state.rho;
    jac.set(N13, H1, scratch);

    scratch = -screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(N13)*state.rho;
    jac.set(N13, He4, scratch);

    scratch = screened_rates(k_p_C12_to_N13_reaclib)*Y(H1)*state.rho;
    jac.set(N13, C12, scratch);

    jac.set(N13, C13, 0.0);

    scratch = -screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(He4)*state.rho - screened_rates(k_N13_to_C13_reaclib) - screened_rates(k_N13_to_p_C12_derived) - screened_rates(k_p_N13_to_O14_reaclib)*Y(H1)*state.rho;
    jac.set(N13, N13, scratch);

    jac.set(N13, N14, 0.0);

    jac.set(N13, N15, 0.0);

    scratch = screened_rates(k_O14_to_p_N13_derived);
    jac.set(N13, O14, scratch);

    jac.set(N13, O15, 0.0);

    scratch = screened_rates(k_p_O16_to_He4_N13_derived)*Y(H1)*state.rho;
    jac.set(N13, O16, scratch);

    jac.set(N13, O17, 0.0);

    jac.set(N13, O18, 0.0);

    jac.set(N13, F17, 0.0);

    jac.set(N13, F18, 0.0);

    jac.set(N13, F19, 0.0);

    jac.set(N13, Ne18, 0.0);

    jac.set(N13, Ne19, 0.0);

    jac.set(N13, Ne20, 0.0);

    jac.set(N13, Ne21, 0.0);

    jac.set(N13, Na22, 0.0);

    jac.set(N13, Na23, 0.0);

    jac.set(N13, Mg22, 0.0);

    jac.set(N13, Mg24, 0.0);

    jac.set(N13, Al27, 0.0);

    jac.set(N13, Si28, 0.0);

    jac.set(N13, P31, 0.0);

    jac.set(N13, S32, 0.0);

    jac.set(N13, Ar36, 0.0);

    jac.set(N13, Ca40, 0.0);

    jac.set(N13, Ti44, 0.0);

    jac.set(N13, Cr48, 0.0);

    jac.set(N13, Fe52, 0.0);

    jac.set(N13, Ni56, 0.0);

    scratch = screened_rates(k_p_C13_to_N14_reaclib)*Y(C13)*state.rho - screened_rates(k_p_N14_to_O15_reaclib)*Y(N14)*state.rho + screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(O17)*state.rho;
    jac.set(N14, H1, scratch);

    scratch = -screened_rates(k_He4_N14_to_F18_reaclib)*Y(N14)*state.rho - screened_rates(k_He4_N14_to_p_O17_derived)*Y(N14)*state.rho;
    jac.set(N14, He4, scratch);

    jac.set(N14, C12, 0.0);

    scratch = screened_rates(k_p_C13_to_N14_reaclib)*Y(H1)*state.rho;
    jac.set(N14, C13, scratch);

    jac.set(N14, N13, 0.0);

    scratch = -screened_rates(k_He4_N14_to_F18_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_N14_to_p_O17_derived)*Y(He4)*state.rho - screened_rates(k_N14_to_p_C13_derived) - screened_rates(k_p_N14_to_O15_reaclib)*Y(H1)*state.rho;
    jac.set(N14, N14, scratch);

    jac.set(N14, N15, 0.0);

    scratch = screened_rates(k_O14_to_N14_reaclib);
    jac.set(N14, O14, scratch);

    scratch = screened_rates(k_O15_to_p_N14_derived);
    jac.set(N14, O15, scratch);

    jac.set(N14, O16, 0.0);

    scratch = screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(H1)*state.rho;
    jac.set(N14, O17, scratch);

    jac.set(N14, O18, 0.0);

    jac.set(N14, F17, 0.0);

    scratch = screened_rates(k_F18_to_He4_N14_derived);
    jac.set(N14, F18, scratch);

    jac.set(N14, F19, 0.0);

    jac.set(N14, Ne18, 0.0);

    jac.set(N14, Ne19, 0.0);

    jac.set(N14, Ne20, 0.0);

    jac.set(N14, Ne21, 0.0);

    jac.set(N14, Na22, 0.0);

    jac.set(N14, Na23, 0.0);

    jac.set(N14, Mg22, 0.0);

    jac.set(N14, Mg24, 0.0);

    jac.set(N14, Al27, 0.0);

    jac.set(N14, Si28, 0.0);

    jac.set(N14, P31, 0.0);

    jac.set(N14, S32, 0.0);

    jac.set(N14, Ar36, 0.0);

    jac.set(N14, Ca40, 0.0);

    jac.set(N14, Ti44, 0.0);

    jac.set(N14, Cr48, 0.0);

    jac.set(N14, Fe52, 0.0);

    jac.set(N14, Ni56, 0.0);

    scratch = -screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(N15)*state.rho - screened_rates(k_p_N15_to_O16_reaclib)*Y(N15)*state.rho + screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(O18)*state.rho;
    jac.set(N15, H1, scratch);

    scratch = screened_rates(k_He4_C12_to_p_N15_derived)*Y(C12)*state.rho - screened_rates(k_He4_N15_to_F19_reaclib)*Y(N15)*state.rho - screened_rates(k_He4_N15_to_p_O18_derived)*Y(N15)*state.rho;
    jac.set(N15, He4, scratch);

    scratch = screened_rates(k_He4_C12_to_p_N15_derived)*Y(He4)*state.rho;
    jac.set(N15, C12, scratch);

    jac.set(N15, C13, 0.0);

    jac.set(N15, N13, 0.0);

    jac.set(N15, N14, 0.0);

    scratch = -screened_rates(k_He4_N15_to_F19_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_N15_to_p_O18_derived)*Y(He4)*state.rho - screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(H1)*state.rho - screened_rates(k_p_N15_to_O16_reaclib)*Y(H1)*state.rho;
    jac.set(N15, N15, scratch);

    jac.set(N15, O14, 0.0);

    scratch = screened_rates(k_O15_to_N15_reaclib);
    jac.set(N15, O15, scratch);

    scratch = screened_rates(k_O16_to_p_N15_derived);
    jac.set(N15, O16, scratch);

    jac.set(N15, O17, 0.0);

    scratch = screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(H1)*state.rho;
    jac.set(N15, O18, scratch);

    jac.set(N15, F17, 0.0);

    jac.set(N15, F18, 0.0);

    scratch = screened_rates(k_F19_to_He4_N15_derived);
    jac.set(N15, F19, scratch);

    jac.set(N15, Ne18, 0.0);

    jac.set(N15, Ne19, 0.0);

    jac.set(N15, Ne20, 0.0);

    jac.set(N15, Ne21, 0.0);

    jac.set(N15, Na22, 0.0);

    jac.set(N15, Na23, 0.0);

    jac.set(N15, Mg22, 0.0);

    jac.set(N15, Mg24, 0.0);

    jac.set(N15, Al27, 0.0);

    jac.set(N15, Si28, 0.0);

    jac.set(N15, P31, 0.0);

    jac.set(N15, S32, 0.0);

    jac.set(N15, Ar36, 0.0);

    jac.set(N15, Ca40, 0.0);

    jac.set(N15, Ti44, 0.0);

    jac.set(N15, Cr48, 0.0);

    jac.set(N15, Fe52, 0.0);

    jac.set(N15, Ni56, 0.0);

    scratch = screened_rates(k_p_F17_to_He4_O14_derived)*Y(F17)*state.rho + screened_rates(k_p_N13_to_O14_reaclib)*Y(N13)*state.rho;
    jac.set(O14, H1, scratch);

    scratch = -screened_rates(k_He4_O14_to_Ne18_reaclib)*Y(O14)*state.rho - screened_rates(k_He4_O14_to_p_F17_reaclib)*Y(O14)*state.rho;
    jac.set(O14, He4, scratch);

    jac.set(O14, C12, 0.0);

    jac.set(O14, C13, 0.0);

    scratch = screened_rates(k_p_N13_to_O14_reaclib)*Y(H1)*state.rho;
    jac.set(O14, N13, scratch);

    jac.set(O14, N14, 0.0);

    jac.set(O14, N15, 0.0);

    scratch = -screened_rates(k_He4_O14_to_Ne18_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_O14_to_p_F17_reaclib)*Y(He4)*state.rho - screened_rates(k_O14_to_N14_reaclib) - screened_rates(k_O14_to_p_N13_derived);
    jac.set(O14, O14, scratch);

    jac.set(O14, O15, 0.0);

    jac.set(O14, O16, 0.0);

    jac.set(O14, O17, 0.0);

    jac.set(O14, O18, 0.0);

    scratch = screened_rates(k_p_F17_to_He4_O14_derived)*Y(H1)*state.rho;
    jac.set(O14, F17, scratch);

    jac.set(O14, F18, 0.0);

    jac.set(O14, F19, 0.0);

    scratch = screened_rates(k_Ne18_to_He4_O14_derived);
    jac.set(O14, Ne18, scratch);

    jac.set(O14, Ne19, 0.0);

    jac.set(O14, Ne20, 0.0);

    jac.set(O14, Ne21, 0.0);

    jac.set(O14, Na22, 0.0);

    jac.set(O14, Na23, 0.0);

    jac.set(O14, Mg22, 0.0);

    jac.set(O14, Mg24, 0.0);

    jac.set(O14, Al27, 0.0);

    jac.set(O14, Si28, 0.0);

    jac.set(O14, P31, 0.0);

    jac.set(O14, S32, 0.0);

    jac.set(O14, Ar36, 0.0);

    jac.set(O14, Ca40, 0.0);

    jac.set(O14, Ti44, 0.0);

    jac.set(O14, Cr48, 0.0);

    jac.set(O14, Fe52, 0.0);

    jac.set(O14, Ni56, 0.0);

    scratch = screened_rates(k_p_F18_to_He4_O15_reaclib)*Y(F18)*state.rho + screened_rates(k_p_N14_to_O15_reaclib)*Y(N14)*state.rho;
    jac.set(O15, H1, scratch);

    scratch = -screened_rates(k_He4_O15_to_Ne19_reaclib)*Y(O15)*state.rho - screened_rates(k_He4_O15_to_p_F18_derived)*Y(O15)*state.rho;
    jac.set(O15, He4, scratch);

    jac.set(O15, C12, 0.0);

    jac.set(O15, C13, 0.0);

    jac.set(O15, N13, 0.0);

    scratch = screened_rates(k_p_N14_to_O15_reaclib)*Y(H1)*state.rho;
    jac.set(O15, N14, scratch);

    jac.set(O15, N15, 0.0);

    jac.set(O15, O14, 0.0);

    scratch = -screened_rates(k_He4_O15_to_Ne19_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_O15_to_p_F18_derived)*Y(He4)*state.rho - screened_rates(k_O15_to_N15_reaclib) - screened_rates(k_O15_to_p_N14_derived);
    jac.set(O15, O15, scratch);

    jac.set(O15, O16, 0.0);

    jac.set(O15, O17, 0.0);

    jac.set(O15, O18, 0.0);

    jac.set(O15, F17, 0.0);

    scratch = screened_rates(k_p_F18_to_He4_O15_reaclib)*Y(H1)*state.rho;
    jac.set(O15, F18, scratch);

    jac.set(O15, F19, 0.0);

    jac.set(O15, Ne18, 0.0);

    scratch = screened_rates(k_Ne19_to_He4_O15_derived);
    jac.set(O15, Ne19, scratch);

    jac.set(O15, Ne20, 0.0);

    jac.set(O15, Ne21, 0.0);

    jac.set(O15, Na22, 0.0);

    jac.set(O15, Na23, 0.0);

    jac.set(O15, Mg22, 0.0);

    jac.set(O15, Mg24, 0.0);

    jac.set(O15, Al27, 0.0);

    jac.set(O15, Si28, 0.0);

    jac.set(O15, P31, 0.0);

    jac.set(O15, S32, 0.0);

    jac.set(O15, Ar36, 0.0);

    jac.set(O15, Ca40, 0.0);

    jac.set(O15, Ti44, 0.0);

    jac.set(O15, Cr48, 0.0);

    jac.set(O15, Fe52, 0.0);

    jac.set(O15, Ni56, 0.0);

    scratch = screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(F19)*state.rho + screened_rates(k_p_N15_to_O16_reaclib)*Y(N15)*state.rho - screened_rates(k_p_O16_to_F17_reaclib)*Y(O16)*state.rho - screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*state.rho;
    jac.set(O16, H1, scratch);

    scratch = screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*state.rho + screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(N13)*state.rho - screened_rates(k_He4_O16_to_Ne20_reaclib)*Y(O16)*state.rho - screened_rates(k_He4_O16_to_p_F19_derived)*Y(O16)*state.rho;
    jac.set(O16, He4, scratch);

    scratch = -screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(O16)*state.rho - screened_rates(k_C12_O16_to_Si28_modified)*Y(O16)*state.rho - screened_rates(k_C12_O16_to_p_Al27_reaclib)*Y(O16)*state.rho + screened_rates(k_He4_C12_to_O16_reaclib)*Y(He4)*state.rho;
    jac.set(O16, C12, scratch);

    jac.set(O16, C13, 0.0);

    scratch = screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(He4)*state.rho;
    jac.set(O16, N13, scratch);

    jac.set(O16, N14, 0.0);

    scratch = screened_rates(k_p_N15_to_O16_reaclib)*Y(H1)*state.rho;
    jac.set(O16, N15, scratch);

    jac.set(O16, O14, 0.0);

    jac.set(O16, O15, 0.0);

    scratch = -screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(C12)*state.rho - screened_rates(k_C12_O16_to_Si28_modified)*Y(C12)*state.rho - screened_rates(k_C12_O16_to_p_Al27_reaclib)*Y(C12)*state.rho - screened_rates(k_He4_O16_to_Ne20_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_O16_to_p_F19_derived)*Y(He4)*state.rho - 2.0*screened_rates(k_O16_O16_to_He4_Si28_reaclib)*Y(O16)*state.rho - 2.0*screened_rates(k_O16_O16_to_S32_modified)*Y(O16)*state.rho - 2.0*screened_rates(k_O16_O16_to_p_P31_reaclib)*Y(O16)*state.rho - screened_rates(k_O16_to_He4_C12_derived) - screened_rates(k_O16_to_p_N15_derived) - screened_rates(k_p_O16_to_F17_reaclib)*Y(H1)*state.rho - screened_rates(k_p_O16_to_He4_N13_derived)*Y(H1)*state.rho;
    jac.set(O16, O16, scratch);

    jac.set(O16, O17, 0.0);

    jac.set(O16, O18, 0.0);

    scratch = screened_rates(k_F17_to_p_O16_derived);
    jac.set(O16, F17, scratch);

    jac.set(O16, F18, 0.0);

    scratch = screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(H1)*state.rho;
    jac.set(O16, F19, scratch);

    jac.set(O16, Ne18, 0.0);

    jac.set(O16, Ne19, 0.0);

    scratch = screened_rates(k_Ne20_to_He4_O16_derived);
    jac.set(O16, Ne20, scratch);

    jac.set(O16, Ne21, 0.0);

    jac.set(O16, Na22, 0.0);

    jac.set(O16, Na23, 0.0);

    jac.set(O16, Mg22, 0.0);

    jac.set(O16, Mg24, 0.0);

    jac.set(O16, Al27, 0.0);

    jac.set(O16, Si28, 0.0);

    jac.set(O16, P31, 0.0);

    jac.set(O16, S32, 0.0);

    jac.set(O16, Ar36, 0.0);

    jac.set(O16, Ca40, 0.0);

    jac.set(O16, Ti44, 0.0);

    jac.set(O16, Cr48, 0.0);

    jac.set(O16, Fe52, 0.0);

    jac.set(O16, Ni56, 0.0);

    scratch = -screened_rates(k_p_O17_to_F18_reaclib)*Y(O17)*state.rho - screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(O17)*state.rho;
    jac.set(O17, H1, scratch);

    scratch = screened_rates(k_He4_N14_to_p_O17_derived)*Y(N14)*state.rho - screened_rates(k_He4_O17_to_Ne21_reaclib)*Y(O17)*state.rho;
    jac.set(O17, He4, scratch);

    jac.set(O17, C12, 0.0);

    jac.set(O17, C13, 0.0);

    jac.set(O17, N13, 0.0);

    scratch = screened_rates(k_He4_N14_to_p_O17_derived)*Y(He4)*state.rho;
    jac.set(O17, N14, scratch);

    jac.set(O17, N15, 0.0);

    jac.set(O17, O14, 0.0);

    jac.set(O17, O15, 0.0);

    jac.set(O17, O16, 0.0);

    scratch = -screened_rates(k_He4_O17_to_Ne21_reaclib)*Y(He4)*state.rho - screened_rates(k_O17_to_F17_weaktab) - screened_rates(k_p_O17_to_F18_reaclib)*Y(H1)*state.rho - screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(H1)*state.rho;
    jac.set(O17, O17, scratch);

    jac.set(O17, O18, 0.0);

    scratch = screened_rates(k_F17_to_O17_weaktab);
    jac.set(O17, F17, scratch);

    scratch = screened_rates(k_F18_to_p_O17_derived);
    jac.set(O17, F18, scratch);

    jac.set(O17, F19, 0.0);

    jac.set(O17, Ne18, 0.0);

    jac.set(O17, Ne19, 0.0);

    jac.set(O17, Ne20, 0.0);

    scratch = screened_rates(k_Ne21_to_He4_O17_derived);
    jac.set(O17, Ne21, scratch);

    jac.set(O17, Na22, 0.0);

    jac.set(O17, Na23, 0.0);

    jac.set(O17, Mg22, 0.0);

    jac.set(O17, Mg24, 0.0);

    jac.set(O17, Al27, 0.0);

    jac.set(O17, Si28, 0.0);

    jac.set(O17, P31, 0.0);

    jac.set(O17, S32, 0.0);

    jac.set(O17, Ar36, 0.0);

    jac.set(O17, Ca40, 0.0);

    jac.set(O17, Ti44, 0.0);

    jac.set(O17, Cr48, 0.0);

    jac.set(O17, Fe52, 0.0);

    jac.set(O17, Ni56, 0.0);

    scratch = -screened_rates(k_p_O18_to_F19_reaclib)*Y(O18)*state.rho - screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(O18)*state.rho;
    jac.set(O18, H1, scratch);

    scratch = screened_rates(k_He4_N15_to_p_O18_derived)*Y(N15)*state.rho;
    jac.set(O18, He4, scratch);

    jac.set(O18, C12, 0.0);

    jac.set(O18, C13, 0.0);

    jac.set(O18, N13, 0.0);

    jac.set(O18, N14, 0.0);

    scratch = screened_rates(k_He4_N15_to_p_O18_derived)*Y(He4)*state.rho;
    jac.set(O18, N15, scratch);

    jac.set(O18, O14, 0.0);

    jac.set(O18, O15, 0.0);

    jac.set(O18, O16, 0.0);

    jac.set(O18, O17, 0.0);

    scratch = -screened_rates(k_O18_to_F18_weaktab) - screened_rates(k_p_O18_to_F19_reaclib)*Y(H1)*state.rho - screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(H1)*state.rho;
    jac.set(O18, O18, scratch);

    jac.set(O18, F17, 0.0);

    scratch = screened_rates(k_F18_to_O18_weaktab);
    jac.set(O18, F18, scratch);

    scratch = screened_rates(k_F19_to_p_O18_derived);
    jac.set(O18, F19, scratch);

    jac.set(O18, Ne18, 0.0);

    jac.set(O18, Ne19, 0.0);

    jac.set(O18, Ne20, 0.0);

    jac.set(O18, Ne21, 0.0);

    jac.set(O18, Na22, 0.0);

    jac.set(O18, Na23, 0.0);

    jac.set(O18, Mg22, 0.0);

    jac.set(O18, Mg24, 0.0);

    jac.set(O18, Al27, 0.0);

    jac.set(O18, Si28, 0.0);

    jac.set(O18, P31, 0.0);

    jac.set(O18, S32, 0.0);

    jac.set(O18, Ar36, 0.0);

    jac.set(O18, Ca40, 0.0);

    jac.set(O18, Ti44, 0.0);

    jac.set(O18, Cr48, 0.0);

    jac.set(O18, Fe52, 0.0);

    jac.set(O18, Ni56, 0.0);

    scratch = -screened_rates(k_p_F17_to_He4_O14_derived)*Y(F17)*state.rho - screened_rates(k_p_F17_to_Ne18_reaclib)*Y(F17)*state.rho + screened_rates(k_p_Ne20_to_He4_F17_reaclib)*Y(Ne20)*state.rho + screened_rates(k_p_O16_to_F17_reaclib)*Y(O16)*state.rho;
    jac.set(F17, H1, scratch);

    scratch = -screened_rates(k_He4_F17_to_p_Ne20_derived)*Y(F17)*state.rho + screened_rates(k_He4_O14_to_p_F17_reaclib)*Y(O14)*state.rho;
    jac.set(F17, He4, scratch);

    jac.set(F17, C12, 0.0);

    jac.set(F17, C13, 0.0);

    jac.set(F17, N13, 0.0);

    jac.set(F17, N14, 0.0);

    jac.set(F17, N15, 0.0);

    scratch = screened_rates(k_He4_O14_to_p_F17_reaclib)*Y(He4)*state.rho;
    jac.set(F17, O14, scratch);

    jac.set(F17, O15, 0.0);

    scratch = screened_rates(k_p_O16_to_F17_reaclib)*Y(H1)*state.rho;
    jac.set(F17, O16, scratch);

    scratch = screened_rates(k_O17_to_F17_weaktab);
    jac.set(F17, O17, scratch);

    jac.set(F17, O18, 0.0);

    scratch = -screened_rates(k_F17_to_O17_weaktab) - screened_rates(k_F17_to_p_O16_derived) - screened_rates(k_He4_F17_to_p_Ne20_derived)*Y(He4)*state.rho - screened_rates(k_p_F17_to_He4_O14_derived)*Y(H1)*state.rho - screened_rates(k_p_F17_to_Ne18_reaclib)*Y(H1)*state.rho;
    jac.set(F17, F17, scratch);

    jac.set(F17, F18, 0.0);

    jac.set(F17, F19, 0.0);

    scratch = screened_rates(k_Ne18_to_p_F17_derived);
    jac.set(F17, Ne18, scratch);

    jac.set(F17, Ne19, 0.0);

    scratch = screened_rates(k_p_Ne20_to_He4_F17_reaclib)*Y(H1)*state.rho;
    jac.set(F17, Ne20, scratch);

    jac.set(F17, Ne21, 0.0);

    jac.set(F17, Na22, 0.0);

    jac.set(F17, Na23, 0.0);

    jac.set(F17, Mg22, 0.0);

    jac.set(F17, Mg24, 0.0);

    jac.set(F17, Al27, 0.0);

    jac.set(F17, Si28, 0.0);

    jac.set(F17, P31, 0.0);

    jac.set(F17, S32, 0.0);

    jac.set(F17, Ar36, 0.0);

    jac.set(F17, Ca40, 0.0);

    jac.set(F17, Ti44, 0.0);

    jac.set(F17, Cr48, 0.0);

    jac.set(F17, Fe52, 0.0);

    jac.set(F17, Ni56, 0.0);

    scratch = -screened_rates(k_p_F18_to_He4_O15_reaclib)*Y(F18)*state.rho - screened_rates(k_p_F18_to_Ne19_reaclib)*Y(F18)*state.rho + screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(Ne21)*state.rho + screened_rates(k_p_O17_to_F18_reaclib)*Y(O17)*state.rho;
    jac.set(F18, H1, scratch);

    scratch = -screened_rates(k_He4_F18_to_Na22_reaclib)*Y(F18)*state.rho - screened_rates(k_He4_F18_to_p_Ne21_reaclib)*Y(F18)*state.rho + screened_rates(k_He4_N14_to_F18_reaclib)*Y(N14)*state.rho + screened_rates(k_He4_O15_to_p_F18_derived)*Y(O15)*state.rho;
    jac.set(F18, He4, scratch);

    jac.set(F18, C12, 0.0);

    jac.set(F18, C13, 0.0);

    jac.set(F18, N13, 0.0);

    scratch = screened_rates(k_He4_N14_to_F18_reaclib)*Y(He4)*state.rho;
    jac.set(F18, N14, scratch);

    jac.set(F18, N15, 0.0);

    jac.set(F18, O14, 0.0);

    scratch = screened_rates(k_He4_O15_to_p_F18_derived)*Y(He4)*state.rho;
    jac.set(F18, O15, scratch);

    jac.set(F18, O16, 0.0);

    scratch = screened_rates(k_p_O17_to_F18_reaclib)*Y(H1)*state.rho;
    jac.set(F18, O17, scratch);

    scratch = screened_rates(k_O18_to_F18_weaktab);
    jac.set(F18, O18, scratch);

    jac.set(F18, F17, 0.0);

    scratch = -screened_rates(k_F18_to_He4_N14_derived) - screened_rates(k_F18_to_Ne18_weaktab) - screened_rates(k_F18_to_O18_weaktab) - screened_rates(k_F18_to_p_O17_derived) - screened_rates(k_He4_F18_to_Na22_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_F18_to_p_Ne21_reaclib)*Y(He4)*state.rho - screened_rates(k_p_F18_to_He4_O15_reaclib)*Y(H1)*state.rho - screened_rates(k_p_F18_to_Ne19_reaclib)*Y(H1)*state.rho;
    jac.set(F18, F18, scratch);

    jac.set(F18, F19, 0.0);

    scratch = screened_rates(k_Ne18_to_F18_weaktab);
    jac.set(F18, Ne18, scratch);

    scratch = screened_rates(k_Ne19_to_p_F18_derived);
    jac.set(F18, Ne19, scratch);

    jac.set(F18, Ne20, 0.0);

    scratch = screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(H1)*state.rho;
    jac.set(F18, Ne21, scratch);

    scratch = screened_rates(k_Na22_to_He4_F18_derived);
    jac.set(F18, Na22, scratch);

    jac.set(F18, Na23, 0.0);

    jac.set(F18, Mg22, 0.0);

    jac.set(F18, Mg24, 0.0);

    jac.set(F18, Al27, 0.0);

    jac.set(F18, Si28, 0.0);

    jac.set(F18, P31, 0.0);

    jac.set(F18, S32, 0.0);

    jac.set(F18, Ar36, 0.0);

    jac.set(F18, Ca40, 0.0);

    jac.set(F18, Ti44, 0.0);

    jac.set(F18, Cr48, 0.0);

    jac.set(F18, Fe52, 0.0);

    jac.set(F18, Ni56, 0.0);

    scratch = -screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(F19)*state.rho - screened_rates(k_p_F19_to_Ne20_reaclib)*Y(F19)*state.rho + screened_rates(k_p_O18_to_F19_reaclib)*Y(O18)*state.rho;
    jac.set(F19, H1, scratch);

    scratch = -screened_rates(k_He4_F19_to_Na23_reaclib)*Y(F19)*state.rho + screened_rates(k_He4_N15_to_F19_reaclib)*Y(N15)*state.rho + screened_rates(k_He4_O16_to_p_F19_derived)*Y(O16)*state.rho;
    jac.set(F19, He4, scratch);

    jac.set(F19, C12, 0.0);

    jac.set(F19, C13, 0.0);

    jac.set(F19, N13, 0.0);

    jac.set(F19, N14, 0.0);

    scratch = screened_rates(k_He4_N15_to_F19_reaclib)*Y(He4)*state.rho;
    jac.set(F19, N15, scratch);

    jac.set(F19, O14, 0.0);

    jac.set(F19, O15, 0.0);

    scratch = screened_rates(k_He4_O16_to_p_F19_derived)*Y(He4)*state.rho;
    jac.set(F19, O16, scratch);

    jac.set(F19, O17, 0.0);

    scratch = screened_rates(k_p_O18_to_F19_reaclib)*Y(H1)*state.rho;
    jac.set(F19, O18, scratch);

    jac.set(F19, F17, 0.0);

    jac.set(F19, F18, 0.0);

    scratch = -screened_rates(k_F19_to_He4_N15_derived) - screened_rates(k_F19_to_Ne19_weaktab) - screened_rates(k_F19_to_p_O18_derived) - screened_rates(k_He4_F19_to_Na23_reaclib)*Y(He4)*state.rho - screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(H1)*state.rho - screened_rates(k_p_F19_to_Ne20_reaclib)*Y(H1)*state.rho;
    jac.set(F19, F19, scratch);

    jac.set(F19, Ne18, 0.0);

    scratch = screened_rates(k_Ne19_to_F19_weaktab);
    jac.set(F19, Ne19, scratch);

    scratch = screened_rates(k_Ne20_to_p_F19_derived);
    jac.set(F19, Ne20, scratch);

    jac.set(F19, Ne21, 0.0);

    jac.set(F19, Na22, 0.0);

    scratch = screened_rates(k_Na23_to_He4_F19_derived);
    jac.set(F19, Na23, scratch);

    jac.set(F19, Mg22, 0.0);

    jac.set(F19, Mg24, 0.0);

    jac.set(F19, Al27, 0.0);

    jac.set(F19, Si28, 0.0);

    jac.set(F19, P31, 0.0);

    jac.set(F19, S32, 0.0);

    jac.set(F19, Ar36, 0.0);

    jac.set(F19, Ca40, 0.0);

    jac.set(F19, Ti44, 0.0);

    jac.set(F19, Cr48, 0.0);

    jac.set(F19, Fe52, 0.0);

    jac.set(F19, Ni56, 0.0);

    scratch = screened_rates(k_p_F17_to_Ne18_reaclib)*Y(F17)*state.rho;
    jac.set(Ne18, H1, scratch);

    scratch = -screened_rates(k_He4_Ne18_to_Mg22_reaclib)*Y(Ne18)*state.rho + screened_rates(k_He4_O14_to_Ne18_reaclib)*Y(O14)*state.rho;
    jac.set(Ne18, He4, scratch);

    jac.set(Ne18, C12, 0.0);

    jac.set(Ne18, C13, 0.0);

    jac.set(Ne18, N13, 0.0);

    jac.set(Ne18, N14, 0.0);

    jac.set(Ne18, N15, 0.0);

    scratch = screened_rates(k_He4_O14_to_Ne18_reaclib)*Y(He4)*state.rho;
    jac.set(Ne18, O14, scratch);

    jac.set(Ne18, O15, 0.0);

    jac.set(Ne18, O16, 0.0);

    jac.set(Ne18, O17, 0.0);

    jac.set(Ne18, O18, 0.0);

    scratch = screened_rates(k_p_F17_to_Ne18_reaclib)*Y(H1)*state.rho;
    jac.set(Ne18, F17, scratch);

    scratch = screened_rates(k_F18_to_Ne18_weaktab);
    jac.set(Ne18, F18, scratch);

    jac.set(Ne18, F19, 0.0);

    scratch = -screened_rates(k_He4_Ne18_to_Mg22_reaclib)*Y(He4)*state.rho - screened_rates(k_Ne18_to_F18_weaktab) - screened_rates(k_Ne18_to_He4_O14_derived) - screened_rates(k_Ne18_to_p_F17_derived);
    jac.set(Ne18, Ne18, scratch);

    jac.set(Ne18, Ne19, 0.0);

    jac.set(Ne18, Ne20, 0.0);

    jac.set(Ne18, Ne21, 0.0);

    jac.set(Ne18, Na22, 0.0);

    jac.set(Ne18, Na23, 0.0);

    scratch = screened_rates(k_Mg22_to_He4_Ne18_derived);
    jac.set(Ne18, Mg22, scratch);

    jac.set(Ne18, Mg24, 0.0);

    jac.set(Ne18, Al27, 0.0);

    jac.set(Ne18, Si28, 0.0);

    jac.set(Ne18, P31, 0.0);

    jac.set(Ne18, S32, 0.0);

    jac.set(Ne18, Ar36, 0.0);

    jac.set(Ne18, Ca40, 0.0);

    jac.set(Ne18, Ti44, 0.0);

    jac.set(Ne18, Cr48, 0.0);

    jac.set(Ne18, Fe52, 0.0);

    jac.set(Ne18, Ni56, 0.0);

    scratch = screened_rates(k_p_F18_to_Ne19_reaclib)*Y(F18)*state.rho + screened_rates(k_p_Na22_to_He4_Ne19_derived)*Y(Na22)*state.rho;
    jac.set(Ne19, H1, scratch);

    scratch = -screened_rates(k_He4_Ne19_to_p_Na22_reaclib)*Y(Ne19)*state.rho + screened_rates(k_He4_O15_to_Ne19_reaclib)*Y(O15)*state.rho;
    jac.set(Ne19, He4, scratch);

    jac.set(Ne19, C12, 0.0);

    jac.set(Ne19, C13, 0.0);

    jac.set(Ne19, N13, 0.0);

    jac.set(Ne19, N14, 0.0);

    jac.set(Ne19, N15, 0.0);

    jac.set(Ne19, O14, 0.0);

    scratch = screened_rates(k_He4_O15_to_Ne19_reaclib)*Y(He4)*state.rho;
    jac.set(Ne19, O15, scratch);

    jac.set(Ne19, O16, 0.0);

    jac.set(Ne19, O17, 0.0);

    jac.set(Ne19, O18, 0.0);

    jac.set(Ne19, F17, 0.0);

    scratch = screened_rates(k_p_F18_to_Ne19_reaclib)*Y(H1)*state.rho;
    jac.set(Ne19, F18, scratch);

    scratch = screened_rates(k_F19_to_Ne19_weaktab);
    jac.set(Ne19, F19, scratch);

    jac.set(Ne19, Ne18, 0.0);

    scratch = -screened_rates(k_He4_Ne19_to_p_Na22_reaclib)*Y(He4)*state.rho - screened_rates(k_Ne19_to_F19_weaktab) - screened_rates(k_Ne19_to_He4_O15_derived) - screened_rates(k_Ne19_to_p_F18_derived);
    jac.set(Ne19, Ne19, scratch);

    jac.set(Ne19, Ne20, 0.0);

    jac.set(Ne19, Ne21, 0.0);

    scratch = screened_rates(k_p_Na22_to_He4_Ne19_derived)*Y(H1)*state.rho;
    jac.set(Ne19, Na22, scratch);

    jac.set(Ne19, Na23, 0.0);

    jac.set(Ne19, Mg22, 0.0);

    jac.set(Ne19, Mg24, 0.0);

    jac.set(Ne19, Al27, 0.0);

    jac.set(Ne19, Si28, 0.0);

    jac.set(Ne19, P31, 0.0);

    jac.set(Ne19, S32, 0.0);

    jac.set(Ne19, Ar36, 0.0);

    jac.set(Ne19, Ca40, 0.0);

    jac.set(Ne19, Ti44, 0.0);

    jac.set(Ne19, Cr48, 0.0);

    jac.set(Ne19, Fe52, 0.0);

    jac.set(Ne19, Ni56, 0.0);

    scratch = screened_rates(k_p_F19_to_Ne20_reaclib)*Y(F19)*state.rho + screened_rates(k_p_Na23_to_He4_Ne20_reaclib)*Y(Na23)*state.rho - screened_rates(k_p_Ne20_to_He4_F17_reaclib)*Y(Ne20)*state.rho;
    jac.set(Ne20, H1, scratch);

    scratch = screened_rates(k_He4_F17_to_p_Ne20_derived)*Y(F17)*state.rho - screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(Ne20)*state.rho - screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(Ne20)*state.rho + screened_rates(k_He4_O16_to_Ne20_reaclib)*Y(O16)*state.rho;
    jac.set(Ne20, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*Y(C12)*state.rho;
    jac.set(Ne20, C12, scratch);

    jac.set(Ne20, C13, 0.0);

    jac.set(Ne20, N13, 0.0);

    jac.set(Ne20, N14, 0.0);

    jac.set(Ne20, N15, 0.0);

    jac.set(Ne20, O14, 0.0);

    jac.set(Ne20, O15, 0.0);

    scratch = screened_rates(k_He4_O16_to_Ne20_reaclib)*Y(He4)*state.rho;
    jac.set(Ne20, O16, scratch);

    jac.set(Ne20, O17, 0.0);

    jac.set(Ne20, O18, 0.0);

    scratch = screened_rates(k_He4_F17_to_p_Ne20_derived)*Y(He4)*state.rho;
    jac.set(Ne20, F17, scratch);

    jac.set(Ne20, F18, 0.0);

    scratch = screened_rates(k_p_F19_to_Ne20_reaclib)*Y(H1)*state.rho;
    jac.set(Ne20, F19, scratch);

    jac.set(Ne20, Ne18, 0.0);

    jac.set(Ne20, Ne19, 0.0);

    scratch = -screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(He4)*state.rho - screened_rates(k_Ne20_to_He4_O16_derived) - screened_rates(k_Ne20_to_p_F19_derived) - screened_rates(k_p_Ne20_to_He4_F17_reaclib)*Y(H1)*state.rho;
    jac.set(Ne20, Ne20, scratch);

    jac.set(Ne20, Ne21, 0.0);

    jac.set(Ne20, Na22, 0.0);

    scratch = screened_rates(k_p_Na23_to_He4_Ne20_reaclib)*Y(H1)*state.rho;
    jac.set(Ne20, Na23, scratch);

    jac.set(Ne20, Mg22, 0.0);

    scratch = screened_rates(k_Mg24_to_He4_Ne20_derived);
    jac.set(Ne20, Mg24, scratch);

    jac.set(Ne20, Al27, 0.0);

    jac.set(Ne20, Si28, 0.0);

    jac.set(Ne20, P31, 0.0);

    jac.set(Ne20, S32, 0.0);

    jac.set(Ne20, Ar36, 0.0);

    jac.set(Ne20, Ca40, 0.0);

    jac.set(Ne20, Ti44, 0.0);

    jac.set(Ne20, Cr48, 0.0);

    jac.set(Ne20, Fe52, 0.0);

    jac.set(Ne20, Ni56, 0.0);

    scratch = -screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(Ne21)*state.rho - screened_rates(k_p_Ne21_to_Na22_reaclib)*Y(Ne21)*state.rho;
    jac.set(Ne21, H1, scratch);

    scratch = screened_rates(k_He4_F18_to_p_Ne21_reaclib)*Y(F18)*state.rho + screened_rates(k_He4_O17_to_Ne21_reaclib)*Y(O17)*state.rho;
    jac.set(Ne21, He4, scratch);

    jac.set(Ne21, C12, 0.0);

    jac.set(Ne21, C13, 0.0);

    jac.set(Ne21, N13, 0.0);

    jac.set(Ne21, N14, 0.0);

    jac.set(Ne21, N15, 0.0);

    jac.set(Ne21, O14, 0.0);

    jac.set(Ne21, O15, 0.0);

    jac.set(Ne21, O16, 0.0);

    scratch = screened_rates(k_He4_O17_to_Ne21_reaclib)*Y(He4)*state.rho;
    jac.set(Ne21, O17, scratch);

    jac.set(Ne21, O18, 0.0);

    jac.set(Ne21, F17, 0.0);

    scratch = screened_rates(k_He4_F18_to_p_Ne21_reaclib)*Y(He4)*state.rho;
    jac.set(Ne21, F18, scratch);

    jac.set(Ne21, F19, 0.0);

    jac.set(Ne21, Ne18, 0.0);

    jac.set(Ne21, Ne19, 0.0);

    jac.set(Ne21, Ne20, 0.0);

    scratch = -screened_rates(k_Ne21_to_He4_O17_derived) - screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(H1)*state.rho - screened_rates(k_p_Ne21_to_Na22_reaclib)*Y(H1)*state.rho;
    jac.set(Ne21, Ne21, scratch);

    scratch = screened_rates(k_Na22_to_p_Ne21_derived);
    jac.set(Ne21, Na22, scratch);

    jac.set(Ne21, Na23, 0.0);

    jac.set(Ne21, Mg22, 0.0);

    jac.set(Ne21, Mg24, 0.0);

    jac.set(Ne21, Al27, 0.0);

    jac.set(Ne21, Si28, 0.0);

    jac.set(Ne21, P31, 0.0);

    jac.set(Ne21, S32, 0.0);

    jac.set(Ne21, Ar36, 0.0);

    jac.set(Ne21, Ca40, 0.0);

    jac.set(Ne21, Ti44, 0.0);

    jac.set(Ne21, Cr48, 0.0);

    jac.set(Ne21, Fe52, 0.0);

    jac.set(Ne21, Ni56, 0.0);

    scratch = -screened_rates(k_p_Na22_to_He4_Ne19_derived)*Y(Na22)*state.rho + screened_rates(k_p_Ne21_to_Na22_reaclib)*Y(Ne21)*state.rho;
    jac.set(Na22, H1, scratch);

    scratch = screened_rates(k_He4_F18_to_Na22_reaclib)*Y(F18)*state.rho + screened_rates(k_He4_Ne19_to_p_Na22_reaclib)*Y(Ne19)*state.rho;
    jac.set(Na22, He4, scratch);

    jac.set(Na22, C12, 0.0);

    jac.set(Na22, C13, 0.0);

    jac.set(Na22, N13, 0.0);

    jac.set(Na22, N14, 0.0);

    jac.set(Na22, N15, 0.0);

    jac.set(Na22, O14, 0.0);

    jac.set(Na22, O15, 0.0);

    jac.set(Na22, O16, 0.0);

    jac.set(Na22, O17, 0.0);

    jac.set(Na22, O18, 0.0);

    jac.set(Na22, F17, 0.0);

    scratch = screened_rates(k_He4_F18_to_Na22_reaclib)*Y(He4)*state.rho;
    jac.set(Na22, F18, scratch);

    jac.set(Na22, F19, 0.0);

    jac.set(Na22, Ne18, 0.0);

    scratch = screened_rates(k_He4_Ne19_to_p_Na22_reaclib)*Y(He4)*state.rho;
    jac.set(Na22, Ne19, scratch);

    jac.set(Na22, Ne20, 0.0);

    scratch = screened_rates(k_p_Ne21_to_Na22_reaclib)*Y(H1)*state.rho;
    jac.set(Na22, Ne21, scratch);

    scratch = -screened_rates(k_Na22_to_He4_F18_derived) - screened_rates(k_Na22_to_Mg22_weaktab) - screened_rates(k_Na22_to_p_Ne21_derived) - screened_rates(k_p_Na22_to_He4_Ne19_derived)*Y(H1)*state.rho;
    jac.set(Na22, Na22, scratch);

    jac.set(Na22, Na23, 0.0);

    scratch = screened_rates(k_Mg22_to_Na22_weaktab);
    jac.set(Na22, Mg22, scratch);

    jac.set(Na22, Mg24, 0.0);

    jac.set(Na22, Al27, 0.0);

    jac.set(Na22, Si28, 0.0);

    jac.set(Na22, P31, 0.0);

    jac.set(Na22, S32, 0.0);

    jac.set(Na22, Ar36, 0.0);

    jac.set(Na22, Ca40, 0.0);

    jac.set(Na22, Ti44, 0.0);

    jac.set(Na22, Cr48, 0.0);

    jac.set(Na22, Fe52, 0.0);

    jac.set(Na22, Ni56, 0.0);

    scratch = -screened_rates(k_p_Na23_to_He4_Ne20_reaclib)*Y(Na23)*state.rho - screened_rates(k_p_Na23_to_Mg24_reaclib)*Y(Na23)*state.rho;
    jac.set(Na23, H1, scratch);

    scratch = screened_rates(k_He4_F19_to_Na23_reaclib)*Y(F19)*state.rho + screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(Ne20)*state.rho;
    jac.set(Na23, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_p_Na23_reaclib)*Y(C12)*state.rho;
    jac.set(Na23, C12, scratch);

    jac.set(Na23, C13, 0.0);

    jac.set(Na23, N13, 0.0);

    jac.set(Na23, N14, 0.0);

    jac.set(Na23, N15, 0.0);

    jac.set(Na23, O14, 0.0);

    jac.set(Na23, O15, 0.0);

    jac.set(Na23, O16, 0.0);

    jac.set(Na23, O17, 0.0);

    jac.set(Na23, O18, 0.0);

    jac.set(Na23, F17, 0.0);

    jac.set(Na23, F18, 0.0);

    scratch = screened_rates(k_He4_F19_to_Na23_reaclib)*Y(He4)*state.rho;
    jac.set(Na23, F19, scratch);

    jac.set(Na23, Ne18, 0.0);

    jac.set(Na23, Ne19, 0.0);

    scratch = screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(He4)*state.rho;
    jac.set(Na23, Ne20, scratch);

    jac.set(Na23, Ne21, 0.0);

    jac.set(Na23, Na22, 0.0);

    scratch = -screened_rates(k_Na23_to_He4_F19_derived) - screened_rates(k_p_Na23_to_He4_Ne20_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Na23_to_Mg24_reaclib)*Y(H1)*state.rho;
    jac.set(Na23, Na23, scratch);

    jac.set(Na23, Mg22, 0.0);

    scratch = screened_rates(k_Mg24_to_p_Na23_derived);
    jac.set(Na23, Mg24, scratch);

    jac.set(Na23, Al27, 0.0);

    jac.set(Na23, Si28, 0.0);

    jac.set(Na23, P31, 0.0);

    jac.set(Na23, S32, 0.0);

    jac.set(Na23, Ar36, 0.0);

    jac.set(Na23, Ca40, 0.0);

    jac.set(Na23, Ti44, 0.0);

    jac.set(Na23, Cr48, 0.0);

    jac.set(Na23, Fe52, 0.0);

    jac.set(Na23, Ni56, 0.0);

    jac.set(Mg22, H1, 0.0);

    scratch = screened_rates(k_He4_Ne18_to_Mg22_reaclib)*Y(Ne18)*state.rho;
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

    scratch = screened_rates(k_He4_Ne18_to_Mg22_reaclib)*Y(He4)*state.rho;
    jac.set(Mg22, Ne18, scratch);

    jac.set(Mg22, Ne19, 0.0);

    jac.set(Mg22, Ne20, 0.0);

    jac.set(Mg22, Ne21, 0.0);

    scratch = screened_rates(k_Na22_to_Mg22_weaktab);
    jac.set(Mg22, Na22, scratch);

    jac.set(Mg22, Na23, 0.0);

    scratch = -screened_rates(k_Mg22_to_He4_Ne18_derived) - screened_rates(k_Mg22_to_Na22_weaktab);
    jac.set(Mg22, Mg22, scratch);

    jac.set(Mg22, Mg24, 0.0);

    jac.set(Mg22, Al27, 0.0);

    jac.set(Mg22, Si28, 0.0);

    jac.set(Mg22, P31, 0.0);

    jac.set(Mg22, S32, 0.0);

    jac.set(Mg22, Ar36, 0.0);

    jac.set(Mg22, Ca40, 0.0);

    jac.set(Mg22, Ti44, 0.0);

    jac.set(Mg22, Cr48, 0.0);

    jac.set(Mg22, Fe52, 0.0);

    jac.set(Mg22, Ni56, 0.0);

    scratch = screened_rates(k_p_Al27_to_He4_Mg24_reaclib)*Y(Al27)*state.rho + screened_rates(k_p_Na23_to_Mg24_reaclib)*Y(Na23)*state.rho;
    jac.set(Mg24, H1, scratch);

    scratch = -screened_rates(k_He4_Mg24_to_Si28_reaclib)*Y(Mg24)*state.rho - screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(Mg24)*state.rho + screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(Ne20)*state.rho;
    jac.set(Mg24, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_Mg24_modified)*Y(C12)*state.rho + screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(O16)*state.rho;
    jac.set(Mg24, C12, scratch);

    jac.set(Mg24, C13, 0.0);

    jac.set(Mg24, N13, 0.0);

    jac.set(Mg24, N14, 0.0);

    jac.set(Mg24, N15, 0.0);

    jac.set(Mg24, O14, 0.0);

    jac.set(Mg24, O15, 0.0);

    scratch = screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(C12)*state.rho;
    jac.set(Mg24, O16, scratch);

    jac.set(Mg24, O17, 0.0);

    jac.set(Mg24, O18, 0.0);

    jac.set(Mg24, F17, 0.0);

    jac.set(Mg24, F18, 0.0);

    jac.set(Mg24, F19, 0.0);

    jac.set(Mg24, Ne18, 0.0);

    jac.set(Mg24, Ne19, 0.0);

    scratch = screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(He4)*state.rho;
    jac.set(Mg24, Ne20, scratch);

    jac.set(Mg24, Ne21, 0.0);

    jac.set(Mg24, Na22, 0.0);

    scratch = screened_rates(k_p_Na23_to_Mg24_reaclib)*Y(H1)*state.rho;
    jac.set(Mg24, Na23, scratch);

    jac.set(Mg24, Mg22, 0.0);

    scratch = -screened_rates(k_He4_Mg24_to_Si28_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(He4)*state.rho - screened_rates(k_Mg24_to_He4_Ne20_derived) - screened_rates(k_Mg24_to_p_Na23_derived);
    jac.set(Mg24, Mg24, scratch);

    scratch = screened_rates(k_p_Al27_to_He4_Mg24_reaclib)*Y(H1)*state.rho;
    jac.set(Mg24, Al27, scratch);

    scratch = screened_rates(k_Si28_to_He4_Mg24_derived);
    jac.set(Mg24, Si28, scratch);

    jac.set(Mg24, P31, 0.0);

    jac.set(Mg24, S32, 0.0);

    jac.set(Mg24, Ar36, 0.0);

    jac.set(Mg24, Ca40, 0.0);

    jac.set(Mg24, Ti44, 0.0);

    jac.set(Mg24, Cr48, 0.0);

    jac.set(Mg24, Fe52, 0.0);

    jac.set(Mg24, Ni56, 0.0);

    scratch = -screened_rates(k_p_Al27_to_He4_Mg24_reaclib)*Y(Al27)*state.rho - screened_rates(k_p_Al27_to_Si28_reaclib)*Y(Al27)*state.rho;
    jac.set(Al27, H1, scratch);

    scratch = screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(Mg24)*state.rho;
    jac.set(Al27, He4, scratch);

    scratch = screened_rates(k_C12_O16_to_p_Al27_reaclib)*Y(O16)*state.rho;
    jac.set(Al27, C12, scratch);

    jac.set(Al27, C13, 0.0);

    jac.set(Al27, N13, 0.0);

    jac.set(Al27, N14, 0.0);

    jac.set(Al27, N15, 0.0);

    jac.set(Al27, O14, 0.0);

    jac.set(Al27, O15, 0.0);

    scratch = screened_rates(k_C12_O16_to_p_Al27_reaclib)*Y(C12)*state.rho;
    jac.set(Al27, O16, scratch);

    jac.set(Al27, O17, 0.0);

    jac.set(Al27, O18, 0.0);

    jac.set(Al27, F17, 0.0);

    jac.set(Al27, F18, 0.0);

    jac.set(Al27, F19, 0.0);

    jac.set(Al27, Ne18, 0.0);

    jac.set(Al27, Ne19, 0.0);

    jac.set(Al27, Ne20, 0.0);

    jac.set(Al27, Ne21, 0.0);

    jac.set(Al27, Na22, 0.0);

    jac.set(Al27, Na23, 0.0);

    jac.set(Al27, Mg22, 0.0);

    scratch = screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(He4)*state.rho;
    jac.set(Al27, Mg24, scratch);

    scratch = -screened_rates(k_p_Al27_to_He4_Mg24_reaclib)*Y(H1)*state.rho - screened_rates(k_p_Al27_to_Si28_reaclib)*Y(H1)*state.rho;
    jac.set(Al27, Al27, scratch);

    scratch = screened_rates(k_Si28_to_p_Al27_derived);
    jac.set(Al27, Si28, scratch);

    jac.set(Al27, P31, 0.0);

    jac.set(Al27, S32, 0.0);

    jac.set(Al27, Ar36, 0.0);

    jac.set(Al27, Ca40, 0.0);

    jac.set(Al27, Ti44, 0.0);

    jac.set(Al27, Cr48, 0.0);

    jac.set(Al27, Fe52, 0.0);

    jac.set(Al27, Ni56, 0.0);

    scratch = screened_rates(k_p_Al27_to_Si28_reaclib)*Y(Al27)*state.rho + screened_rates(k_p_P31_to_He4_Si28_reaclib)*Y(P31)*state.rho;
    jac.set(Si28, H1, scratch);

    scratch = screened_rates(k_He4_Mg24_to_Si28_reaclib)*Y(Mg24)*state.rho - screened_rates(k_He4_Si28_to_S32_reaclib)*Y(Si28)*state.rho - screened_rates(k_He4_Si28_to_p_P31_derived)*Y(Si28)*state.rho;
    jac.set(Si28, He4, scratch);

    scratch = screened_rates(k_C12_O16_to_Si28_modified)*Y(O16)*state.rho;
    jac.set(Si28, C12, scratch);

    jac.set(Si28, C13, 0.0);

    jac.set(Si28, N13, 0.0);

    jac.set(Si28, N14, 0.0);

    jac.set(Si28, N15, 0.0);

    jac.set(Si28, O14, 0.0);

    jac.set(Si28, O15, 0.0);

    scratch = screened_rates(k_C12_O16_to_Si28_modified)*Y(C12)*state.rho + 1.0*screened_rates(k_O16_O16_to_He4_Si28_reaclib)*Y(O16)*state.rho;
    jac.set(Si28, O16, scratch);

    jac.set(Si28, O17, 0.0);

    jac.set(Si28, O18, 0.0);

    jac.set(Si28, F17, 0.0);

    jac.set(Si28, F18, 0.0);

    jac.set(Si28, F19, 0.0);

    jac.set(Si28, Ne18, 0.0);

    jac.set(Si28, Ne19, 0.0);

    jac.set(Si28, Ne20, 0.0);

    jac.set(Si28, Ne21, 0.0);

    jac.set(Si28, Na22, 0.0);

    jac.set(Si28, Na23, 0.0);

    jac.set(Si28, Mg22, 0.0);

    scratch = screened_rates(k_He4_Mg24_to_Si28_reaclib)*Y(He4)*state.rho;
    jac.set(Si28, Mg24, scratch);

    scratch = screened_rates(k_p_Al27_to_Si28_reaclib)*Y(H1)*state.rho;
    jac.set(Si28, Al27, scratch);

    scratch = -screened_rates(k_He4_Si28_to_S32_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_Si28_to_p_P31_derived)*Y(He4)*state.rho - screened_rates(k_Si28_to_He4_Mg24_derived) - screened_rates(k_Si28_to_p_Al27_derived);
    jac.set(Si28, Si28, scratch);

    scratch = screened_rates(k_p_P31_to_He4_Si28_reaclib)*Y(H1)*state.rho;
    jac.set(Si28, P31, scratch);

    scratch = screened_rates(k_S32_to_He4_Si28_derived);
    jac.set(Si28, S32, scratch);

    jac.set(Si28, Ar36, 0.0);

    jac.set(Si28, Ca40, 0.0);

    jac.set(Si28, Ti44, 0.0);

    jac.set(Si28, Cr48, 0.0);

    jac.set(Si28, Fe52, 0.0);

    jac.set(Si28, Ni56, 0.0);

    scratch = -screened_rates(k_p_P31_to_He4_Si28_reaclib)*Y(P31)*state.rho - screened_rates(k_p_P31_to_S32_reaclib)*Y(P31)*state.rho;
    jac.set(P31, H1, scratch);

    scratch = screened_rates(k_He4_Si28_to_p_P31_derived)*Y(Si28)*state.rho;
    jac.set(P31, He4, scratch);

    jac.set(P31, C12, 0.0);

    jac.set(P31, C13, 0.0);

    jac.set(P31, N13, 0.0);

    jac.set(P31, N14, 0.0);

    jac.set(P31, N15, 0.0);

    jac.set(P31, O14, 0.0);

    jac.set(P31, O15, 0.0);

    scratch = 1.0*screened_rates(k_O16_O16_to_p_P31_reaclib)*Y(O16)*state.rho;
    jac.set(P31, O16, scratch);

    jac.set(P31, O17, 0.0);

    jac.set(P31, O18, 0.0);

    jac.set(P31, F17, 0.0);

    jac.set(P31, F18, 0.0);

    jac.set(P31, F19, 0.0);

    jac.set(P31, Ne18, 0.0);

    jac.set(P31, Ne19, 0.0);

    jac.set(P31, Ne20, 0.0);

    jac.set(P31, Ne21, 0.0);

    jac.set(P31, Na22, 0.0);

    jac.set(P31, Na23, 0.0);

    jac.set(P31, Mg22, 0.0);

    jac.set(P31, Mg24, 0.0);

    jac.set(P31, Al27, 0.0);

    scratch = screened_rates(k_He4_Si28_to_p_P31_derived)*Y(He4)*state.rho;
    jac.set(P31, Si28, scratch);

    scratch = -screened_rates(k_p_P31_to_He4_Si28_reaclib)*Y(H1)*state.rho - screened_rates(k_p_P31_to_S32_reaclib)*Y(H1)*state.rho;
    jac.set(P31, P31, scratch);

    scratch = screened_rates(k_S32_to_p_P31_derived);
    jac.set(P31, S32, scratch);

    jac.set(P31, Ar36, 0.0);

    jac.set(P31, Ca40, 0.0);

    jac.set(P31, Ti44, 0.0);

    jac.set(P31, Cr48, 0.0);

    jac.set(P31, Fe52, 0.0);

    jac.set(P31, Ni56, 0.0);

    scratch = screened_rates(k_p_P31_to_S32_reaclib)*Y(P31)*state.rho;
    jac.set(S32, H1, scratch);

    scratch = screened_rates(k_He4_Si28_to_S32_reaclib)*Y(Si28)*state.rho - screened_rates(k_S32_He4_to_Ar36_approx)*Y(S32)*state.rho;
    jac.set(S32, He4, scratch);

    jac.set(S32, C12, 0.0);

    jac.set(S32, C13, 0.0);

    jac.set(S32, N13, 0.0);

    jac.set(S32, N14, 0.0);

    jac.set(S32, N15, 0.0);

    jac.set(S32, O14, 0.0);

    jac.set(S32, O15, 0.0);

    scratch = 1.0*screened_rates(k_O16_O16_to_S32_modified)*Y(O16)*state.rho;
    jac.set(S32, O16, scratch);

    jac.set(S32, O17, 0.0);

    jac.set(S32, O18, 0.0);

    jac.set(S32, F17, 0.0);

    jac.set(S32, F18, 0.0);

    jac.set(S32, F19, 0.0);

    jac.set(S32, Ne18, 0.0);

    jac.set(S32, Ne19, 0.0);

    jac.set(S32, Ne20, 0.0);

    jac.set(S32, Ne21, 0.0);

    jac.set(S32, Na22, 0.0);

    jac.set(S32, Na23, 0.0);

    jac.set(S32, Mg22, 0.0);

    jac.set(S32, Mg24, 0.0);

    jac.set(S32, Al27, 0.0);

    scratch = screened_rates(k_He4_Si28_to_S32_reaclib)*Y(He4)*state.rho;
    jac.set(S32, Si28, scratch);

    scratch = screened_rates(k_p_P31_to_S32_reaclib)*Y(H1)*state.rho;
    jac.set(S32, P31, scratch);

    scratch = -screened_rates(k_S32_He4_to_Ar36_approx)*Y(He4)*state.rho - screened_rates(k_S32_to_He4_Si28_derived) - screened_rates(k_S32_to_p_P31_derived);
    jac.set(S32, S32, scratch);

    scratch = screened_rates(k_Ar36_to_S32_He4_approx);
    jac.set(S32, Ar36, scratch);

    jac.set(S32, Ca40, 0.0);

    jac.set(S32, Ti44, 0.0);

    jac.set(S32, Cr48, 0.0);

    jac.set(S32, Fe52, 0.0);

    jac.set(S32, Ni56, 0.0);

    jac.set(Ar36, H1, 0.0);

    scratch = -screened_rates(k_Ar36_He4_to_Ca40_approx)*Y(Ar36)*state.rho + screened_rates(k_S32_He4_to_Ar36_approx)*Y(S32)*state.rho;
    jac.set(Ar36, He4, scratch);

    jac.set(Ar36, C12, 0.0);

    jac.set(Ar36, C13, 0.0);

    jac.set(Ar36, N13, 0.0);

    jac.set(Ar36, N14, 0.0);

    jac.set(Ar36, N15, 0.0);

    jac.set(Ar36, O14, 0.0);

    jac.set(Ar36, O15, 0.0);

    jac.set(Ar36, O16, 0.0);

    jac.set(Ar36, O17, 0.0);

    jac.set(Ar36, O18, 0.0);

    jac.set(Ar36, F17, 0.0);

    jac.set(Ar36, F18, 0.0);

    jac.set(Ar36, F19, 0.0);

    jac.set(Ar36, Ne18, 0.0);

    jac.set(Ar36, Ne19, 0.0);

    jac.set(Ar36, Ne20, 0.0);

    jac.set(Ar36, Ne21, 0.0);

    jac.set(Ar36, Na22, 0.0);

    jac.set(Ar36, Na23, 0.0);

    jac.set(Ar36, Mg22, 0.0);

    jac.set(Ar36, Mg24, 0.0);

    jac.set(Ar36, Al27, 0.0);

    jac.set(Ar36, Si28, 0.0);

    jac.set(Ar36, P31, 0.0);

    scratch = screened_rates(k_S32_He4_to_Ar36_approx)*Y(He4)*state.rho;
    jac.set(Ar36, S32, scratch);

    scratch = -screened_rates(k_Ar36_He4_to_Ca40_approx)*Y(He4)*state.rho - screened_rates(k_Ar36_to_S32_He4_approx);
    jac.set(Ar36, Ar36, scratch);

    scratch = screened_rates(k_Ca40_to_Ar36_He4_approx);
    jac.set(Ar36, Ca40, scratch);

    jac.set(Ar36, Ti44, 0.0);

    jac.set(Ar36, Cr48, 0.0);

    jac.set(Ar36, Fe52, 0.0);

    jac.set(Ar36, Ni56, 0.0);

    jac.set(Ca40, H1, 0.0);

    scratch = screened_rates(k_Ar36_He4_to_Ca40_approx)*Y(Ar36)*state.rho - screened_rates(k_Ca40_He4_to_Ti44_approx)*Y(Ca40)*state.rho;
    jac.set(Ca40, He4, scratch);

    jac.set(Ca40, C12, 0.0);

    jac.set(Ca40, C13, 0.0);

    jac.set(Ca40, N13, 0.0);

    jac.set(Ca40, N14, 0.0);

    jac.set(Ca40, N15, 0.0);

    jac.set(Ca40, O14, 0.0);

    jac.set(Ca40, O15, 0.0);

    jac.set(Ca40, O16, 0.0);

    jac.set(Ca40, O17, 0.0);

    jac.set(Ca40, O18, 0.0);

    jac.set(Ca40, F17, 0.0);

    jac.set(Ca40, F18, 0.0);

    jac.set(Ca40, F19, 0.0);

    jac.set(Ca40, Ne18, 0.0);

    jac.set(Ca40, Ne19, 0.0);

    jac.set(Ca40, Ne20, 0.0);

    jac.set(Ca40, Ne21, 0.0);

    jac.set(Ca40, Na22, 0.0);

    jac.set(Ca40, Na23, 0.0);

    jac.set(Ca40, Mg22, 0.0);

    jac.set(Ca40, Mg24, 0.0);

    jac.set(Ca40, Al27, 0.0);

    jac.set(Ca40, Si28, 0.0);

    jac.set(Ca40, P31, 0.0);

    jac.set(Ca40, S32, 0.0);

    scratch = screened_rates(k_Ar36_He4_to_Ca40_approx)*Y(He4)*state.rho;
    jac.set(Ca40, Ar36, scratch);

    scratch = -screened_rates(k_Ca40_He4_to_Ti44_approx)*Y(He4)*state.rho - screened_rates(k_Ca40_to_Ar36_He4_approx);
    jac.set(Ca40, Ca40, scratch);

    scratch = screened_rates(k_Ti44_to_Ca40_He4_approx);
    jac.set(Ca40, Ti44, scratch);

    jac.set(Ca40, Cr48, 0.0);

    jac.set(Ca40, Fe52, 0.0);

    jac.set(Ca40, Ni56, 0.0);

    jac.set(Ti44, H1, 0.0);

    scratch = screened_rates(k_Ca40_He4_to_Ti44_approx)*Y(Ca40)*state.rho - screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(Ti44)*state.rho;
    jac.set(Ti44, He4, scratch);

    jac.set(Ti44, C12, 0.0);

    jac.set(Ti44, C13, 0.0);

    jac.set(Ti44, N13, 0.0);

    jac.set(Ti44, N14, 0.0);

    jac.set(Ti44, N15, 0.0);

    jac.set(Ti44, O14, 0.0);

    jac.set(Ti44, O15, 0.0);

    jac.set(Ti44, O16, 0.0);

    jac.set(Ti44, O17, 0.0);

    jac.set(Ti44, O18, 0.0);

    jac.set(Ti44, F17, 0.0);

    jac.set(Ti44, F18, 0.0);

    jac.set(Ti44, F19, 0.0);

    jac.set(Ti44, Ne18, 0.0);

    jac.set(Ti44, Ne19, 0.0);

    jac.set(Ti44, Ne20, 0.0);

    jac.set(Ti44, Ne21, 0.0);

    jac.set(Ti44, Na22, 0.0);

    jac.set(Ti44, Na23, 0.0);

    jac.set(Ti44, Mg22, 0.0);

    jac.set(Ti44, Mg24, 0.0);

    jac.set(Ti44, Al27, 0.0);

    jac.set(Ti44, Si28, 0.0);

    jac.set(Ti44, P31, 0.0);

    jac.set(Ti44, S32, 0.0);

    jac.set(Ti44, Ar36, 0.0);

    scratch = screened_rates(k_Ca40_He4_to_Ti44_approx)*Y(He4)*state.rho;
    jac.set(Ti44, Ca40, scratch);

    scratch = -screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(He4)*state.rho - screened_rates(k_Ti44_to_Ca40_He4_approx);
    jac.set(Ti44, Ti44, scratch);

    scratch = screened_rates(k_Cr48_to_Ti44_He4_approx);
    jac.set(Ti44, Cr48, scratch);

    jac.set(Ti44, Fe52, 0.0);

    jac.set(Ti44, Ni56, 0.0);

    jac.set(Cr48, H1, 0.0);

    scratch = -screened_rates(k_Cr48_He4_to_Fe52_approx)*Y(Cr48)*state.rho + screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(Ti44)*state.rho;
    jac.set(Cr48, He4, scratch);

    jac.set(Cr48, C12, 0.0);

    jac.set(Cr48, C13, 0.0);

    jac.set(Cr48, N13, 0.0);

    jac.set(Cr48, N14, 0.0);

    jac.set(Cr48, N15, 0.0);

    jac.set(Cr48, O14, 0.0);

    jac.set(Cr48, O15, 0.0);

    jac.set(Cr48, O16, 0.0);

    jac.set(Cr48, O17, 0.0);

    jac.set(Cr48, O18, 0.0);

    jac.set(Cr48, F17, 0.0);

    jac.set(Cr48, F18, 0.0);

    jac.set(Cr48, F19, 0.0);

    jac.set(Cr48, Ne18, 0.0);

    jac.set(Cr48, Ne19, 0.0);

    jac.set(Cr48, Ne20, 0.0);

    jac.set(Cr48, Ne21, 0.0);

    jac.set(Cr48, Na22, 0.0);

    jac.set(Cr48, Na23, 0.0);

    jac.set(Cr48, Mg22, 0.0);

    jac.set(Cr48, Mg24, 0.0);

    jac.set(Cr48, Al27, 0.0);

    jac.set(Cr48, Si28, 0.0);

    jac.set(Cr48, P31, 0.0);

    jac.set(Cr48, S32, 0.0);

    jac.set(Cr48, Ar36, 0.0);

    jac.set(Cr48, Ca40, 0.0);

    scratch = screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(He4)*state.rho;
    jac.set(Cr48, Ti44, scratch);

    scratch = -screened_rates(k_Cr48_He4_to_Fe52_approx)*Y(He4)*state.rho - screened_rates(k_Cr48_to_Ti44_He4_approx);
    jac.set(Cr48, Cr48, scratch);

    scratch = screened_rates(k_Fe52_to_Cr48_He4_approx);
    jac.set(Cr48, Fe52, scratch);

    jac.set(Cr48, Ni56, 0.0);

    jac.set(Fe52, H1, 0.0);

    scratch = screened_rates(k_Cr48_He4_to_Fe52_approx)*Y(Cr48)*state.rho - screened_rates(k_Fe52_He4_to_Ni56_approx)*Y(Fe52)*state.rho;
    jac.set(Fe52, He4, scratch);

    jac.set(Fe52, C12, 0.0);

    jac.set(Fe52, C13, 0.0);

    jac.set(Fe52, N13, 0.0);

    jac.set(Fe52, N14, 0.0);

    jac.set(Fe52, N15, 0.0);

    jac.set(Fe52, O14, 0.0);

    jac.set(Fe52, O15, 0.0);

    jac.set(Fe52, O16, 0.0);

    jac.set(Fe52, O17, 0.0);

    jac.set(Fe52, O18, 0.0);

    jac.set(Fe52, F17, 0.0);

    jac.set(Fe52, F18, 0.0);

    jac.set(Fe52, F19, 0.0);

    jac.set(Fe52, Ne18, 0.0);

    jac.set(Fe52, Ne19, 0.0);

    jac.set(Fe52, Ne20, 0.0);

    jac.set(Fe52, Ne21, 0.0);

    jac.set(Fe52, Na22, 0.0);

    jac.set(Fe52, Na23, 0.0);

    jac.set(Fe52, Mg22, 0.0);

    jac.set(Fe52, Mg24, 0.0);

    jac.set(Fe52, Al27, 0.0);

    jac.set(Fe52, Si28, 0.0);

    jac.set(Fe52, P31, 0.0);

    jac.set(Fe52, S32, 0.0);

    jac.set(Fe52, Ar36, 0.0);

    jac.set(Fe52, Ca40, 0.0);

    jac.set(Fe52, Ti44, 0.0);

    scratch = screened_rates(k_Cr48_He4_to_Fe52_approx)*Y(He4)*state.rho;
    jac.set(Fe52, Cr48, scratch);

    scratch = -screened_rates(k_Fe52_He4_to_Ni56_approx)*Y(He4)*state.rho - screened_rates(k_Fe52_to_Cr48_He4_approx);
    jac.set(Fe52, Fe52, scratch);

    scratch = screened_rates(k_Ni56_to_Fe52_He4_approx);
    jac.set(Fe52, Ni56, scratch);

    jac.set(Ni56, H1, 0.0);

    scratch = screened_rates(k_Fe52_He4_to_Ni56_approx)*Y(Fe52)*state.rho;
    jac.set(Ni56, He4, scratch);

    jac.set(Ni56, C12, 0.0);

    jac.set(Ni56, C13, 0.0);

    jac.set(Ni56, N13, 0.0);

    jac.set(Ni56, N14, 0.0);

    jac.set(Ni56, N15, 0.0);

    jac.set(Ni56, O14, 0.0);

    jac.set(Ni56, O15, 0.0);

    jac.set(Ni56, O16, 0.0);

    jac.set(Ni56, O17, 0.0);

    jac.set(Ni56, O18, 0.0);

    jac.set(Ni56, F17, 0.0);

    jac.set(Ni56, F18, 0.0);

    jac.set(Ni56, F19, 0.0);

    jac.set(Ni56, Ne18, 0.0);

    jac.set(Ni56, Ne19, 0.0);

    jac.set(Ni56, Ne20, 0.0);

    jac.set(Ni56, Ne21, 0.0);

    jac.set(Ni56, Na22, 0.0);

    jac.set(Ni56, Na23, 0.0);

    jac.set(Ni56, Mg22, 0.0);

    jac.set(Ni56, Mg24, 0.0);

    jac.set(Ni56, Al27, 0.0);

    jac.set(Ni56, Si28, 0.0);

    jac.set(Ni56, P31, 0.0);

    jac.set(Ni56, S32, 0.0);

    jac.set(Ni56, Ar36, 0.0);

    jac.set(Ni56, Ca40, 0.0);

    jac.set(Ni56, Ti44, 0.0);

    jac.set(Ni56, Cr48, 0.0);

    scratch = screened_rates(k_Fe52_He4_to_Ni56_approx)*Y(He4)*state.rho;
    jac.set(Ni56, Fe52, scratch);

    scratch = -screened_rates(k_Ni56_to_Fe52_He4_approx);
    jac.set(Ni56, Ni56, scratch);


}



template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
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

    init_tabular();

}


#endif
