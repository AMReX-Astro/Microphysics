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

    ratraw = rate_eval.screened_rates(k_He4_Ne20_to_C12_C12_derived);
    rate_eval.screened_rates(k_He4_Ne20_to_C12_C12_derived) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ne20_to_C12_C12_derived);
        rate_eval.dscreened_rates_dT(k_He4_Ne20_to_C12_C12_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(6.0_rt, 12.0_rt, 6.0_rt, 12.0_rt);
        static_assert(scn_fac.z1 == 6.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_C12_C12_to_He4_Ne20_reaclib);
    rate_eval.screened_rates(k_C12_C12_to_He4_Ne20_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_C12_to_He4_Ne20_reaclib);
        rate_eval.dscreened_rates_dT(k_C12_C12_to_He4_Ne20_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
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

    ratraw = rate_eval.screened_rates(k_C12_O16_to_He4_Mg24_reaclib);
    rate_eval.screened_rates(k_C12_O16_to_He4_Mg24_reaclib) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_O16_to_He4_Mg24_reaclib);
        rate_eval.dscreened_rates_dT(k_C12_O16_to_He4_Mg24_reaclib) = ratraw * dscor_dt + dratraw_dT * scor;
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
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 12.0_rt, 24.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }

    ratraw = rate_eval.screened_rates(k_He4_Mg24_to_C12_O16_derived);
    rate_eval.screened_rates(k_He4_Mg24_to_C12_O16_derived) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mg24_to_C12_O16_derived);
        rate_eval.dscreened_rates_dT(k_He4_Mg24_to_C12_O16_derived) = ratraw * dscor_dt + dratraw_dT * scor;
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

    ydot_nuc(Mg22) = 0.0_rt;

    ydot_nuc(Mg24) = 0.0_rt;

    ydot_nuc(Fe56) = 0.0_rt;

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
        (screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(He4)*Y(N13)*state.rho + -screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*Y(H1)*state.rho) +
        (-screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(N15)*Y(H1)*state.rho + screened_rates(k_He4_C12_to_p_N15_derived)*Y(C12)*Y(He4)*state.rho) +
        (screened_rates(k_He4_O14_to_p_F17_reaclib)*Y(He4)*Y(O14)*state.rho + -screened_rates(k_p_F17_to_He4_O14_derived)*Y(F17)*Y(H1)*state.rho) +
        (-screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(O17)*Y(H1)*state.rho + screened_rates(k_He4_N14_to_p_O17_derived)*Y(He4)*Y(N14)*state.rho) +
        (-screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(O18)*Y(H1)*state.rho + screened_rates(k_He4_N15_to_p_O18_derived)*Y(He4)*Y(N15)*state.rho) +
        (-screened_rates(k_p_F18_to_He4_O15_reaclib)*Y(F18)*Y(H1)*state.rho + screened_rates(k_He4_O15_to_p_F18_derived)*Y(He4)*Y(O15)*state.rho) +
        (-screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(F19)*Y(H1)*state.rho + screened_rates(k_He4_O16_to_p_F19_derived)*Y(He4)*Y(O16)*state.rho) +
        (screened_rates(k_He4_F17_to_p_Ne20_derived)*Y(F17)*Y(He4)*state.rho + -screened_rates(k_p_Ne20_to_He4_F17_reaclib)*Y(Ne20)*Y(H1)*state.rho);

    ydot_nuc(He4) =
        (-screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*Y(He4)*state.rho + screened_rates(k_O16_to_He4_C12_derived)*Y(O16)) +
        (-screened_rates(k_He4_N14_to_F18_reaclib)*Y(He4)*Y(N14)*state.rho + screened_rates(k_F18_to_He4_N14_derived)*Y(F18)) +
        (-screened_rates(k_He4_N15_to_F19_reaclib)*Y(He4)*Y(N15)*state.rho + screened_rates(k_F19_to_He4_N15_derived)*Y(F19)) +
        (-screened_rates(k_He4_O14_to_Ne18_reaclib)*Y(He4)*Y(O14)*state.rho + screened_rates(k_Ne18_to_He4_O14_derived)*Y(Ne18)) +
        (-screened_rates(k_He4_O15_to_Ne19_reaclib)*Y(He4)*Y(O15)*state.rho + screened_rates(k_Ne19_to_He4_O15_derived)*Y(Ne19)) +
        (-screened_rates(k_He4_O16_to_Ne20_reaclib)*Y(He4)*Y(O16)*state.rho + screened_rates(k_Ne20_to_He4_O16_derived)*Y(Ne20)) +
        (-screened_rates(k_He4_Ne18_to_Mg22_reaclib)*Y(He4)*Y(Ne18)*state.rho + screened_rates(k_Mg22_to_He4_Ne18_derived)*Y(Mg22)) +
        (-screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(He4)*Y(Ne20)*state.rho + screened_rates(k_Mg24_to_He4_Ne20_derived)*Y(Mg24)) +
        (0.5*screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho + -screened_rates(k_He4_Ne20_to_C12_C12_derived)*Y(He4)*Y(Ne20)*state.rho) +
        (-screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(He4)*Y(N13)*state.rho + screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*Y(H1)*state.rho) +
        (screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(N15)*Y(H1)*state.rho + -screened_rates(k_He4_C12_to_p_N15_derived)*Y(C12)*Y(He4)*state.rho) +
        (-screened_rates(k_He4_O14_to_p_F17_reaclib)*Y(He4)*Y(O14)*state.rho + screened_rates(k_p_F17_to_He4_O14_derived)*Y(F17)*Y(H1)*state.rho) +
        (screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(C12)*Y(O16)*state.rho + -screened_rates(k_He4_Mg24_to_C12_O16_derived)*Y(He4)*Y(Mg24)*state.rho) +
        (screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(O17)*Y(H1)*state.rho + -screened_rates(k_He4_N14_to_p_O17_derived)*Y(He4)*Y(N14)*state.rho) +
        (screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(O18)*Y(H1)*state.rho + -screened_rates(k_He4_N15_to_p_O18_derived)*Y(He4)*Y(N15)*state.rho) +
        (screened_rates(k_p_F18_to_He4_O15_reaclib)*Y(F18)*Y(H1)*state.rho + -screened_rates(k_He4_O15_to_p_F18_derived)*Y(He4)*Y(O15)*state.rho) +
        (screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(F19)*Y(H1)*state.rho + -screened_rates(k_He4_O16_to_p_F19_derived)*Y(He4)*Y(O16)*state.rho) +
        (-0.5*screened_rates(k_He4_He4_He4_to_C12_reaclib)*amrex::Math::powi<3>(Y(He4))*amrex::Math::powi<2>(state.rho) + 3.0*screened_rates(k_C12_to_He4_He4_He4_derived)*Y(C12)) +
        (-screened_rates(k_He4_F17_to_p_Ne20_derived)*Y(F17)*Y(He4)*state.rho + screened_rates(k_p_Ne20_to_He4_F17_reaclib)*Y(Ne20)*Y(H1)*state.rho);

    ydot_nuc(C12) =
        (-screened_rates(k_p_C12_to_N13_reaclib)*Y(C12)*Y(H1)*state.rho + screened_rates(k_N13_to_p_C12_derived)*Y(N13)) +
        (-screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*Y(He4)*state.rho + screened_rates(k_O16_to_He4_C12_derived)*Y(O16)) +
        (-screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho + 2.0*screened_rates(k_He4_Ne20_to_C12_C12_derived)*Y(He4)*Y(Ne20)*state.rho) +
        (screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(N15)*Y(H1)*state.rho + -screened_rates(k_He4_C12_to_p_N15_derived)*Y(C12)*Y(He4)*state.rho) +
        (-screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(C12)*Y(O16)*state.rho + screened_rates(k_He4_Mg24_to_C12_O16_derived)*Y(He4)*Y(Mg24)*state.rho) +
        (0.16666666666666667*screened_rates(k_He4_He4_He4_to_C12_reaclib)*amrex::Math::powi<3>(Y(He4))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_C12_to_He4_He4_He4_derived)*Y(C12));

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
        (-screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(C12)*Y(O16)*state.rho + screened_rates(k_He4_Mg24_to_C12_O16_derived)*Y(He4)*Y(Mg24)*state.rho) +
        (screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(F19)*Y(H1)*state.rho + -screened_rates(k_He4_O16_to_p_F19_derived)*Y(He4)*Y(O16)*state.rho);

    ydot_nuc(O17) =
        (-screened_rates(k_p_O17_to_F18_reaclib)*Y(O17)*Y(H1)*state.rho + screened_rates(k_F18_to_p_O17_derived)*Y(F18)) +
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
        (-screened_rates(k_F17_to_O17_weaktab)*Y(F17) + screened_rates(k_O17_to_F17_weaktab)*Y(O17)) +
        (-screened_rates(k_He4_F17_to_p_Ne20_derived)*Y(F17)*Y(He4)*state.rho + screened_rates(k_p_Ne20_to_He4_F17_reaclib)*Y(Ne20)*Y(H1)*state.rho);

    ydot_nuc(F18) =
        (screened_rates(k_He4_N14_to_F18_reaclib)*Y(He4)*Y(N14)*state.rho + -screened_rates(k_F18_to_He4_N14_derived)*Y(F18)) +
        (screened_rates(k_p_O17_to_F18_reaclib)*Y(O17)*Y(H1)*state.rho + -screened_rates(k_F18_to_p_O17_derived)*Y(F18)) +
        (-screened_rates(k_p_F18_to_Ne19_reaclib)*Y(F18)*Y(H1)*state.rho + screened_rates(k_Ne19_to_p_F18_derived)*Y(Ne19)) +
        (-screened_rates(k_p_F18_to_He4_O15_reaclib)*Y(F18)*Y(H1)*state.rho + screened_rates(k_He4_O15_to_p_F18_derived)*Y(He4)*Y(O15)*state.rho) +
        (-screened_rates(k_F18_to_O18_weaktab)*Y(F18) + screened_rates(k_O18_to_F18_weaktab)*Y(O18)) +
        (screened_rates(k_Ne18_to_F18_weaktab)*Y(Ne18) + -screened_rates(k_F18_to_Ne18_weaktab)*Y(F18));

    ydot_nuc(F19) =
        (screened_rates(k_He4_N15_to_F19_reaclib)*Y(He4)*Y(N15)*state.rho + -screened_rates(k_F19_to_He4_N15_derived)*Y(F19)) +
        (screened_rates(k_p_O18_to_F19_reaclib)*Y(O18)*Y(H1)*state.rho + -screened_rates(k_F19_to_p_O18_derived)*Y(F19)) +
        (-screened_rates(k_p_F19_to_Ne20_reaclib)*Y(F19)*Y(H1)*state.rho + screened_rates(k_Ne20_to_p_F19_derived)*Y(Ne20)) +
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
        (-screened_rates(k_Ne19_to_F19_weaktab)*Y(Ne19) + screened_rates(k_F19_to_Ne19_weaktab)*Y(F19));

    ydot_nuc(Ne20) =
        (screened_rates(k_He4_O16_to_Ne20_reaclib)*Y(He4)*Y(O16)*state.rho + -screened_rates(k_Ne20_to_He4_O16_derived)*Y(Ne20)) +
        (screened_rates(k_p_F19_to_Ne20_reaclib)*Y(F19)*Y(H1)*state.rho + -screened_rates(k_Ne20_to_p_F19_derived)*Y(Ne20)) +
        (-screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(He4)*Y(Ne20)*state.rho + screened_rates(k_Mg24_to_He4_Ne20_derived)*Y(Mg24)) +
        (0.5*screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho + -screened_rates(k_He4_Ne20_to_C12_C12_derived)*Y(He4)*Y(Ne20)*state.rho) +
        (screened_rates(k_He4_F17_to_p_Ne20_derived)*Y(F17)*Y(He4)*state.rho + -screened_rates(k_p_Ne20_to_He4_F17_reaclib)*Y(Ne20)*Y(H1)*state.rho);

    ydot_nuc(Mg22) =
        (screened_rates(k_He4_Ne18_to_Mg22_reaclib)*Y(He4)*Y(Ne18)*state.rho + -screened_rates(k_Mg22_to_He4_Ne18_derived)*Y(Mg22));

    ydot_nuc(Mg24) =
        (screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(He4)*Y(Ne20)*state.rho + -screened_rates(k_Mg24_to_He4_Ne20_derived)*Y(Mg24)) +
        (screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(C12)*Y(O16)*state.rho + -screened_rates(k_He4_Mg24_to_C12_O16_derived)*Y(He4)*Y(Mg24)*state.rho);

    ydot_nuc(Fe56) = 0.0_rt;

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

    scratch = -screened_rates(k_p_C12_to_N13_reaclib)*Y(C12)*state.rho - screened_rates(k_p_C13_to_N14_reaclib)*Y(C13)*state.rho - screened_rates(k_p_F17_to_He4_O14_derived)*Y(F17)*state.rho - screened_rates(k_p_F17_to_Ne18_reaclib)*Y(F17)*state.rho - screened_rates(k_p_F18_to_He4_O15_reaclib)*Y(F18)*state.rho - screened_rates(k_p_F18_to_Ne19_reaclib)*Y(F18)*state.rho - screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(F19)*state.rho - screened_rates(k_p_F19_to_Ne20_reaclib)*Y(F19)*state.rho - screened_rates(k_p_N13_to_O14_reaclib)*Y(N13)*state.rho - screened_rates(k_p_N14_to_O15_reaclib)*Y(N14)*state.rho - screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(N15)*state.rho - screened_rates(k_p_N15_to_O16_reaclib)*Y(N15)*state.rho - screened_rates(k_p_Ne20_to_He4_F17_reaclib)*Y(Ne20)*state.rho - screened_rates(k_p_O16_to_F17_reaclib)*Y(O16)*state.rho - screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*state.rho - screened_rates(k_p_O17_to_F18_reaclib)*Y(O17)*state.rho - screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(O17)*state.rho - screened_rates(k_p_O18_to_F19_reaclib)*Y(O18)*state.rho - screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(O18)*state.rho;
    jac.set(H1, H1, scratch);

    scratch = screened_rates(k_He4_C12_to_p_N15_derived)*Y(C12)*state.rho + screened_rates(k_He4_F17_to_p_Ne20_derived)*Y(F17)*state.rho + screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(N13)*state.rho + screened_rates(k_He4_N14_to_p_O17_derived)*Y(N14)*state.rho + screened_rates(k_He4_N15_to_p_O18_derived)*Y(N15)*state.rho + screened_rates(k_He4_O14_to_p_F17_reaclib)*Y(O14)*state.rho + screened_rates(k_He4_O15_to_p_F18_derived)*Y(O15)*state.rho + screened_rates(k_He4_O16_to_p_F19_derived)*Y(O16)*state.rho;
    jac.set(H1, He4, scratch);

    scratch = screened_rates(k_He4_C12_to_p_N15_derived)*Y(He4)*state.rho - screened_rates(k_p_C12_to_N13_reaclib)*Y(H1)*state.rho;
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

    scratch = screened_rates(k_He4_O16_to_p_F19_derived)*Y(He4)*state.rho + screened_rates(k_O16_to_p_N15_derived) - screened_rates(k_p_O16_to_F17_reaclib)*Y(H1)*state.rho - screened_rates(k_p_O16_to_He4_N13_derived)*Y(H1)*state.rho;
    jac.set(H1, O16, scratch);

    scratch = -screened_rates(k_p_O17_to_F18_reaclib)*Y(H1)*state.rho - screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(H1)*state.rho;
    jac.set(H1, O17, scratch);

    scratch = -screened_rates(k_p_O18_to_F19_reaclib)*Y(H1)*state.rho - screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(H1)*state.rho;
    jac.set(H1, O18, scratch);

    scratch = screened_rates(k_F17_to_p_O16_derived) + screened_rates(k_He4_F17_to_p_Ne20_derived)*Y(He4)*state.rho - screened_rates(k_p_F17_to_He4_O14_derived)*Y(H1)*state.rho - screened_rates(k_p_F17_to_Ne18_reaclib)*Y(H1)*state.rho;
    jac.set(H1, F17, scratch);

    scratch = screened_rates(k_F18_to_p_O17_derived) - screened_rates(k_p_F18_to_He4_O15_reaclib)*Y(H1)*state.rho - screened_rates(k_p_F18_to_Ne19_reaclib)*Y(H1)*state.rho;
    jac.set(H1, F18, scratch);

    scratch = screened_rates(k_F19_to_p_O18_derived) - screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(H1)*state.rho - screened_rates(k_p_F19_to_Ne20_reaclib)*Y(H1)*state.rho;
    jac.set(H1, F19, scratch);

    scratch = screened_rates(k_Ne18_to_p_F17_derived);
    jac.set(H1, Ne18, scratch);

    scratch = screened_rates(k_Ne19_to_p_F18_derived);
    jac.set(H1, Ne19, scratch);

    scratch = screened_rates(k_Ne20_to_p_F19_derived) - screened_rates(k_p_Ne20_to_He4_F17_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Ne20, scratch);

    jac.set(H1, Mg22, 0.0);

    jac.set(H1, Mg24, 0.0);

    jac.set(H1, Fe56, 0.0);

    scratch = screened_rates(k_p_F17_to_He4_O14_derived)*Y(F17)*state.rho + screened_rates(k_p_F18_to_He4_O15_reaclib)*Y(F18)*state.rho + screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(F19)*state.rho + screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(N15)*state.rho + screened_rates(k_p_Ne20_to_He4_F17_reaclib)*Y(Ne20)*state.rho + screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*state.rho + screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(O17)*state.rho + screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(O18)*state.rho;
    jac.set(He4, H1, scratch);

    scratch = -screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*state.rho - screened_rates(k_He4_C12_to_p_N15_derived)*Y(C12)*state.rho - screened_rates(k_He4_F17_to_p_Ne20_derived)*Y(F17)*state.rho - 1.5*screened_rates(k_He4_He4_He4_to_C12_reaclib)*amrex::Math::powi<2>(Y(He4))*amrex::Math::powi<2>(state.rho) - screened_rates(k_He4_Mg24_to_C12_O16_derived)*Y(Mg24)*state.rho - screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(N13)*state.rho - screened_rates(k_He4_N14_to_F18_reaclib)*Y(N14)*state.rho - screened_rates(k_He4_N14_to_p_O17_derived)*Y(N14)*state.rho - screened_rates(k_He4_N15_to_F19_reaclib)*Y(N15)*state.rho - screened_rates(k_He4_N15_to_p_O18_derived)*Y(N15)*state.rho - screened_rates(k_He4_Ne18_to_Mg22_reaclib)*Y(Ne18)*state.rho - screened_rates(k_He4_Ne20_to_C12_C12_derived)*Y(Ne20)*state.rho - screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(Ne20)*state.rho - screened_rates(k_He4_O14_to_Ne18_reaclib)*Y(O14)*state.rho - screened_rates(k_He4_O14_to_p_F17_reaclib)*Y(O14)*state.rho - screened_rates(k_He4_O15_to_Ne19_reaclib)*Y(O15)*state.rho - screened_rates(k_He4_O15_to_p_F18_derived)*Y(O15)*state.rho - screened_rates(k_He4_O16_to_Ne20_reaclib)*Y(O16)*state.rho - screened_rates(k_He4_O16_to_p_F19_derived)*Y(O16)*state.rho;
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

    scratch = screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(C12)*state.rho - screened_rates(k_He4_O16_to_Ne20_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_O16_to_p_F19_derived)*Y(He4)*state.rho + screened_rates(k_O16_to_He4_C12_derived) + screened_rates(k_p_O16_to_He4_N13_derived)*Y(H1)*state.rho;
    jac.set(He4, O16, scratch);

    scratch = screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(H1)*state.rho;
    jac.set(He4, O17, scratch);

    scratch = screened_rates(k_p_O18_to_He4_N15_reaclib)*Y(H1)*state.rho;
    jac.set(He4, O18, scratch);

    scratch = -screened_rates(k_He4_F17_to_p_Ne20_derived)*Y(He4)*state.rho + screened_rates(k_p_F17_to_He4_O14_derived)*Y(H1)*state.rho;
    jac.set(He4, F17, scratch);

    scratch = screened_rates(k_F18_to_He4_N14_derived) + screened_rates(k_p_F18_to_He4_O15_reaclib)*Y(H1)*state.rho;
    jac.set(He4, F18, scratch);

    scratch = screened_rates(k_F19_to_He4_N15_derived) + screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(H1)*state.rho;
    jac.set(He4, F19, scratch);

    scratch = -screened_rates(k_He4_Ne18_to_Mg22_reaclib)*Y(He4)*state.rho + screened_rates(k_Ne18_to_He4_O14_derived);
    jac.set(He4, Ne18, scratch);

    scratch = screened_rates(k_Ne19_to_He4_O15_derived);
    jac.set(He4, Ne19, scratch);

    scratch = -screened_rates(k_He4_Ne20_to_C12_C12_derived)*Y(He4)*state.rho - screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(He4)*state.rho + screened_rates(k_Ne20_to_He4_O16_derived) + screened_rates(k_p_Ne20_to_He4_F17_reaclib)*Y(H1)*state.rho;
    jac.set(He4, Ne20, scratch);

    scratch = screened_rates(k_Mg22_to_He4_Ne18_derived);
    jac.set(He4, Mg22, scratch);

    scratch = -screened_rates(k_He4_Mg24_to_C12_O16_derived)*Y(He4)*state.rho + screened_rates(k_Mg24_to_He4_Ne20_derived);
    jac.set(He4, Mg24, scratch);

    jac.set(He4, Fe56, 0.0);

    scratch = -screened_rates(k_p_C12_to_N13_reaclib)*Y(C12)*state.rho + screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(N15)*state.rho;
    jac.set(C12, H1, scratch);

    scratch = -screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*state.rho - screened_rates(k_He4_C12_to_p_N15_derived)*Y(C12)*state.rho + 0.5*screened_rates(k_He4_He4_He4_to_C12_reaclib)*amrex::Math::powi<2>(Y(He4))*amrex::Math::powi<2>(state.rho) + screened_rates(k_He4_Mg24_to_C12_O16_derived)*Y(Mg24)*state.rho + 2.0*screened_rates(k_He4_Ne20_to_C12_C12_derived)*Y(Ne20)*state.rho;
    jac.set(C12, He4, scratch);

    scratch = -2.0*screened_rates(k_C12_C12_to_He4_Ne20_reaclib)*Y(C12)*state.rho - screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(O16)*state.rho - screened_rates(k_C12_to_He4_He4_He4_derived) - screened_rates(k_He4_C12_to_O16_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_C12_to_p_N15_derived)*Y(He4)*state.rho - screened_rates(k_p_C12_to_N13_reaclib)*Y(H1)*state.rho;
    jac.set(C12, C12, scratch);

    jac.set(C12, C13, 0.0);

    scratch = screened_rates(k_N13_to_p_C12_derived);
    jac.set(C12, N13, scratch);

    jac.set(C12, N14, 0.0);

    scratch = screened_rates(k_p_N15_to_He4_C12_reaclib)*Y(H1)*state.rho;
    jac.set(C12, N15, scratch);

    jac.set(C12, O14, 0.0);

    jac.set(C12, O15, 0.0);

    scratch = -screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(C12)*state.rho + screened_rates(k_O16_to_He4_C12_derived);
    jac.set(C12, O16, scratch);

    jac.set(C12, O17, 0.0);

    jac.set(C12, O18, 0.0);

    jac.set(C12, F17, 0.0);

    jac.set(C12, F18, 0.0);

    jac.set(C12, F19, 0.0);

    jac.set(C12, Ne18, 0.0);

    jac.set(C12, Ne19, 0.0);

    scratch = 2.0*screened_rates(k_He4_Ne20_to_C12_C12_derived)*Y(He4)*state.rho;
    jac.set(C12, Ne20, scratch);

    jac.set(C12, Mg22, 0.0);

    scratch = screened_rates(k_He4_Mg24_to_C12_O16_derived)*Y(He4)*state.rho;
    jac.set(C12, Mg24, scratch);

    jac.set(C12, Fe56, 0.0);

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

    jac.set(C13, Mg22, 0.0);

    jac.set(C13, Mg24, 0.0);

    jac.set(C13, Fe56, 0.0);

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

    jac.set(N13, Mg22, 0.0);

    jac.set(N13, Mg24, 0.0);

    jac.set(N13, Fe56, 0.0);

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

    jac.set(N14, Mg22, 0.0);

    jac.set(N14, Mg24, 0.0);

    jac.set(N14, Fe56, 0.0);

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

    jac.set(N15, Mg22, 0.0);

    jac.set(N15, Mg24, 0.0);

    jac.set(N15, Fe56, 0.0);

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

    jac.set(O14, Mg22, 0.0);

    jac.set(O14, Mg24, 0.0);

    jac.set(O14, Fe56, 0.0);

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

    jac.set(O15, Mg22, 0.0);

    jac.set(O15, Mg24, 0.0);

    jac.set(O15, Fe56, 0.0);

    scratch = screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(F19)*state.rho + screened_rates(k_p_N15_to_O16_reaclib)*Y(N15)*state.rho - screened_rates(k_p_O16_to_F17_reaclib)*Y(O16)*state.rho - screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*state.rho;
    jac.set(O16, H1, scratch);

    scratch = screened_rates(k_He4_C12_to_O16_reaclib)*Y(C12)*state.rho + screened_rates(k_He4_Mg24_to_C12_O16_derived)*Y(Mg24)*state.rho + screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(N13)*state.rho - screened_rates(k_He4_O16_to_Ne20_reaclib)*Y(O16)*state.rho - screened_rates(k_He4_O16_to_p_F19_derived)*Y(O16)*state.rho;
    jac.set(O16, He4, scratch);

    scratch = -screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(O16)*state.rho + screened_rates(k_He4_C12_to_O16_reaclib)*Y(He4)*state.rho;
    jac.set(O16, C12, scratch);

    jac.set(O16, C13, 0.0);

    scratch = screened_rates(k_He4_N13_to_p_O16_reaclib)*Y(He4)*state.rho;
    jac.set(O16, N13, scratch);

    jac.set(O16, N14, 0.0);

    scratch = screened_rates(k_p_N15_to_O16_reaclib)*Y(H1)*state.rho;
    jac.set(O16, N15, scratch);

    jac.set(O16, O14, 0.0);

    jac.set(O16, O15, 0.0);

    scratch = -screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(C12)*state.rho - screened_rates(k_He4_O16_to_Ne20_reaclib)*Y(He4)*state.rho - screened_rates(k_He4_O16_to_p_F19_derived)*Y(He4)*state.rho - screened_rates(k_O16_to_He4_C12_derived) - screened_rates(k_O16_to_p_N15_derived) - screened_rates(k_p_O16_to_F17_reaclib)*Y(H1)*state.rho - screened_rates(k_p_O16_to_He4_N13_derived)*Y(H1)*state.rho;
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

    jac.set(O16, Mg22, 0.0);

    scratch = screened_rates(k_He4_Mg24_to_C12_O16_derived)*Y(He4)*state.rho;
    jac.set(O16, Mg24, scratch);

    jac.set(O16, Fe56, 0.0);

    scratch = -screened_rates(k_p_O17_to_F18_reaclib)*Y(O17)*state.rho - screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(O17)*state.rho;
    jac.set(O17, H1, scratch);

    scratch = screened_rates(k_He4_N14_to_p_O17_derived)*Y(N14)*state.rho;
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

    scratch = -screened_rates(k_O17_to_F17_weaktab) - screened_rates(k_p_O17_to_F18_reaclib)*Y(H1)*state.rho - screened_rates(k_p_O17_to_He4_N14_reaclib)*Y(H1)*state.rho;
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

    jac.set(O17, Mg22, 0.0);

    jac.set(O17, Mg24, 0.0);

    jac.set(O17, Fe56, 0.0);

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

    jac.set(O18, Mg22, 0.0);

    jac.set(O18, Mg24, 0.0);

    jac.set(O18, Fe56, 0.0);

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

    jac.set(F17, Mg22, 0.0);

    jac.set(F17, Mg24, 0.0);

    jac.set(F17, Fe56, 0.0);

    scratch = -screened_rates(k_p_F18_to_He4_O15_reaclib)*Y(F18)*state.rho - screened_rates(k_p_F18_to_Ne19_reaclib)*Y(F18)*state.rho + screened_rates(k_p_O17_to_F18_reaclib)*Y(O17)*state.rho;
    jac.set(F18, H1, scratch);

    scratch = screened_rates(k_He4_N14_to_F18_reaclib)*Y(N14)*state.rho + screened_rates(k_He4_O15_to_p_F18_derived)*Y(O15)*state.rho;
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

    scratch = -screened_rates(k_F18_to_He4_N14_derived) - screened_rates(k_F18_to_Ne18_weaktab) - screened_rates(k_F18_to_O18_weaktab) - screened_rates(k_F18_to_p_O17_derived) - screened_rates(k_p_F18_to_He4_O15_reaclib)*Y(H1)*state.rho - screened_rates(k_p_F18_to_Ne19_reaclib)*Y(H1)*state.rho;
    jac.set(F18, F18, scratch);

    jac.set(F18, F19, 0.0);

    scratch = screened_rates(k_Ne18_to_F18_weaktab);
    jac.set(F18, Ne18, scratch);

    scratch = screened_rates(k_Ne19_to_p_F18_derived);
    jac.set(F18, Ne19, scratch);

    jac.set(F18, Ne20, 0.0);

    jac.set(F18, Mg22, 0.0);

    jac.set(F18, Mg24, 0.0);

    jac.set(F18, Fe56, 0.0);

    scratch = -screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(F19)*state.rho - screened_rates(k_p_F19_to_Ne20_reaclib)*Y(F19)*state.rho + screened_rates(k_p_O18_to_F19_reaclib)*Y(O18)*state.rho;
    jac.set(F19, H1, scratch);

    scratch = screened_rates(k_He4_N15_to_F19_reaclib)*Y(N15)*state.rho + screened_rates(k_He4_O16_to_p_F19_derived)*Y(O16)*state.rho;
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

    scratch = -screened_rates(k_F19_to_He4_N15_derived) - screened_rates(k_F19_to_Ne19_weaktab) - screened_rates(k_F19_to_p_O18_derived) - screened_rates(k_p_F19_to_He4_O16_reaclib)*Y(H1)*state.rho - screened_rates(k_p_F19_to_Ne20_reaclib)*Y(H1)*state.rho;
    jac.set(F19, F19, scratch);

    jac.set(F19, Ne18, 0.0);

    scratch = screened_rates(k_Ne19_to_F19_weaktab);
    jac.set(F19, Ne19, scratch);

    scratch = screened_rates(k_Ne20_to_p_F19_derived);
    jac.set(F19, Ne20, scratch);

    jac.set(F19, Mg22, 0.0);

    jac.set(F19, Mg24, 0.0);

    jac.set(F19, Fe56, 0.0);

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

    scratch = screened_rates(k_Mg22_to_He4_Ne18_derived);
    jac.set(Ne18, Mg22, scratch);

    jac.set(Ne18, Mg24, 0.0);

    jac.set(Ne18, Fe56, 0.0);

    scratch = screened_rates(k_p_F18_to_Ne19_reaclib)*Y(F18)*state.rho;
    jac.set(Ne19, H1, scratch);

    scratch = screened_rates(k_He4_O15_to_Ne19_reaclib)*Y(O15)*state.rho;
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

    scratch = -screened_rates(k_Ne19_to_F19_weaktab) - screened_rates(k_Ne19_to_He4_O15_derived) - screened_rates(k_Ne19_to_p_F18_derived);
    jac.set(Ne19, Ne19, scratch);

    jac.set(Ne19, Ne20, 0.0);

    jac.set(Ne19, Mg22, 0.0);

    jac.set(Ne19, Mg24, 0.0);

    jac.set(Ne19, Fe56, 0.0);

    scratch = screened_rates(k_p_F19_to_Ne20_reaclib)*Y(F19)*state.rho - screened_rates(k_p_Ne20_to_He4_F17_reaclib)*Y(Ne20)*state.rho;
    jac.set(Ne20, H1, scratch);

    scratch = screened_rates(k_He4_F17_to_p_Ne20_derived)*Y(F17)*state.rho - screened_rates(k_He4_Ne20_to_C12_C12_derived)*Y(Ne20)*state.rho - screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(Ne20)*state.rho + screened_rates(k_He4_O16_to_Ne20_reaclib)*Y(O16)*state.rho;
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

    scratch = -screened_rates(k_He4_Ne20_to_C12_C12_derived)*Y(He4)*state.rho - screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(He4)*state.rho - screened_rates(k_Ne20_to_He4_O16_derived) - screened_rates(k_Ne20_to_p_F19_derived) - screened_rates(k_p_Ne20_to_He4_F17_reaclib)*Y(H1)*state.rho;
    jac.set(Ne20, Ne20, scratch);

    jac.set(Ne20, Mg22, 0.0);

    scratch = screened_rates(k_Mg24_to_He4_Ne20_derived);
    jac.set(Ne20, Mg24, scratch);

    jac.set(Ne20, Fe56, 0.0);

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

    scratch = -screened_rates(k_Mg22_to_He4_Ne18_derived);
    jac.set(Mg22, Mg22, scratch);

    jac.set(Mg22, Mg24, 0.0);

    jac.set(Mg22, Fe56, 0.0);

    jac.set(Mg24, H1, 0.0);

    scratch = -screened_rates(k_He4_Mg24_to_C12_O16_derived)*Y(Mg24)*state.rho + screened_rates(k_He4_Ne20_to_Mg24_reaclib)*Y(Ne20)*state.rho;
    jac.set(Mg24, He4, scratch);

    scratch = screened_rates(k_C12_O16_to_He4_Mg24_reaclib)*Y(O16)*state.rho;
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

    jac.set(Mg24, Mg22, 0.0);

    scratch = -screened_rates(k_He4_Mg24_to_C12_O16_derived)*Y(He4)*state.rho - screened_rates(k_Mg24_to_He4_Ne20_derived);
    jac.set(Mg24, Mg24, scratch);

    jac.set(Mg24, Fe56, 0.0);

    jac.set(Fe56, H1, 0.0);

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
