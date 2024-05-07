#ifndef actual_rhs_H
#define actual_rhs_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <extern_parameters.H>
#include <actual_network.H>
#include <burn_type.H>
#include <jacobian_utilities.H>
#include <screen.H>
#include <sneut5.H>
#include <reaclib_rates.H>
#include <table_rates.H>

using namespace amrex;
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

    for (int n = 1; n <= NumSpec; ++n) {
        enuc += dydt(n) * network::mion(n);
    }

    enuc *= C::Legacy::enuc_conv2;
}


template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void evaluate_rates(const burn_t& state, T& rate_eval) {


    // create molar fractions

    amrex::Array1D<amrex::Real, 1, NumSpec> Y;
    for (int n = 1; n <= NumSpec; ++n) {
        Y(n) = state.xn[n-1] * aion_inv[n-1];
    }

    [[maybe_unused]] amrex::Real rhoy = state.rho * state.y_e;

    // Calculate Reaclib rates

    plasma_state_t pstate{};
    fill_plasma_state(pstate, state.T, state.rho, Y);

    tf_t tfactors = evaluate_tfactors(state.T);

    fill_reaclib_rates<do_T_derivatives, T>(tfactors, rate_eval);



    // Evaluate screening factors

    amrex::Real ratraw, dratraw_dT;
    amrex::Real scor, dscor_dt;
    [[maybe_unused]] amrex::Real scor2, dscor2_dt;


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 6.0_rt, 12.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_C12_to_N13);
    rate_eval.screened_rates(k_p_C12_to_N13) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_C12_to_N13);
        rate_eval.dscreened_rates_dT(k_p_C12_to_N13) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 6.0_rt, 12.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_C12_to_O16);
    rate_eval.screened_rates(k_He4_C12_to_O16) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_C12_to_O16);
        rate_eval.dscreened_rates_dT(k_He4_C12_to_O16) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_C12_to_p_N15);
    rate_eval.screened_rates(k_He4_C12_to_p_N15) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_C12_to_p_N15);
        rate_eval.dscreened_rates_dT(k_He4_C12_to_p_N15) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 6.0_rt, 13.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_C13_to_N14);
    rate_eval.screened_rates(k_p_C13_to_N14) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_C13_to_N14);
        rate_eval.dscreened_rates_dT(k_p_C13_to_N14) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 7.0_rt, 13.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_N13_to_O14);
    rate_eval.screened_rates(k_p_N13_to_O14) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_N13_to_O14);
        rate_eval.dscreened_rates_dT(k_p_N13_to_O14) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 7.0_rt, 14.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_N14_to_O15);
    rate_eval.screened_rates(k_p_N14_to_O15) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_N14_to_O15);
        rate_eval.dscreened_rates_dT(k_p_N14_to_O15) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 7.0_rt, 14.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_N14_to_F18);
    rate_eval.screened_rates(k_He4_N14_to_F18) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_N14_to_F18);
        rate_eval.dscreened_rates_dT(k_He4_N14_to_F18) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_N14_to_p_O17);
    rate_eval.screened_rates(k_He4_N14_to_p_O17) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_N14_to_p_O17);
        rate_eval.dscreened_rates_dT(k_He4_N14_to_p_O17) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 7.0_rt, 15.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_N15_to_O16);
    rate_eval.screened_rates(k_p_N15_to_O16) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_N15_to_O16);
        rate_eval.dscreened_rates_dT(k_p_N15_to_O16) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_N15_to_He4_C12);
    rate_eval.screened_rates(k_p_N15_to_He4_C12) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_N15_to_He4_C12);
        rate_eval.dscreened_rates_dT(k_p_N15_to_He4_C12) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 7.0_rt, 15.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_N15_to_F19);
    rate_eval.screened_rates(k_He4_N15_to_F19) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_N15_to_F19);
        rate_eval.dscreened_rates_dT(k_He4_N15_to_F19) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_N15_to_p_O18);
    rate_eval.screened_rates(k_He4_N15_to_p_O18) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_N15_to_p_O18);
        rate_eval.dscreened_rates_dT(k_He4_N15_to_p_O18) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 8.0_rt, 14.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_O14_to_Ne18);
    rate_eval.screened_rates(k_He4_O14_to_Ne18) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_O14_to_Ne18);
        rate_eval.dscreened_rates_dT(k_He4_O14_to_Ne18) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_O14_to_p_F17);
    rate_eval.screened_rates(k_He4_O14_to_p_F17) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_O14_to_p_F17);
        rate_eval.dscreened_rates_dT(k_He4_O14_to_p_F17) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 8.0_rt, 15.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_O15_to_Ne19);
    rate_eval.screened_rates(k_He4_O15_to_Ne19) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_O15_to_Ne19);
        rate_eval.dscreened_rates_dT(k_He4_O15_to_Ne19) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_O15_to_p_F18);
    rate_eval.screened_rates(k_He4_O15_to_p_F18) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_O15_to_p_F18);
        rate_eval.dscreened_rates_dT(k_He4_O15_to_p_F18) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 8.0_rt, 16.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_O16_to_F17);
    rate_eval.screened_rates(k_p_O16_to_F17) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_O16_to_F17);
        rate_eval.dscreened_rates_dT(k_p_O16_to_F17) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_O16_to_He4_N13);
    rate_eval.screened_rates(k_p_O16_to_He4_N13) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_O16_to_He4_N13);
        rate_eval.dscreened_rates_dT(k_p_O16_to_He4_N13) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 8.0_rt, 16.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_O16_to_Ne20);
    rate_eval.screened_rates(k_He4_O16_to_Ne20) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_O16_to_Ne20);
        rate_eval.dscreened_rates_dT(k_He4_O16_to_Ne20) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_O16_to_p_F19);
    rate_eval.screened_rates(k_He4_O16_to_p_F19) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_O16_to_p_F19);
        rate_eval.dscreened_rates_dT(k_He4_O16_to_p_F19) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 8.0_rt, 17.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_O17_to_F18);
    rate_eval.screened_rates(k_p_O17_to_F18) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_O17_to_F18);
        rate_eval.dscreened_rates_dT(k_p_O17_to_F18) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_O17_to_He4_N14);
    rate_eval.screened_rates(k_p_O17_to_He4_N14) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_O17_to_He4_N14);
        rate_eval.dscreened_rates_dT(k_p_O17_to_He4_N14) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 8.0_rt, 18.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_O18_to_F19);
    rate_eval.screened_rates(k_p_O18_to_F19) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_O18_to_F19);
        rate_eval.dscreened_rates_dT(k_p_O18_to_F19) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_O18_to_He4_N15);
    rate_eval.screened_rates(k_p_O18_to_He4_N15) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_O18_to_He4_N15);
        rate_eval.dscreened_rates_dT(k_p_O18_to_He4_N15) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 9.0_rt, 17.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_F17_to_Ne18);
    rate_eval.screened_rates(k_p_F17_to_Ne18) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_F17_to_Ne18);
        rate_eval.dscreened_rates_dT(k_p_F17_to_Ne18) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_F17_to_He4_O14);
    rate_eval.screened_rates(k_p_F17_to_He4_O14) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_F17_to_He4_O14);
        rate_eval.dscreened_rates_dT(k_p_F17_to_He4_O14) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 9.0_rt, 18.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_F18_to_Ne19);
    rate_eval.screened_rates(k_p_F18_to_Ne19) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_F18_to_Ne19);
        rate_eval.dscreened_rates_dT(k_p_F18_to_Ne19) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_F18_to_He4_O15);
    rate_eval.screened_rates(k_p_F18_to_He4_O15) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_F18_to_He4_O15);
        rate_eval.dscreened_rates_dT(k_p_F18_to_He4_O15) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 9.0_rt, 19.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_F19_to_Ne20);
    rate_eval.screened_rates(k_p_F19_to_Ne20) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_F19_to_Ne20);
        rate_eval.dscreened_rates_dT(k_p_F19_to_Ne20) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_F19_to_He4_O16);
    rate_eval.screened_rates(k_p_F19_to_He4_O16) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_F19_to_He4_O16);
        rate_eval.dscreened_rates_dT(k_p_F19_to_He4_O16) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 10.0_rt, 18.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ne18_to_Mg22);
    rate_eval.screened_rates(k_He4_Ne18_to_Mg22) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ne18_to_Mg22);
        rate_eval.dscreened_rates_dT(k_He4_Ne18_to_Mg22) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 10.0_rt, 20.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ne20_to_Mg24);
    rate_eval.screened_rates(k_He4_Ne20_to_Mg24) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ne20_to_Mg24);
        rate_eval.dscreened_rates_dT(k_He4_Ne20_to_Mg24) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ne20_to_C12_C12);
    rate_eval.screened_rates(k_He4_Ne20_to_C12_C12) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ne20_to_C12_C12);
        rate_eval.dscreened_rates_dT(k_He4_Ne20_to_C12_C12) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(6.0_rt, 12.0_rt, 6.0_rt, 12.0_rt);


        static_assert(scn_fac.z1 == 6.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_C12_C12_to_He4_Ne20);
    rate_eval.screened_rates(k_C12_C12_to_He4_Ne20) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_C12_to_He4_Ne20);
        rate_eval.dscreened_rates_dT(k_C12_C12_to_He4_Ne20) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 7.0_rt, 13.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_N13_to_p_O16);
    rate_eval.screened_rates(k_He4_N13_to_p_O16) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_N13_to_p_O16);
        rate_eval.dscreened_rates_dT(k_He4_N13_to_p_O16) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(6.0_rt, 12.0_rt, 8.0_rt, 16.0_rt);


        static_assert(scn_fac.z1 == 6.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_C12_O16_to_He4_Mg24);
    rate_eval.screened_rates(k_C12_O16_to_He4_Mg24) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_O16_to_He4_Mg24);
        rate_eval.dscreened_rates_dT(k_C12_O16_to_He4_Mg24) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 9.0_rt, 17.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_F17_to_p_Ne20);
    rate_eval.screened_rates(k_He4_F17_to_p_Ne20) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_F17_to_p_Ne20);
        rate_eval.dscreened_rates_dT(k_He4_F17_to_p_Ne20) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 10.0_rt, 20.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ne20_to_He4_F17);
    rate_eval.screened_rates(k_p_Ne20_to_He4_F17) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ne20_to_He4_F17);
        rate_eval.dscreened_rates_dT(k_p_Ne20_to_He4_F17) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 12.0_rt, 24.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Mg24_to_C12_O16);
    rate_eval.screened_rates(k_He4_Mg24_to_C12_O16) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mg24_to_C12_O16);
        rate_eval.dscreened_rates_dT(k_He4_Mg24_to_C12_O16) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 2.0_rt, 4.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }



    {
        constexpr auto scn_fac2 = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 4.0_rt, 8.0_rt);


        static_assert(scn_fac2.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac2, scor2, dscor2_dt);

    }


    ratraw = rate_eval.screened_rates(k_He4_He4_He4_to_C12);
    rate_eval.screened_rates(k_He4_He4_He4_to_C12) *= scor * scor2;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_He4_He4_to_C12);
        rate_eval.dscreened_rates_dT(k_He4_He4_He4_to_C12) = ratraw * (scor * dscor2_dt + dscor_dt * scor2) + dratraw_dT * scor * scor2;
    }


    // Fill approximate rates

    fill_approx_rates<do_T_derivatives, T>(tfactors, rate_eval);

    // Calculate tabular rates

    [[maybe_unused]] amrex::Real rate, drate_dt, edot_nu, edot_gamma;

    rate_eval.enuc_weak = 0.0;


}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rhs_nuc(const burn_t& state,
             amrex::Array1D<amrex::Real, 1, neqs>& ydot_nuc,
             const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
             const amrex::Array1D<amrex::Real, 1, NumRates>& screened_rates) {

    using namespace Rates;

    ydot_nuc(H1) =
        (-screened_rates(k_p_C12_to_N13)*Y(C12)*Y(H1)*state.rho + screened_rates(k_N13_to_p_C12)*Y(N13)) +
        (-screened_rates(k_p_C13_to_N14)*Y(C13)*Y(H1)*state.rho + screened_rates(k_N14_to_p_C13)*Y(N14)) +
        (-screened_rates(k_p_N13_to_O14)*Y(N13)*Y(H1)*state.rho + screened_rates(k_O14_to_p_N13)*Y(O14)) +
        (-screened_rates(k_p_N14_to_O15)*Y(N14)*Y(H1)*state.rho + screened_rates(k_O15_to_p_N14)*Y(O15)) +
        (-screened_rates(k_p_N15_to_O16)*Y(N15)*Y(H1)*state.rho + screened_rates(k_O16_to_p_N15)*Y(O16)) +
        (-screened_rates(k_p_O16_to_F17)*Y(O16)*Y(H1)*state.rho + screened_rates(k_F17_to_p_O16)*Y(F17)) +
        (-screened_rates(k_p_O17_to_F18)*Y(O17)*Y(H1)*state.rho + screened_rates(k_F18_to_p_O17)*Y(F18)) +
        (-screened_rates(k_p_O18_to_F19)*Y(O18)*Y(H1)*state.rho + screened_rates(k_F19_to_p_O18)*Y(F19)) +
        (-screened_rates(k_p_F17_to_Ne18)*Y(F17)*Y(H1)*state.rho + screened_rates(k_Ne18_to_p_F17)*Y(Ne18)) +
        (-screened_rates(k_p_F18_to_Ne19)*Y(F18)*Y(H1)*state.rho + screened_rates(k_Ne19_to_p_F18)*Y(Ne19)) +
        (-screened_rates(k_p_F19_to_Ne20)*Y(F19)*Y(H1)*state.rho + screened_rates(k_Ne20_to_p_F19)*Y(Ne20)) +
        (screened_rates(k_He4_N13_to_p_O16)*Y(He4)*Y(N13)*state.rho + -screened_rates(k_p_O16_to_He4_N13)*Y(O16)*Y(H1)*state.rho) +
        (-screened_rates(k_p_N15_to_He4_C12)*Y(N15)*Y(H1)*state.rho + screened_rates(k_He4_C12_to_p_N15)*Y(C12)*Y(He4)*state.rho) +
        (screened_rates(k_He4_O14_to_p_F17)*Y(He4)*Y(O14)*state.rho + -screened_rates(k_p_F17_to_He4_O14)*Y(F17)*Y(H1)*state.rho) +
        (-screened_rates(k_p_O17_to_He4_N14)*Y(O17)*Y(H1)*state.rho + screened_rates(k_He4_N14_to_p_O17)*Y(He4)*Y(N14)*state.rho) +
        (-screened_rates(k_p_O18_to_He4_N15)*Y(O18)*Y(H1)*state.rho + screened_rates(k_He4_N15_to_p_O18)*Y(He4)*Y(N15)*state.rho) +
        (screened_rates(k_He4_F17_to_p_Ne20)*Y(F17)*Y(He4)*state.rho + -screened_rates(k_p_Ne20_to_He4_F17)*Y(Ne20)*Y(H1)*state.rho) +
        (-screened_rates(k_p_F18_to_He4_O15)*Y(F18)*Y(H1)*state.rho + screened_rates(k_He4_O15_to_p_F18)*Y(He4)*Y(O15)*state.rho) +
        (-screened_rates(k_p_F19_to_He4_O16)*Y(F19)*Y(H1)*state.rho + screened_rates(k_He4_O16_to_p_F19)*Y(He4)*Y(O16)*state.rho);

    ydot_nuc(He4) =
        (-screened_rates(k_He4_C12_to_O16)*Y(C12)*Y(He4)*state.rho + screened_rates(k_O16_to_He4_C12)*Y(O16)) +
        (-screened_rates(k_He4_N14_to_F18)*Y(He4)*Y(N14)*state.rho + screened_rates(k_F18_to_He4_N14)*Y(F18)) +
        (-screened_rates(k_He4_N15_to_F19)*Y(He4)*Y(N15)*state.rho + screened_rates(k_F19_to_He4_N15)*Y(F19)) +
        (-screened_rates(k_He4_O14_to_Ne18)*Y(He4)*Y(O14)*state.rho + screened_rates(k_Ne18_to_He4_O14)*Y(Ne18)) +
        (-screened_rates(k_He4_O15_to_Ne19)*Y(He4)*Y(O15)*state.rho + screened_rates(k_Ne19_to_He4_O15)*Y(Ne19)) +
        (-screened_rates(k_He4_O16_to_Ne20)*Y(He4)*Y(O16)*state.rho + screened_rates(k_Ne20_to_He4_O16)*Y(Ne20)) +
        (-screened_rates(k_He4_Ne18_to_Mg22)*Y(He4)*Y(Ne18)*state.rho + screened_rates(k_Mg22_to_He4_Ne18)*Y(Mg22)) +
        (-screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*Y(Ne20)*state.rho + screened_rates(k_Mg24_to_He4_Ne20)*Y(Mg24)) +
        (0.5*screened_rates(k_C12_C12_to_He4_Ne20)*amrex::Math::powi<2>(Y(C12))*state.rho + -screened_rates(k_He4_Ne20_to_C12_C12)*Y(He4)*Y(Ne20)*state.rho) +
        (-screened_rates(k_He4_N13_to_p_O16)*Y(He4)*Y(N13)*state.rho + screened_rates(k_p_O16_to_He4_N13)*Y(O16)*Y(H1)*state.rho) +
        (screened_rates(k_p_N15_to_He4_C12)*Y(N15)*Y(H1)*state.rho + -screened_rates(k_He4_C12_to_p_N15)*Y(C12)*Y(He4)*state.rho) +
        (-screened_rates(k_He4_O14_to_p_F17)*Y(He4)*Y(O14)*state.rho + screened_rates(k_p_F17_to_He4_O14)*Y(F17)*Y(H1)*state.rho) +
        (screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*Y(O16)*state.rho + -screened_rates(k_He4_Mg24_to_C12_O16)*Y(He4)*Y(Mg24)*state.rho) +
        (screened_rates(k_p_O17_to_He4_N14)*Y(O17)*Y(H1)*state.rho + -screened_rates(k_He4_N14_to_p_O17)*Y(He4)*Y(N14)*state.rho) +
        (screened_rates(k_p_O18_to_He4_N15)*Y(O18)*Y(H1)*state.rho + -screened_rates(k_He4_N15_to_p_O18)*Y(He4)*Y(N15)*state.rho) +
        (-screened_rates(k_He4_F17_to_p_Ne20)*Y(F17)*Y(He4)*state.rho + screened_rates(k_p_Ne20_to_He4_F17)*Y(Ne20)*Y(H1)*state.rho) +
        (screened_rates(k_p_F18_to_He4_O15)*Y(F18)*Y(H1)*state.rho + -screened_rates(k_He4_O15_to_p_F18)*Y(He4)*Y(O15)*state.rho) +
        (screened_rates(k_p_F19_to_He4_O16)*Y(F19)*Y(H1)*state.rho + -screened_rates(k_He4_O16_to_p_F19)*Y(He4)*Y(O16)*state.rho) +
        (-0.5*screened_rates(k_He4_He4_He4_to_C12)*amrex::Math::powi<3>(Y(He4))*amrex::Math::powi<2>(state.rho) + 3.0*screened_rates(k_C12_to_He4_He4_He4)*Y(C12));

    ydot_nuc(C12) =
        (-screened_rates(k_p_C12_to_N13)*Y(C12)*Y(H1)*state.rho + screened_rates(k_N13_to_p_C12)*Y(N13)) +
        (-screened_rates(k_He4_C12_to_O16)*Y(C12)*Y(He4)*state.rho + screened_rates(k_O16_to_He4_C12)*Y(O16)) +
        (-screened_rates(k_C12_C12_to_He4_Ne20)*amrex::Math::powi<2>(Y(C12))*state.rho + 2.0*screened_rates(k_He4_Ne20_to_C12_C12)*Y(He4)*Y(Ne20)*state.rho) +
        (screened_rates(k_p_N15_to_He4_C12)*Y(N15)*Y(H1)*state.rho + -screened_rates(k_He4_C12_to_p_N15)*Y(C12)*Y(He4)*state.rho) +
        (-screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*Y(O16)*state.rho + screened_rates(k_He4_Mg24_to_C12_O16)*Y(He4)*Y(Mg24)*state.rho) +
        (0.16666666666666667*screened_rates(k_He4_He4_He4_to_C12)*amrex::Math::powi<3>(Y(He4))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_C12_to_He4_He4_He4)*Y(C12));

    ydot_nuc(C13) =
        screened_rates(k_N13_to_C13_weak_wc12)*Y(N13) +
        (-screened_rates(k_p_C13_to_N14)*Y(C13)*Y(H1)*state.rho + screened_rates(k_N14_to_p_C13)*Y(N14));

    ydot_nuc(N13) =
        -screened_rates(k_N13_to_C13_weak_wc12)*Y(N13) +
        (screened_rates(k_p_C12_to_N13)*Y(C12)*Y(H1)*state.rho + -screened_rates(k_N13_to_p_C12)*Y(N13)) +
        (-screened_rates(k_p_N13_to_O14)*Y(N13)*Y(H1)*state.rho + screened_rates(k_O14_to_p_N13)*Y(O14)) +
        (-screened_rates(k_He4_N13_to_p_O16)*Y(He4)*Y(N13)*state.rho + screened_rates(k_p_O16_to_He4_N13)*Y(O16)*Y(H1)*state.rho);

    ydot_nuc(N14) =
        screened_rates(k_O14_to_N14_weak_wc12)*Y(O14) +
        (screened_rates(k_p_C13_to_N14)*Y(C13)*Y(H1)*state.rho + -screened_rates(k_N14_to_p_C13)*Y(N14)) +
        (-screened_rates(k_p_N14_to_O15)*Y(N14)*Y(H1)*state.rho + screened_rates(k_O15_to_p_N14)*Y(O15)) +
        (-screened_rates(k_He4_N14_to_F18)*Y(He4)*Y(N14)*state.rho + screened_rates(k_F18_to_He4_N14)*Y(F18)) +
        (screened_rates(k_p_O17_to_He4_N14)*Y(O17)*Y(H1)*state.rho + -screened_rates(k_He4_N14_to_p_O17)*Y(He4)*Y(N14)*state.rho);

    ydot_nuc(N15) =
        screened_rates(k_O15_to_N15_weak_wc12)*Y(O15) +
        (-screened_rates(k_p_N15_to_O16)*Y(N15)*Y(H1)*state.rho + screened_rates(k_O16_to_p_N15)*Y(O16)) +
        (-screened_rates(k_He4_N15_to_F19)*Y(He4)*Y(N15)*state.rho + screened_rates(k_F19_to_He4_N15)*Y(F19)) +
        (-screened_rates(k_p_N15_to_He4_C12)*Y(N15)*Y(H1)*state.rho + screened_rates(k_He4_C12_to_p_N15)*Y(C12)*Y(He4)*state.rho) +
        (screened_rates(k_p_O18_to_He4_N15)*Y(O18)*Y(H1)*state.rho + -screened_rates(k_He4_N15_to_p_O18)*Y(He4)*Y(N15)*state.rho);

    ydot_nuc(O14) =
        -screened_rates(k_O14_to_N14_weak_wc12)*Y(O14) +
        (screened_rates(k_p_N13_to_O14)*Y(N13)*Y(H1)*state.rho + -screened_rates(k_O14_to_p_N13)*Y(O14)) +
        (-screened_rates(k_He4_O14_to_Ne18)*Y(He4)*Y(O14)*state.rho + screened_rates(k_Ne18_to_He4_O14)*Y(Ne18)) +
        (-screened_rates(k_He4_O14_to_p_F17)*Y(He4)*Y(O14)*state.rho + screened_rates(k_p_F17_to_He4_O14)*Y(F17)*Y(H1)*state.rho);

    ydot_nuc(O15) =
        -screened_rates(k_O15_to_N15_weak_wc12)*Y(O15) +
        (screened_rates(k_p_N14_to_O15)*Y(N14)*Y(H1)*state.rho + -screened_rates(k_O15_to_p_N14)*Y(O15)) +
        (-screened_rates(k_He4_O15_to_Ne19)*Y(He4)*Y(O15)*state.rho + screened_rates(k_Ne19_to_He4_O15)*Y(Ne19)) +
        (screened_rates(k_p_F18_to_He4_O15)*Y(F18)*Y(H1)*state.rho + -screened_rates(k_He4_O15_to_p_F18)*Y(He4)*Y(O15)*state.rho);

    ydot_nuc(O16) =
        (screened_rates(k_He4_C12_to_O16)*Y(C12)*Y(He4)*state.rho + -screened_rates(k_O16_to_He4_C12)*Y(O16)) +
        (screened_rates(k_p_N15_to_O16)*Y(N15)*Y(H1)*state.rho + -screened_rates(k_O16_to_p_N15)*Y(O16)) +
        (-screened_rates(k_p_O16_to_F17)*Y(O16)*Y(H1)*state.rho + screened_rates(k_F17_to_p_O16)*Y(F17)) +
        (-screened_rates(k_He4_O16_to_Ne20)*Y(He4)*Y(O16)*state.rho + screened_rates(k_Ne20_to_He4_O16)*Y(Ne20)) +
        (screened_rates(k_He4_N13_to_p_O16)*Y(He4)*Y(N13)*state.rho + -screened_rates(k_p_O16_to_He4_N13)*Y(O16)*Y(H1)*state.rho) +
        (-screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*Y(O16)*state.rho + screened_rates(k_He4_Mg24_to_C12_O16)*Y(He4)*Y(Mg24)*state.rho) +
        (screened_rates(k_p_F19_to_He4_O16)*Y(F19)*Y(H1)*state.rho + -screened_rates(k_He4_O16_to_p_F19)*Y(He4)*Y(O16)*state.rho);

    ydot_nuc(O17) =
        screened_rates(k_F17_to_O17_weak_wc12)*Y(F17) +
        (-screened_rates(k_p_O17_to_F18)*Y(O17)*Y(H1)*state.rho + screened_rates(k_F18_to_p_O17)*Y(F18)) +
        (-screened_rates(k_p_O17_to_He4_N14)*Y(O17)*Y(H1)*state.rho + screened_rates(k_He4_N14_to_p_O17)*Y(He4)*Y(N14)*state.rho);

    ydot_nuc(O18) =
        screened_rates(k_F18_to_O18_weak_wc12)*Y(F18) +
        (-screened_rates(k_p_O18_to_F19)*Y(O18)*Y(H1)*state.rho + screened_rates(k_F19_to_p_O18)*Y(F19)) +
        (-screened_rates(k_p_O18_to_He4_N15)*Y(O18)*Y(H1)*state.rho + screened_rates(k_He4_N15_to_p_O18)*Y(He4)*Y(N15)*state.rho);

    ydot_nuc(F17) =
        -screened_rates(k_F17_to_O17_weak_wc12)*Y(F17) +
        (screened_rates(k_p_O16_to_F17)*Y(O16)*Y(H1)*state.rho + -screened_rates(k_F17_to_p_O16)*Y(F17)) +
        (-screened_rates(k_p_F17_to_Ne18)*Y(F17)*Y(H1)*state.rho + screened_rates(k_Ne18_to_p_F17)*Y(Ne18)) +
        (screened_rates(k_He4_O14_to_p_F17)*Y(He4)*Y(O14)*state.rho + -screened_rates(k_p_F17_to_He4_O14)*Y(F17)*Y(H1)*state.rho) +
        (-screened_rates(k_He4_F17_to_p_Ne20)*Y(F17)*Y(He4)*state.rho + screened_rates(k_p_Ne20_to_He4_F17)*Y(Ne20)*Y(H1)*state.rho);

    ydot_nuc(F18) =
        -screened_rates(k_F18_to_O18_weak_wc12)*Y(F18) +
        screened_rates(k_Ne18_to_F18_weak_wc12)*Y(Ne18) +
        (screened_rates(k_He4_N14_to_F18)*Y(He4)*Y(N14)*state.rho + -screened_rates(k_F18_to_He4_N14)*Y(F18)) +
        (screened_rates(k_p_O17_to_F18)*Y(O17)*Y(H1)*state.rho + -screened_rates(k_F18_to_p_O17)*Y(F18)) +
        (-screened_rates(k_p_F18_to_Ne19)*Y(F18)*Y(H1)*state.rho + screened_rates(k_Ne19_to_p_F18)*Y(Ne19)) +
        (-screened_rates(k_p_F18_to_He4_O15)*Y(F18)*Y(H1)*state.rho + screened_rates(k_He4_O15_to_p_F18)*Y(He4)*Y(O15)*state.rho);

    ydot_nuc(F19) =
        screened_rates(k_Ne19_to_F19_weak_wc12)*Y(Ne19) +
        (screened_rates(k_He4_N15_to_F19)*Y(He4)*Y(N15)*state.rho + -screened_rates(k_F19_to_He4_N15)*Y(F19)) +
        (screened_rates(k_p_O18_to_F19)*Y(O18)*Y(H1)*state.rho + -screened_rates(k_F19_to_p_O18)*Y(F19)) +
        (-screened_rates(k_p_F19_to_Ne20)*Y(F19)*Y(H1)*state.rho + screened_rates(k_Ne20_to_p_F19)*Y(Ne20)) +
        (-screened_rates(k_p_F19_to_He4_O16)*Y(F19)*Y(H1)*state.rho + screened_rates(k_He4_O16_to_p_F19)*Y(He4)*Y(O16)*state.rho);

    ydot_nuc(Ne18) =
        -screened_rates(k_Ne18_to_F18_weak_wc12)*Y(Ne18) +
        (screened_rates(k_He4_O14_to_Ne18)*Y(He4)*Y(O14)*state.rho + -screened_rates(k_Ne18_to_He4_O14)*Y(Ne18)) +
        (screened_rates(k_p_F17_to_Ne18)*Y(F17)*Y(H1)*state.rho + -screened_rates(k_Ne18_to_p_F17)*Y(Ne18)) +
        (-screened_rates(k_He4_Ne18_to_Mg22)*Y(He4)*Y(Ne18)*state.rho + screened_rates(k_Mg22_to_He4_Ne18)*Y(Mg22));

    ydot_nuc(Ne19) =
        -screened_rates(k_Ne19_to_F19_weak_wc12)*Y(Ne19) +
        (screened_rates(k_He4_O15_to_Ne19)*Y(He4)*Y(O15)*state.rho + -screened_rates(k_Ne19_to_He4_O15)*Y(Ne19)) +
        (screened_rates(k_p_F18_to_Ne19)*Y(F18)*Y(H1)*state.rho + -screened_rates(k_Ne19_to_p_F18)*Y(Ne19));

    ydot_nuc(Ne20) =
        (screened_rates(k_He4_O16_to_Ne20)*Y(He4)*Y(O16)*state.rho + -screened_rates(k_Ne20_to_He4_O16)*Y(Ne20)) +
        (screened_rates(k_p_F19_to_Ne20)*Y(F19)*Y(H1)*state.rho + -screened_rates(k_Ne20_to_p_F19)*Y(Ne20)) +
        (-screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*Y(Ne20)*state.rho + screened_rates(k_Mg24_to_He4_Ne20)*Y(Mg24)) +
        (0.5*screened_rates(k_C12_C12_to_He4_Ne20)*amrex::Math::powi<2>(Y(C12))*state.rho + -screened_rates(k_He4_Ne20_to_C12_C12)*Y(He4)*Y(Ne20)*state.rho) +
        (screened_rates(k_He4_F17_to_p_Ne20)*Y(F17)*Y(He4)*state.rho + -screened_rates(k_p_Ne20_to_He4_F17)*Y(Ne20)*Y(H1)*state.rho);

    ydot_nuc(Mg22) =
        (screened_rates(k_He4_Ne18_to_Mg22)*Y(He4)*Y(Ne18)*state.rho + -screened_rates(k_Mg22_to_He4_Ne18)*Y(Mg22));

    ydot_nuc(Mg24) =
        (screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*Y(Ne20)*state.rho + -screened_rates(k_Mg24_to_He4_Ne20)*Y(Mg24)) +
        (screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*Y(O16)*state.rho + -screened_rates(k_He4_Mg24_to_C12_O16)*Y(He4)*Y(Mg24)*state.rho);

    ydot_nuc(Fe56) = 0.0;

}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_rhs (burn_t& state, amrex::Array1D<amrex::Real, 1, neqs>& ydot)
{
    for (int i = 1; i <= neqs; ++i) {
        ydot(i) = 0.0_rt;
    }


    // Set molar abundances
    amrex::Array1D<amrex::Real, 1, NumSpec> Y;
    for (int i = 1; i <= NumSpec; ++i) {
        Y(i) = state.xn[i-1] * aion_inv[i-1];
    }

    // build the rates

    rate_t rate_eval;

    constexpr int do_T_derivatives = 0;

    evaluate_rates<do_T_derivatives, rate_t>(state, rate_eval);

    rhs_nuc(state, ydot, Y, rate_eval.screened_rates);

    // ion binding energy contributions

    amrex::Real enuc;
    ener_gener_rate(ydot, enuc);

    // include any weak rate neutrino losses
    enuc += rate_eval.enuc_weak;

    // Get the thermal neutrino losses

    amrex::Real sneut, dsneutdt, dsneutdd, dsnuda, dsnudz;
    constexpr int do_derivatives{0};
    sneut5<do_derivatives>(state.T, state.rho, state.abar, state.zbar, sneut, dsneutdt, dsneutdd, dsnuda, dsnudz);

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

    scratch = -screened_rates(k_p_C12_to_N13)*Y(C12)*state.rho - screened_rates(k_p_C13_to_N14)*Y(C13)*state.rho - screened_rates(k_p_F17_to_He4_O14)*Y(F17)*state.rho - screened_rates(k_p_F17_to_Ne18)*Y(F17)*state.rho - screened_rates(k_p_F18_to_He4_O15)*Y(F18)*state.rho - screened_rates(k_p_F18_to_Ne19)*Y(F18)*state.rho - screened_rates(k_p_F19_to_He4_O16)*Y(F19)*state.rho - screened_rates(k_p_F19_to_Ne20)*Y(F19)*state.rho - screened_rates(k_p_N13_to_O14)*Y(N13)*state.rho - screened_rates(k_p_N14_to_O15)*Y(N14)*state.rho - screened_rates(k_p_N15_to_He4_C12)*Y(N15)*state.rho - screened_rates(k_p_N15_to_O16)*Y(N15)*state.rho - screened_rates(k_p_Ne20_to_He4_F17)*Y(Ne20)*state.rho - screened_rates(k_p_O16_to_F17)*Y(O16)*state.rho - screened_rates(k_p_O16_to_He4_N13)*Y(O16)*state.rho - screened_rates(k_p_O17_to_F18)*Y(O17)*state.rho - screened_rates(k_p_O17_to_He4_N14)*Y(O17)*state.rho - screened_rates(k_p_O18_to_F19)*Y(O18)*state.rho - screened_rates(k_p_O18_to_He4_N15)*Y(O18)*state.rho;
    jac.set(H1, H1, scratch);

    scratch = screened_rates(k_He4_C12_to_p_N15)*Y(C12)*state.rho + screened_rates(k_He4_F17_to_p_Ne20)*Y(F17)*state.rho + screened_rates(k_He4_N13_to_p_O16)*Y(N13)*state.rho + screened_rates(k_He4_N14_to_p_O17)*Y(N14)*state.rho + screened_rates(k_He4_N15_to_p_O18)*Y(N15)*state.rho + screened_rates(k_He4_O14_to_p_F17)*Y(O14)*state.rho + screened_rates(k_He4_O15_to_p_F18)*Y(O15)*state.rho + screened_rates(k_He4_O16_to_p_F19)*Y(O16)*state.rho;
    jac.set(H1, He4, scratch);

    scratch = screened_rates(k_He4_C12_to_p_N15)*Y(He4)*state.rho - screened_rates(k_p_C12_to_N13)*Y(H1)*state.rho;
    jac.set(H1, C12, scratch);

    scratch = -screened_rates(k_p_C13_to_N14)*Y(H1)*state.rho;
    jac.set(H1, C13, scratch);

    scratch = screened_rates(k_He4_N13_to_p_O16)*Y(He4)*state.rho + screened_rates(k_N13_to_p_C12) - screened_rates(k_p_N13_to_O14)*Y(H1)*state.rho;
    jac.set(H1, N13, scratch);

    scratch = screened_rates(k_He4_N14_to_p_O17)*Y(He4)*state.rho + screened_rates(k_N14_to_p_C13) - screened_rates(k_p_N14_to_O15)*Y(H1)*state.rho;
    jac.set(H1, N14, scratch);

    scratch = screened_rates(k_He4_N15_to_p_O18)*Y(He4)*state.rho - screened_rates(k_p_N15_to_He4_C12)*Y(H1)*state.rho - screened_rates(k_p_N15_to_O16)*Y(H1)*state.rho;
    jac.set(H1, N15, scratch);

    scratch = screened_rates(k_He4_O14_to_p_F17)*Y(He4)*state.rho + screened_rates(k_O14_to_p_N13);
    jac.set(H1, O14, scratch);

    scratch = screened_rates(k_He4_O15_to_p_F18)*Y(He4)*state.rho + screened_rates(k_O15_to_p_N14);
    jac.set(H1, O15, scratch);

    scratch = screened_rates(k_He4_O16_to_p_F19)*Y(He4)*state.rho + screened_rates(k_O16_to_p_N15) - screened_rates(k_p_O16_to_F17)*Y(H1)*state.rho - screened_rates(k_p_O16_to_He4_N13)*Y(H1)*state.rho;
    jac.set(H1, O16, scratch);

    scratch = -screened_rates(k_p_O17_to_F18)*Y(H1)*state.rho - screened_rates(k_p_O17_to_He4_N14)*Y(H1)*state.rho;
    jac.set(H1, O17, scratch);

    scratch = -screened_rates(k_p_O18_to_F19)*Y(H1)*state.rho - screened_rates(k_p_O18_to_He4_N15)*Y(H1)*state.rho;
    jac.set(H1, O18, scratch);

    scratch = screened_rates(k_F17_to_p_O16) + screened_rates(k_He4_F17_to_p_Ne20)*Y(He4)*state.rho - screened_rates(k_p_F17_to_He4_O14)*Y(H1)*state.rho - screened_rates(k_p_F17_to_Ne18)*Y(H1)*state.rho;
    jac.set(H1, F17, scratch);

    scratch = screened_rates(k_F18_to_p_O17) - screened_rates(k_p_F18_to_He4_O15)*Y(H1)*state.rho - screened_rates(k_p_F18_to_Ne19)*Y(H1)*state.rho;
    jac.set(H1, F18, scratch);

    scratch = screened_rates(k_F19_to_p_O18) - screened_rates(k_p_F19_to_He4_O16)*Y(H1)*state.rho - screened_rates(k_p_F19_to_Ne20)*Y(H1)*state.rho;
    jac.set(H1, F19, scratch);

    scratch = screened_rates(k_Ne18_to_p_F17);
    jac.set(H1, Ne18, scratch);

    scratch = screened_rates(k_Ne19_to_p_F18);
    jac.set(H1, Ne19, scratch);

    scratch = screened_rates(k_Ne20_to_p_F19) - screened_rates(k_p_Ne20_to_He4_F17)*Y(H1)*state.rho;
    jac.set(H1, Ne20, scratch);

    scratch = screened_rates(k_p_F17_to_He4_O14)*Y(F17)*state.rho + screened_rates(k_p_F18_to_He4_O15)*Y(F18)*state.rho + screened_rates(k_p_F19_to_He4_O16)*Y(F19)*state.rho + screened_rates(k_p_N15_to_He4_C12)*Y(N15)*state.rho + screened_rates(k_p_Ne20_to_He4_F17)*Y(Ne20)*state.rho + screened_rates(k_p_O16_to_He4_N13)*Y(O16)*state.rho + screened_rates(k_p_O17_to_He4_N14)*Y(O17)*state.rho + screened_rates(k_p_O18_to_He4_N15)*Y(O18)*state.rho;
    jac.set(He4, H1, scratch);

    scratch = -screened_rates(k_He4_C12_to_O16)*Y(C12)*state.rho - screened_rates(k_He4_C12_to_p_N15)*Y(C12)*state.rho - screened_rates(k_He4_F17_to_p_Ne20)*Y(F17)*state.rho - 1.5*screened_rates(k_He4_He4_He4_to_C12)*amrex::Math::powi<2>(Y(He4))*amrex::Math::powi<2>(state.rho) - screened_rates(k_He4_Mg24_to_C12_O16)*Y(Mg24)*state.rho - screened_rates(k_He4_N13_to_p_O16)*Y(N13)*state.rho - screened_rates(k_He4_N14_to_F18)*Y(N14)*state.rho - screened_rates(k_He4_N14_to_p_O17)*Y(N14)*state.rho - screened_rates(k_He4_N15_to_F19)*Y(N15)*state.rho - screened_rates(k_He4_N15_to_p_O18)*Y(N15)*state.rho - screened_rates(k_He4_Ne18_to_Mg22)*Y(Ne18)*state.rho - screened_rates(k_He4_Ne20_to_C12_C12)*Y(Ne20)*state.rho - screened_rates(k_He4_Ne20_to_Mg24)*Y(Ne20)*state.rho - screened_rates(k_He4_O14_to_Ne18)*Y(O14)*state.rho - screened_rates(k_He4_O14_to_p_F17)*Y(O14)*state.rho - screened_rates(k_He4_O15_to_Ne19)*Y(O15)*state.rho - screened_rates(k_He4_O15_to_p_F18)*Y(O15)*state.rho - screened_rates(k_He4_O16_to_Ne20)*Y(O16)*state.rho - screened_rates(k_He4_O16_to_p_F19)*Y(O16)*state.rho;
    jac.set(He4, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_He4_Ne20)*Y(C12)*state.rho + screened_rates(k_C12_O16_to_He4_Mg24)*Y(O16)*state.rho + 3.0*screened_rates(k_C12_to_He4_He4_He4) - screened_rates(k_He4_C12_to_O16)*Y(He4)*state.rho - screened_rates(k_He4_C12_to_p_N15)*Y(He4)*state.rho;
    jac.set(He4, C12, scratch);

    scratch = -screened_rates(k_He4_N13_to_p_O16)*Y(He4)*state.rho;
    jac.set(He4, N13, scratch);

    scratch = -screened_rates(k_He4_N14_to_F18)*Y(He4)*state.rho - screened_rates(k_He4_N14_to_p_O17)*Y(He4)*state.rho;
    jac.set(He4, N14, scratch);

    scratch = -screened_rates(k_He4_N15_to_F19)*Y(He4)*state.rho - screened_rates(k_He4_N15_to_p_O18)*Y(He4)*state.rho + screened_rates(k_p_N15_to_He4_C12)*Y(H1)*state.rho;
    jac.set(He4, N15, scratch);

    scratch = -screened_rates(k_He4_O14_to_Ne18)*Y(He4)*state.rho - screened_rates(k_He4_O14_to_p_F17)*Y(He4)*state.rho;
    jac.set(He4, O14, scratch);

    scratch = -screened_rates(k_He4_O15_to_Ne19)*Y(He4)*state.rho - screened_rates(k_He4_O15_to_p_F18)*Y(He4)*state.rho;
    jac.set(He4, O15, scratch);

    scratch = screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*state.rho - screened_rates(k_He4_O16_to_Ne20)*Y(He4)*state.rho - screened_rates(k_He4_O16_to_p_F19)*Y(He4)*state.rho + screened_rates(k_O16_to_He4_C12) + screened_rates(k_p_O16_to_He4_N13)*Y(H1)*state.rho;
    jac.set(He4, O16, scratch);

    scratch = screened_rates(k_p_O17_to_He4_N14)*Y(H1)*state.rho;
    jac.set(He4, O17, scratch);

    scratch = screened_rates(k_p_O18_to_He4_N15)*Y(H1)*state.rho;
    jac.set(He4, O18, scratch);

    scratch = -screened_rates(k_He4_F17_to_p_Ne20)*Y(He4)*state.rho + screened_rates(k_p_F17_to_He4_O14)*Y(H1)*state.rho;
    jac.set(He4, F17, scratch);

    scratch = screened_rates(k_F18_to_He4_N14) + screened_rates(k_p_F18_to_He4_O15)*Y(H1)*state.rho;
    jac.set(He4, F18, scratch);

    scratch = screened_rates(k_F19_to_He4_N15) + screened_rates(k_p_F19_to_He4_O16)*Y(H1)*state.rho;
    jac.set(He4, F19, scratch);

    scratch = -screened_rates(k_He4_Ne18_to_Mg22)*Y(He4)*state.rho + screened_rates(k_Ne18_to_He4_O14);
    jac.set(He4, Ne18, scratch);

    scratch = screened_rates(k_Ne19_to_He4_O15);
    jac.set(He4, Ne19, scratch);

    scratch = -screened_rates(k_He4_Ne20_to_C12_C12)*Y(He4)*state.rho - screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*state.rho + screened_rates(k_Ne20_to_He4_O16) + screened_rates(k_p_Ne20_to_He4_F17)*Y(H1)*state.rho;
    jac.set(He4, Ne20, scratch);

    scratch = screened_rates(k_Mg22_to_He4_Ne18);
    jac.set(He4, Mg22, scratch);

    scratch = -screened_rates(k_He4_Mg24_to_C12_O16)*Y(He4)*state.rho + screened_rates(k_Mg24_to_He4_Ne20);
    jac.set(He4, Mg24, scratch);

    scratch = -screened_rates(k_p_C12_to_N13)*Y(C12)*state.rho + screened_rates(k_p_N15_to_He4_C12)*Y(N15)*state.rho;
    jac.set(C12, H1, scratch);

    scratch = -screened_rates(k_He4_C12_to_O16)*Y(C12)*state.rho - screened_rates(k_He4_C12_to_p_N15)*Y(C12)*state.rho + 0.5*screened_rates(k_He4_He4_He4_to_C12)*amrex::Math::powi<2>(Y(He4))*amrex::Math::powi<2>(state.rho) + screened_rates(k_He4_Mg24_to_C12_O16)*Y(Mg24)*state.rho + 2.0*screened_rates(k_He4_Ne20_to_C12_C12)*Y(Ne20)*state.rho;
    jac.set(C12, He4, scratch);

    scratch = -2.0*screened_rates(k_C12_C12_to_He4_Ne20)*Y(C12)*state.rho - screened_rates(k_C12_O16_to_He4_Mg24)*Y(O16)*state.rho - screened_rates(k_C12_to_He4_He4_He4) - screened_rates(k_He4_C12_to_O16)*Y(He4)*state.rho - screened_rates(k_He4_C12_to_p_N15)*Y(He4)*state.rho - screened_rates(k_p_C12_to_N13)*Y(H1)*state.rho;
    jac.set(C12, C12, scratch);

    scratch = screened_rates(k_N13_to_p_C12);
    jac.set(C12, N13, scratch);

    scratch = screened_rates(k_p_N15_to_He4_C12)*Y(H1)*state.rho;
    jac.set(C12, N15, scratch);

    scratch = -screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*state.rho + screened_rates(k_O16_to_He4_C12);
    jac.set(C12, O16, scratch);

    scratch = 2.0*screened_rates(k_He4_Ne20_to_C12_C12)*Y(He4)*state.rho;
    jac.set(C12, Ne20, scratch);

    scratch = screened_rates(k_He4_Mg24_to_C12_O16)*Y(He4)*state.rho;
    jac.set(C12, Mg24, scratch);

    scratch = -screened_rates(k_p_C13_to_N14)*Y(C13)*state.rho;
    jac.set(C13, H1, scratch);

    scratch = -screened_rates(k_p_C13_to_N14)*Y(H1)*state.rho;
    jac.set(C13, C13, scratch);

    scratch = screened_rates(k_N13_to_C13_weak_wc12);
    jac.set(C13, N13, scratch);

    scratch = screened_rates(k_N14_to_p_C13);
    jac.set(C13, N14, scratch);

    scratch = screened_rates(k_p_C12_to_N13)*Y(C12)*state.rho - screened_rates(k_p_N13_to_O14)*Y(N13)*state.rho + screened_rates(k_p_O16_to_He4_N13)*Y(O16)*state.rho;
    jac.set(N13, H1, scratch);

    scratch = -screened_rates(k_He4_N13_to_p_O16)*Y(N13)*state.rho;
    jac.set(N13, He4, scratch);

    scratch = screened_rates(k_p_C12_to_N13)*Y(H1)*state.rho;
    jac.set(N13, C12, scratch);

    scratch = -screened_rates(k_He4_N13_to_p_O16)*Y(He4)*state.rho - screened_rates(k_N13_to_C13_weak_wc12) - screened_rates(k_N13_to_p_C12) - screened_rates(k_p_N13_to_O14)*Y(H1)*state.rho;
    jac.set(N13, N13, scratch);

    scratch = screened_rates(k_O14_to_p_N13);
    jac.set(N13, O14, scratch);

    scratch = screened_rates(k_p_O16_to_He4_N13)*Y(H1)*state.rho;
    jac.set(N13, O16, scratch);

    scratch = screened_rates(k_p_C13_to_N14)*Y(C13)*state.rho - screened_rates(k_p_N14_to_O15)*Y(N14)*state.rho + screened_rates(k_p_O17_to_He4_N14)*Y(O17)*state.rho;
    jac.set(N14, H1, scratch);

    scratch = -screened_rates(k_He4_N14_to_F18)*Y(N14)*state.rho - screened_rates(k_He4_N14_to_p_O17)*Y(N14)*state.rho;
    jac.set(N14, He4, scratch);

    scratch = screened_rates(k_p_C13_to_N14)*Y(H1)*state.rho;
    jac.set(N14, C13, scratch);

    scratch = -screened_rates(k_He4_N14_to_F18)*Y(He4)*state.rho - screened_rates(k_He4_N14_to_p_O17)*Y(He4)*state.rho - screened_rates(k_N14_to_p_C13) - screened_rates(k_p_N14_to_O15)*Y(H1)*state.rho;
    jac.set(N14, N14, scratch);

    scratch = screened_rates(k_O14_to_N14_weak_wc12);
    jac.set(N14, O14, scratch);

    scratch = screened_rates(k_O15_to_p_N14);
    jac.set(N14, O15, scratch);

    scratch = screened_rates(k_p_O17_to_He4_N14)*Y(H1)*state.rho;
    jac.set(N14, O17, scratch);

    scratch = screened_rates(k_F18_to_He4_N14);
    jac.set(N14, F18, scratch);

    scratch = -screened_rates(k_p_N15_to_He4_C12)*Y(N15)*state.rho - screened_rates(k_p_N15_to_O16)*Y(N15)*state.rho + screened_rates(k_p_O18_to_He4_N15)*Y(O18)*state.rho;
    jac.set(N15, H1, scratch);

    scratch = screened_rates(k_He4_C12_to_p_N15)*Y(C12)*state.rho - screened_rates(k_He4_N15_to_F19)*Y(N15)*state.rho - screened_rates(k_He4_N15_to_p_O18)*Y(N15)*state.rho;
    jac.set(N15, He4, scratch);

    scratch = screened_rates(k_He4_C12_to_p_N15)*Y(He4)*state.rho;
    jac.set(N15, C12, scratch);

    scratch = -screened_rates(k_He4_N15_to_F19)*Y(He4)*state.rho - screened_rates(k_He4_N15_to_p_O18)*Y(He4)*state.rho - screened_rates(k_p_N15_to_He4_C12)*Y(H1)*state.rho - screened_rates(k_p_N15_to_O16)*Y(H1)*state.rho;
    jac.set(N15, N15, scratch);

    scratch = screened_rates(k_O15_to_N15_weak_wc12);
    jac.set(N15, O15, scratch);

    scratch = screened_rates(k_O16_to_p_N15);
    jac.set(N15, O16, scratch);

    scratch = screened_rates(k_p_O18_to_He4_N15)*Y(H1)*state.rho;
    jac.set(N15, O18, scratch);

    scratch = screened_rates(k_F19_to_He4_N15);
    jac.set(N15, F19, scratch);

    scratch = screened_rates(k_p_F17_to_He4_O14)*Y(F17)*state.rho + screened_rates(k_p_N13_to_O14)*Y(N13)*state.rho;
    jac.set(O14, H1, scratch);

    scratch = -screened_rates(k_He4_O14_to_Ne18)*Y(O14)*state.rho - screened_rates(k_He4_O14_to_p_F17)*Y(O14)*state.rho;
    jac.set(O14, He4, scratch);

    scratch = screened_rates(k_p_N13_to_O14)*Y(H1)*state.rho;
    jac.set(O14, N13, scratch);

    scratch = -screened_rates(k_He4_O14_to_Ne18)*Y(He4)*state.rho - screened_rates(k_He4_O14_to_p_F17)*Y(He4)*state.rho - screened_rates(k_O14_to_N14_weak_wc12) - screened_rates(k_O14_to_p_N13);
    jac.set(O14, O14, scratch);

    scratch = screened_rates(k_p_F17_to_He4_O14)*Y(H1)*state.rho;
    jac.set(O14, F17, scratch);

    scratch = screened_rates(k_Ne18_to_He4_O14);
    jac.set(O14, Ne18, scratch);

    scratch = screened_rates(k_p_F18_to_He4_O15)*Y(F18)*state.rho + screened_rates(k_p_N14_to_O15)*Y(N14)*state.rho;
    jac.set(O15, H1, scratch);

    scratch = -screened_rates(k_He4_O15_to_Ne19)*Y(O15)*state.rho - screened_rates(k_He4_O15_to_p_F18)*Y(O15)*state.rho;
    jac.set(O15, He4, scratch);

    scratch = screened_rates(k_p_N14_to_O15)*Y(H1)*state.rho;
    jac.set(O15, N14, scratch);

    scratch = -screened_rates(k_He4_O15_to_Ne19)*Y(He4)*state.rho - screened_rates(k_He4_O15_to_p_F18)*Y(He4)*state.rho - screened_rates(k_O15_to_N15_weak_wc12) - screened_rates(k_O15_to_p_N14);
    jac.set(O15, O15, scratch);

    scratch = screened_rates(k_p_F18_to_He4_O15)*Y(H1)*state.rho;
    jac.set(O15, F18, scratch);

    scratch = screened_rates(k_Ne19_to_He4_O15);
    jac.set(O15, Ne19, scratch);

    scratch = screened_rates(k_p_F19_to_He4_O16)*Y(F19)*state.rho + screened_rates(k_p_N15_to_O16)*Y(N15)*state.rho - screened_rates(k_p_O16_to_F17)*Y(O16)*state.rho - screened_rates(k_p_O16_to_He4_N13)*Y(O16)*state.rho;
    jac.set(O16, H1, scratch);

    scratch = screened_rates(k_He4_C12_to_O16)*Y(C12)*state.rho + screened_rates(k_He4_Mg24_to_C12_O16)*Y(Mg24)*state.rho + screened_rates(k_He4_N13_to_p_O16)*Y(N13)*state.rho - screened_rates(k_He4_O16_to_Ne20)*Y(O16)*state.rho - screened_rates(k_He4_O16_to_p_F19)*Y(O16)*state.rho;
    jac.set(O16, He4, scratch);

    scratch = -screened_rates(k_C12_O16_to_He4_Mg24)*Y(O16)*state.rho + screened_rates(k_He4_C12_to_O16)*Y(He4)*state.rho;
    jac.set(O16, C12, scratch);

    scratch = screened_rates(k_He4_N13_to_p_O16)*Y(He4)*state.rho;
    jac.set(O16, N13, scratch);

    scratch = screened_rates(k_p_N15_to_O16)*Y(H1)*state.rho;
    jac.set(O16, N15, scratch);

    scratch = -screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*state.rho - screened_rates(k_He4_O16_to_Ne20)*Y(He4)*state.rho - screened_rates(k_He4_O16_to_p_F19)*Y(He4)*state.rho - screened_rates(k_O16_to_He4_C12) - screened_rates(k_O16_to_p_N15) - screened_rates(k_p_O16_to_F17)*Y(H1)*state.rho - screened_rates(k_p_O16_to_He4_N13)*Y(H1)*state.rho;
    jac.set(O16, O16, scratch);

    scratch = screened_rates(k_F17_to_p_O16);
    jac.set(O16, F17, scratch);

    scratch = screened_rates(k_p_F19_to_He4_O16)*Y(H1)*state.rho;
    jac.set(O16, F19, scratch);

    scratch = screened_rates(k_Ne20_to_He4_O16);
    jac.set(O16, Ne20, scratch);

    scratch = screened_rates(k_He4_Mg24_to_C12_O16)*Y(He4)*state.rho;
    jac.set(O16, Mg24, scratch);

    scratch = -screened_rates(k_p_O17_to_F18)*Y(O17)*state.rho - screened_rates(k_p_O17_to_He4_N14)*Y(O17)*state.rho;
    jac.set(O17, H1, scratch);

    scratch = screened_rates(k_He4_N14_to_p_O17)*Y(N14)*state.rho;
    jac.set(O17, He4, scratch);

    scratch = screened_rates(k_He4_N14_to_p_O17)*Y(He4)*state.rho;
    jac.set(O17, N14, scratch);

    scratch = -screened_rates(k_p_O17_to_F18)*Y(H1)*state.rho - screened_rates(k_p_O17_to_He4_N14)*Y(H1)*state.rho;
    jac.set(O17, O17, scratch);

    scratch = screened_rates(k_F17_to_O17_weak_wc12);
    jac.set(O17, F17, scratch);

    scratch = screened_rates(k_F18_to_p_O17);
    jac.set(O17, F18, scratch);

    scratch = -screened_rates(k_p_O18_to_F19)*Y(O18)*state.rho - screened_rates(k_p_O18_to_He4_N15)*Y(O18)*state.rho;
    jac.set(O18, H1, scratch);

    scratch = screened_rates(k_He4_N15_to_p_O18)*Y(N15)*state.rho;
    jac.set(O18, He4, scratch);

    scratch = screened_rates(k_He4_N15_to_p_O18)*Y(He4)*state.rho;
    jac.set(O18, N15, scratch);

    scratch = -screened_rates(k_p_O18_to_F19)*Y(H1)*state.rho - screened_rates(k_p_O18_to_He4_N15)*Y(H1)*state.rho;
    jac.set(O18, O18, scratch);

    scratch = screened_rates(k_F18_to_O18_weak_wc12);
    jac.set(O18, F18, scratch);

    scratch = screened_rates(k_F19_to_p_O18);
    jac.set(O18, F19, scratch);

    scratch = -screened_rates(k_p_F17_to_He4_O14)*Y(F17)*state.rho - screened_rates(k_p_F17_to_Ne18)*Y(F17)*state.rho + screened_rates(k_p_Ne20_to_He4_F17)*Y(Ne20)*state.rho + screened_rates(k_p_O16_to_F17)*Y(O16)*state.rho;
    jac.set(F17, H1, scratch);

    scratch = -screened_rates(k_He4_F17_to_p_Ne20)*Y(F17)*state.rho + screened_rates(k_He4_O14_to_p_F17)*Y(O14)*state.rho;
    jac.set(F17, He4, scratch);

    scratch = screened_rates(k_He4_O14_to_p_F17)*Y(He4)*state.rho;
    jac.set(F17, O14, scratch);

    scratch = screened_rates(k_p_O16_to_F17)*Y(H1)*state.rho;
    jac.set(F17, O16, scratch);

    scratch = -screened_rates(k_F17_to_O17_weak_wc12) - screened_rates(k_F17_to_p_O16) - screened_rates(k_He4_F17_to_p_Ne20)*Y(He4)*state.rho - screened_rates(k_p_F17_to_He4_O14)*Y(H1)*state.rho - screened_rates(k_p_F17_to_Ne18)*Y(H1)*state.rho;
    jac.set(F17, F17, scratch);

    scratch = screened_rates(k_Ne18_to_p_F17);
    jac.set(F17, Ne18, scratch);

    scratch = screened_rates(k_p_Ne20_to_He4_F17)*Y(H1)*state.rho;
    jac.set(F17, Ne20, scratch);

    scratch = -screened_rates(k_p_F18_to_He4_O15)*Y(F18)*state.rho - screened_rates(k_p_F18_to_Ne19)*Y(F18)*state.rho + screened_rates(k_p_O17_to_F18)*Y(O17)*state.rho;
    jac.set(F18, H1, scratch);

    scratch = screened_rates(k_He4_N14_to_F18)*Y(N14)*state.rho + screened_rates(k_He4_O15_to_p_F18)*Y(O15)*state.rho;
    jac.set(F18, He4, scratch);

    scratch = screened_rates(k_He4_N14_to_F18)*Y(He4)*state.rho;
    jac.set(F18, N14, scratch);

    scratch = screened_rates(k_He4_O15_to_p_F18)*Y(He4)*state.rho;
    jac.set(F18, O15, scratch);

    scratch = screened_rates(k_p_O17_to_F18)*Y(H1)*state.rho;
    jac.set(F18, O17, scratch);

    scratch = -screened_rates(k_F18_to_He4_N14) - screened_rates(k_F18_to_O18_weak_wc12) - screened_rates(k_F18_to_p_O17) - screened_rates(k_p_F18_to_He4_O15)*Y(H1)*state.rho - screened_rates(k_p_F18_to_Ne19)*Y(H1)*state.rho;
    jac.set(F18, F18, scratch);

    scratch = screened_rates(k_Ne18_to_F18_weak_wc12);
    jac.set(F18, Ne18, scratch);

    scratch = screened_rates(k_Ne19_to_p_F18);
    jac.set(F18, Ne19, scratch);

    scratch = -screened_rates(k_p_F19_to_He4_O16)*Y(F19)*state.rho - screened_rates(k_p_F19_to_Ne20)*Y(F19)*state.rho + screened_rates(k_p_O18_to_F19)*Y(O18)*state.rho;
    jac.set(F19, H1, scratch);

    scratch = screened_rates(k_He4_N15_to_F19)*Y(N15)*state.rho + screened_rates(k_He4_O16_to_p_F19)*Y(O16)*state.rho;
    jac.set(F19, He4, scratch);

    scratch = screened_rates(k_He4_N15_to_F19)*Y(He4)*state.rho;
    jac.set(F19, N15, scratch);

    scratch = screened_rates(k_He4_O16_to_p_F19)*Y(He4)*state.rho;
    jac.set(F19, O16, scratch);

    scratch = screened_rates(k_p_O18_to_F19)*Y(H1)*state.rho;
    jac.set(F19, O18, scratch);

    scratch = -screened_rates(k_F19_to_He4_N15) - screened_rates(k_F19_to_p_O18) - screened_rates(k_p_F19_to_He4_O16)*Y(H1)*state.rho - screened_rates(k_p_F19_to_Ne20)*Y(H1)*state.rho;
    jac.set(F19, F19, scratch);

    scratch = screened_rates(k_Ne19_to_F19_weak_wc12);
    jac.set(F19, Ne19, scratch);

    scratch = screened_rates(k_Ne20_to_p_F19);
    jac.set(F19, Ne20, scratch);

    scratch = screened_rates(k_p_F17_to_Ne18)*Y(F17)*state.rho;
    jac.set(Ne18, H1, scratch);

    scratch = -screened_rates(k_He4_Ne18_to_Mg22)*Y(Ne18)*state.rho + screened_rates(k_He4_O14_to_Ne18)*Y(O14)*state.rho;
    jac.set(Ne18, He4, scratch);

    scratch = screened_rates(k_He4_O14_to_Ne18)*Y(He4)*state.rho;
    jac.set(Ne18, O14, scratch);

    scratch = screened_rates(k_p_F17_to_Ne18)*Y(H1)*state.rho;
    jac.set(Ne18, F17, scratch);

    scratch = -screened_rates(k_He4_Ne18_to_Mg22)*Y(He4)*state.rho - screened_rates(k_Ne18_to_F18_weak_wc12) - screened_rates(k_Ne18_to_He4_O14) - screened_rates(k_Ne18_to_p_F17);
    jac.set(Ne18, Ne18, scratch);

    scratch = screened_rates(k_Mg22_to_He4_Ne18);
    jac.set(Ne18, Mg22, scratch);

    scratch = screened_rates(k_p_F18_to_Ne19)*Y(F18)*state.rho;
    jac.set(Ne19, H1, scratch);

    scratch = screened_rates(k_He4_O15_to_Ne19)*Y(O15)*state.rho;
    jac.set(Ne19, He4, scratch);

    scratch = screened_rates(k_He4_O15_to_Ne19)*Y(He4)*state.rho;
    jac.set(Ne19, O15, scratch);

    scratch = screened_rates(k_p_F18_to_Ne19)*Y(H1)*state.rho;
    jac.set(Ne19, F18, scratch);

    scratch = -screened_rates(k_Ne19_to_F19_weak_wc12) - screened_rates(k_Ne19_to_He4_O15) - screened_rates(k_Ne19_to_p_F18);
    jac.set(Ne19, Ne19, scratch);

    scratch = screened_rates(k_p_F19_to_Ne20)*Y(F19)*state.rho - screened_rates(k_p_Ne20_to_He4_F17)*Y(Ne20)*state.rho;
    jac.set(Ne20, H1, scratch);

    scratch = screened_rates(k_He4_F17_to_p_Ne20)*Y(F17)*state.rho - screened_rates(k_He4_Ne20_to_C12_C12)*Y(Ne20)*state.rho - screened_rates(k_He4_Ne20_to_Mg24)*Y(Ne20)*state.rho + screened_rates(k_He4_O16_to_Ne20)*Y(O16)*state.rho;
    jac.set(Ne20, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_He4_Ne20)*Y(C12)*state.rho;
    jac.set(Ne20, C12, scratch);

    scratch = screened_rates(k_He4_O16_to_Ne20)*Y(He4)*state.rho;
    jac.set(Ne20, O16, scratch);

    scratch = screened_rates(k_He4_F17_to_p_Ne20)*Y(He4)*state.rho;
    jac.set(Ne20, F17, scratch);

    scratch = screened_rates(k_p_F19_to_Ne20)*Y(H1)*state.rho;
    jac.set(Ne20, F19, scratch);

    scratch = -screened_rates(k_He4_Ne20_to_C12_C12)*Y(He4)*state.rho - screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*state.rho - screened_rates(k_Ne20_to_He4_O16) - screened_rates(k_Ne20_to_p_F19) - screened_rates(k_p_Ne20_to_He4_F17)*Y(H1)*state.rho;
    jac.set(Ne20, Ne20, scratch);

    scratch = screened_rates(k_Mg24_to_He4_Ne20);
    jac.set(Ne20, Mg24, scratch);

    scratch = screened_rates(k_He4_Ne18_to_Mg22)*Y(Ne18)*state.rho;
    jac.set(Mg22, He4, scratch);

    scratch = screened_rates(k_He4_Ne18_to_Mg22)*Y(He4)*state.rho;
    jac.set(Mg22, Ne18, scratch);

    scratch = -screened_rates(k_Mg22_to_He4_Ne18);
    jac.set(Mg22, Mg22, scratch);

    scratch = -screened_rates(k_He4_Mg24_to_C12_O16)*Y(Mg24)*state.rho + screened_rates(k_He4_Ne20_to_Mg24)*Y(Ne20)*state.rho;
    jac.set(Mg24, He4, scratch);

    scratch = screened_rates(k_C12_O16_to_He4_Mg24)*Y(O16)*state.rho;
    jac.set(Mg24, C12, scratch);

    scratch = screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*state.rho;
    jac.set(Mg24, O16, scratch);

    scratch = screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*state.rho;
    jac.set(Mg24, Ne20, scratch);

    scratch = -screened_rates(k_He4_Mg24_to_C12_O16)*Y(He4)*state.rho - screened_rates(k_Mg24_to_He4_Ne20);
    jac.set(Mg24, Mg24, scratch);


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


    jac.zero();

    rate_derivs_t rate_eval;

    constexpr int do_T_derivatives = 1;

    evaluate_rates<do_T_derivatives, rate_derivs_t>(state, rate_eval);

    // Species Jacobian elements with respect to other species

    jac_nuc(state, jac, Y, rate_eval.screened_rates);

    // Energy generation rate Jacobian elements with respect to species

    for (int j = 1; j <= NumSpec; ++j) {
        auto jac_slice_2 = [&](int i) -> amrex::Real { return jac.get(i, j); };
        ener_gener_rate(jac_slice_2, jac(net_ienuc,j));
    }

    // Account for the thermal neutrino losses

    amrex::Real sneut, dsneutdt, dsneutdd, dsnuda, dsnudz;
    constexpr int do_derivatives{1};
    sneut5<do_derivatives>(state.T, state.rho, state.abar, state.zbar, sneut, dsneutdt, dsneutdd, dsnuda, dsnudz);

    for (int j = 1; j <= NumSpec; ++j) {
       amrex::Real b1 = (-state.abar * state.abar * dsnuda + (zion[j-1] - state.zbar) * state.abar * dsnudz);
       jac.add(net_ienuc, j, -b1);
    }


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
