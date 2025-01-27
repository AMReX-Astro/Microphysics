#ifndef actual_rhs_H
#define actual_rhs_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <extern_parameters.H>
#include <actual_network.H>
#include <burn_type.H>
#include <jacobian_utilities.H>
#include <screen.H>
#include <microphysics_autodiff.H>
#ifdef NEUTRINOS
#include <sneut5.H>
#endif
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

    using number_t = std::conditional_t<do_T_derivatives, autodiff::dual, amrex::Real>;
    number_t temp = state.T;
    if constexpr (do_T_derivatives) {
        // seed the dual number for temperature before calculating anything with it
        autodiff::seed(temp);
    }
    plasma_state_t<number_t> pstate{};
    fill_plasma_state(pstate, temp, state.rho, Y);

    tf_t tfactors = evaluate_tfactors(state.T);

    fill_reaclib_rates<do_T_derivatives, T>(tfactors, rate_eval);



    // Evaluate screening factors

    amrex::Real ratraw, dratraw_dT;
    amrex::Real scor, dscor_dt;
    [[maybe_unused]] amrex::Real scor2, dscor2_dt;


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 6.0_rt, 12.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
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


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_C12_to_O16);
    rate_eval.screened_rates(k_He4_C12_to_O16) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_C12_to_O16);
        rate_eval.dscreened_rates_dT(k_He4_C12_to_O16) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 7.0_rt, 14.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_N14_to_F18);
    rate_eval.screened_rates(k_He4_N14_to_F18) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_N14_to_F18);
        rate_eval.dscreened_rates_dT(k_He4_N14_to_F18) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 8.0_rt, 16.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_O16_to_Ne20);
    rate_eval.screened_rates(k_He4_O16_to_Ne20) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_O16_to_Ne20);
        rate_eval.dscreened_rates_dT(k_He4_O16_to_Ne20) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 9.0_rt, 18.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_F18_to_Na22);
    rate_eval.screened_rates(k_He4_F18_to_Na22) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_F18_to_Na22);
        rate_eval.dscreened_rates_dT(k_He4_F18_to_Na22) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_F18_to_p_Ne21);
    rate_eval.screened_rates(k_He4_F18_to_p_Ne21) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_F18_to_p_Ne21);
        rate_eval.dscreened_rates_dT(k_He4_F18_to_p_Ne21) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 10.0_rt, 20.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ne20_to_Mg24);
    rate_eval.screened_rates(k_He4_Ne20_to_Mg24) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ne20_to_Mg24);
        rate_eval.dscreened_rates_dT(k_He4_Ne20_to_Mg24) = ratraw * dscor_dt + dratraw_dT * scor;
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


    ratraw = rate_eval.screened_rates(k_p_Ne21_to_Na22);
    rate_eval.screened_rates(k_p_Ne21_to_Na22) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ne21_to_Na22);
        rate_eval.dscreened_rates_dT(k_p_Ne21_to_Na22) = ratraw * dscor_dt + dratraw_dT * scor;
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


    ratraw = rate_eval.screened_rates(k_p_Na23_to_Mg24);
    rate_eval.screened_rates(k_p_Na23_to_Mg24) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Na23_to_Mg24);
        rate_eval.dscreened_rates_dT(k_p_Na23_to_Mg24) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Na23_to_He4_Ne20);
    rate_eval.screened_rates(k_p_Na23_to_He4_Ne20) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Na23_to_He4_Ne20);
        rate_eval.dscreened_rates_dT(k_p_Na23_to_He4_Ne20) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 12.0_rt, 24.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Mg24_to_Si28);
    rate_eval.screened_rates(k_He4_Mg24_to_Si28) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mg24_to_Si28);
        rate_eval.dscreened_rates_dT(k_He4_Mg24_to_Si28) = ratraw * dscor_dt + dratraw_dT * scor;
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


    ratraw = rate_eval.screened_rates(k_p_Al27_to_Si28);
    rate_eval.screened_rates(k_p_Al27_to_Si28) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Al27_to_Si28);
        rate_eval.dscreened_rates_dT(k_p_Al27_to_Si28) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Al27_to_He4_Mg24);
    rate_eval.screened_rates(k_p_Al27_to_He4_Mg24) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Al27_to_He4_Mg24);
        rate_eval.dscreened_rates_dT(k_p_Al27_to_He4_Mg24) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 14.0_rt, 28.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Si28_to_S32);
    rate_eval.screened_rates(k_He4_Si28_to_S32) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Si28_to_S32);
        rate_eval.dscreened_rates_dT(k_He4_Si28_to_S32) = ratraw * dscor_dt + dratraw_dT * scor;
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


    ratraw = rate_eval.screened_rates(k_p_P31_to_S32);
    rate_eval.screened_rates(k_p_P31_to_S32) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_P31_to_S32);
        rate_eval.dscreened_rates_dT(k_p_P31_to_S32) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_P31_to_He4_Si28);
    rate_eval.screened_rates(k_p_P31_to_He4_Si28) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_P31_to_He4_Si28);
        rate_eval.dscreened_rates_dT(k_p_P31_to_He4_Si28) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 24.0_rt, 48.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Cr48_to_Fe52);
    rate_eval.screened_rates(k_He4_Cr48_to_Fe52) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cr48_to_Fe52);
        rate_eval.dscreened_rates_dT(k_He4_Cr48_to_Fe52) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cr48_to_p_nse_Mn51);
    rate_eval.screened_rates(k_He4_Cr48_to_p_nse_Mn51) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cr48_to_p_nse_Mn51);
        rate_eval.dscreened_rates_dT(k_He4_Cr48_to_p_nse_Mn51) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 25.0_rt, 51.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_nse_Mn51_to_Fe52);
    rate_eval.screened_rates(k_p_nse_Mn51_to_Fe52) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_nse_Mn51_to_Fe52);
        rate_eval.dscreened_rates_dT(k_p_nse_Mn51_to_Fe52) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_nse_Mn51_to_He4_Cr48_derived);
    rate_eval.screened_rates(k_p_nse_Mn51_to_He4_Cr48_derived) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_nse_Mn51_to_He4_Cr48_derived);
        rate_eval.dscreened_rates_dT(k_p_nse_Mn51_to_He4_Cr48_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 25.0_rt, 51.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Mn51_to_Co55);
    rate_eval.screened_rates(k_He4_Mn51_to_Co55) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mn51_to_Co55);
        rate_eval.dscreened_rates_dT(k_He4_Mn51_to_Co55) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Mn51_to_p_nse_Fe54);
    rate_eval.screened_rates(k_He4_Mn51_to_p_nse_Fe54) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mn51_to_p_nse_Fe54);
        rate_eval.dscreened_rates_dT(k_He4_Mn51_to_p_nse_Fe54) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 26.0_rt, 52.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Fe52_to_Ni56);
    rate_eval.screened_rates(k_He4_Fe52_to_Ni56) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe52_to_Ni56);
        rate_eval.dscreened_rates_dT(k_He4_Fe52_to_Ni56) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Fe52_to_p_nse_Co55);
    rate_eval.screened_rates(k_He4_Fe52_to_p_nse_Co55) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe52_to_p_nse_Co55);
        rate_eval.dscreened_rates_dT(k_He4_Fe52_to_p_nse_Co55) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 27.0_rt, 55.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_nse_Co55_to_Ni56);
    rate_eval.screened_rates(k_p_nse_Co55_to_Ni56) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_nse_Co55_to_Ni56);
        rate_eval.dscreened_rates_dT(k_p_nse_Co55_to_Ni56) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_nse_Co55_to_He4_Fe52_derived);
    rate_eval.screened_rates(k_p_nse_Co55_to_He4_Fe52_derived) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_nse_Co55_to_He4_Fe52_derived);
        rate_eval.dscreened_rates_dT(k_p_nse_Co55_to_He4_Fe52_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(6.0_rt, 12.0_rt, 6.0_rt, 12.0_rt);


        static_assert(scn_fac.z1 == 6.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_C12_C12_to_p_Na23);
    rate_eval.screened_rates(k_C12_C12_to_p_Na23) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_C12_to_p_Na23);
        rate_eval.dscreened_rates_dT(k_C12_C12_to_p_Na23) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_C12_C12_to_He4_Ne20);
    rate_eval.screened_rates(k_C12_C12_to_He4_Ne20) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_C12_to_He4_Ne20);
        rate_eval.dscreened_rates_dT(k_C12_C12_to_He4_Ne20) = ratraw * dscor_dt + dratraw_dT * scor;
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


    ratraw = rate_eval.screened_rates(k_He4_N13_to_p_O16);
    rate_eval.screened_rates(k_He4_N13_to_p_O16) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_N13_to_p_O16);
        rate_eval.dscreened_rates_dT(k_He4_N13_to_p_O16) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(6.0_rt, 12.0_rt, 8.0_rt, 16.0_rt);


        static_assert(scn_fac.z1 == 6.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_C12_O16_to_p_Al27);
    rate_eval.screened_rates(k_C12_O16_to_p_Al27) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_O16_to_p_Al27);
        rate_eval.dscreened_rates_dT(k_C12_O16_to_p_Al27) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_C12_O16_to_He4_Mg24);
    rate_eval.screened_rates(k_C12_O16_to_He4_Mg24) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_O16_to_He4_Mg24);
        rate_eval.dscreened_rates_dT(k_C12_O16_to_He4_Mg24) = ratraw * dscor_dt + dratraw_dT * scor;
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


    ratraw = rate_eval.screened_rates(k_O16_O16_to_p_P31);
    rate_eval.screened_rates(k_O16_O16_to_p_P31) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_O16_O16_to_p_P31);
        rate_eval.dscreened_rates_dT(k_O16_O16_to_p_P31) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_O16_O16_to_He4_Si28);
    rate_eval.screened_rates(k_O16_O16_to_He4_Si28) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_O16_O16_to_He4_Si28);
        rate_eval.dscreened_rates_dT(k_O16_O16_to_He4_Si28) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_O16_O16_to_S32_modified);
    rate_eval.screened_rates(k_O16_O16_to_S32_modified) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_O16_O16_to_S32_modified);
        rate_eval.dscreened_rates_dT(k_O16_O16_to_S32_modified) = ratraw * dscor_dt + dratraw_dT * scor;
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


    ratraw = rate_eval.screened_rates(k_He4_He4_He4_to_C12);
    rate_eval.screened_rates(k_He4_He4_He4_to_C12) *= scor * scor2;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_He4_He4_to_C12);
        rate_eval.dscreened_rates_dT(k_He4_He4_He4_to_C12) = ratraw * (scor * dscor2_dt + dscor_dt * scor2) + dratraw_dT * scor * scor2;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 26.0_rt, 54.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_nse_Fe54_to_Co55);
    rate_eval.screened_rates(k_p_nse_Fe54_to_Co55) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_nse_Fe54_to_Co55);
        rate_eval.dscreened_rates_dT(k_p_nse_Fe54_to_Co55) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_nse_Fe54_to_He4_Mn51_derived);
    rate_eval.screened_rates(k_p_nse_Fe54_to_He4_Mn51_derived) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_nse_Fe54_to_He4_Mn51_derived);
        rate_eval.dscreened_rates_dT(k_p_nse_Fe54_to_He4_Mn51_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 26.0_rt, 54.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Fe54_to_Ni58);
    rate_eval.screened_rates(k_He4_Fe54_to_Ni58) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe54_to_Ni58);
        rate_eval.dscreened_rates_dT(k_He4_Fe54_to_Ni58) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Fe54_to_p_nse_Co57_derived);
    rate_eval.screened_rates(k_He4_Fe54_to_p_nse_Co57_derived) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe54_to_p_nse_Co57_derived);
        rate_eval.dscreened_rates_dT(k_He4_Fe54_to_p_nse_Co57_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 26.0_rt, 56.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_nse_Fe56_to_Co57);
    rate_eval.screened_rates(k_p_nse_Fe56_to_Co57) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_nse_Fe56_to_Co57);
        rate_eval.dscreened_rates_dT(k_p_nse_Fe56_to_Co57) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_nse_Fe56_to_n_Co56_derived);
    rate_eval.screened_rates(k_p_nse_Fe56_to_n_Co56_derived) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_nse_Fe56_to_n_Co56_derived);
        rate_eval.dscreened_rates_dT(k_p_nse_Fe56_to_n_Co56_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 27.0_rt, 57.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_nse_Co57_to_Ni58);
    rate_eval.screened_rates(k_p_nse_Co57_to_Ni58) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_nse_Co57_to_Ni58);
        rate_eval.dscreened_rates_dT(k_p_nse_Co57_to_Ni58) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_nse_Co57_to_He4_Fe54);
    rate_eval.screened_rates(k_p_nse_Co57_to_He4_Fe54) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_nse_Co57_to_He4_Fe54);
        rate_eval.dscreened_rates_dT(k_p_nse_Co57_to_He4_Fe54) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 27.0_rt, 55.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Co55_to_p_nse_Ni58);
    rate_eval.screened_rates(k_He4_Co55_to_p_nse_Ni58) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Co55_to_p_nse_Ni58);
        rate_eval.dscreened_rates_dT(k_He4_Co55_to_p_nse_Ni58) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 8.0_rt, 16.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_O16_to_He4_N13_derived);
    rate_eval.screened_rates(k_p_O16_to_He4_N13_derived) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_O16_to_He4_N13_derived);
        rate_eval.dscreened_rates_dT(k_p_O16_to_He4_N13_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 27.0_rt, 56.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_nse_Co56_to_n_Ni56_derived);
    rate_eval.screened_rates(k_p_nse_Co56_to_n_Ni56_derived) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_nse_Co56_to_n_Ni56_derived);
        rate_eval.dscreened_rates_dT(k_p_nse_Co56_to_n_Ni56_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 28.0_rt, 58.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_nse_Ni58_to_He4_Co55_derived);
    rate_eval.screened_rates(k_p_nse_Ni58_to_He4_Co55_derived) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_nse_Ni58_to_He4_Co55_derived);
        rate_eval.dscreened_rates_dT(k_p_nse_Ni58_to_He4_Co55_derived) = ratraw * dscor_dt + dratraw_dT * scor;
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

    ratraw = rate_eval.screened_rates(k_He4_S32_to_p_Cl35_derived_removed);
    rate_eval.screened_rates(k_He4_S32_to_p_Cl35_derived_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_S32_to_p_Cl35_derived_removed);
        rate_eval.dscreened_rates_dT(k_He4_S32_to_p_Cl35_derived_removed) = ratraw * dscor_dt + dratraw_dT * scor;
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

    ratraw = rate_eval.screened_rates(k_He4_Ar36_to_p_K39_derived_removed);
    rate_eval.screened_rates(k_He4_Ar36_to_p_K39_derived_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ar36_to_p_K39_derived_removed);
        rate_eval.dscreened_rates_dT(k_He4_Ar36_to_p_K39_derived_removed) = ratraw * dscor_dt + dratraw_dT * scor;
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

    ratraw = rate_eval.screened_rates(k_He4_Ca40_to_p_Sc43_derived_removed);
    rate_eval.screened_rates(k_He4_Ca40_to_p_Sc43_derived_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca40_to_p_Sc43_derived_removed);
        rate_eval.dscreened_rates_dT(k_He4_Ca40_to_p_Sc43_derived_removed) = ratraw * dscor_dt + dratraw_dT * scor;
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

    ratraw = rate_eval.screened_rates(k_p_V47_to_He4_Ti44_derived_removed);
    rate_eval.screened_rates(k_p_V47_to_He4_Ti44_derived_removed) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_V47_to_He4_Ti44_derived_removed);
        rate_eval.dscreened_rates_dT(k_p_V47_to_He4_Ti44_derived_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    // Fill approximate rates

    fill_approx_rates<do_T_derivatives, T>(tfactors, state.rho, Y, rate_eval);

    // Calculate tabular rates

    [[maybe_unused]] amrex::Real rate, drate_dt, edot_nu, edot_gamma;

    rate_eval.enuc_weak = 0.0_rt;

    tabular_evaluate(j_Co56_Fe56_meta, j_Co56_Fe56_rhoy, j_Co56_Fe56_temp, j_Co56_Fe56_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Co56_to_Fe56) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Co56_to_Fe56) = drate_dt;
    }
    rate_eval.enuc_weak += C::Legacy::n_A * Y(Co56) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Co56_Ni56_meta, j_Co56_Ni56_rhoy, j_Co56_Ni56_temp, j_Co56_Ni56_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Co56_to_Ni56) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Co56_to_Ni56) = drate_dt;
    }
    rate_eval.enuc_weak += C::Legacy::n_A * Y(Co56) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Fe56_Co56_meta, j_Fe56_Co56_rhoy, j_Fe56_Co56_temp, j_Fe56_Co56_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Fe56_to_Co56) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe56_to_Co56) = drate_dt;
    }
    rate_eval.enuc_weak += C::Legacy::n_A * Y(Fe56) * (edot_nu + edot_gamma);

    tabular_evaluate(j_n_p_meta, j_n_p_rhoy, j_n_p_temp, j_n_p_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_n_to_p) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_to_p) = drate_dt;
    }
    rate_eval.enuc_weak += C::Legacy::n_A * Y(N) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Ni56_Co56_meta, j_Ni56_Co56_rhoy, j_Ni56_Co56_temp, j_Ni56_Co56_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Ni56_to_Co56) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni56_to_Co56) = drate_dt;
    }
    rate_eval.enuc_weak += C::Legacy::n_A * Y(Ni56) * (edot_nu + edot_gamma);

    tabular_evaluate(j_p_n_meta, j_p_n_rhoy, j_p_n_temp, j_p_n_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_p_to_n) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_to_n) = drate_dt;
    }
    rate_eval.enuc_weak += C::Legacy::n_A * Y(H1) * (edot_nu + edot_gamma);


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

    tabular_evaluate(j_Co56_Fe56_meta, j_Co56_Fe56_rhoy, j_Co56_Fe56_temp, j_Co56_Fe56_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Co56_to_Fe56) = rate;
    rate_eval.enuc_weak += C::Legacy::n_A * Y(Co56) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Co56_Ni56_meta, j_Co56_Ni56_rhoy, j_Co56_Ni56_temp, j_Co56_Ni56_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Co56_to_Ni56) = rate;
    rate_eval.enuc_weak += C::Legacy::n_A * Y(Co56) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Fe56_Co56_meta, j_Fe56_Co56_rhoy, j_Fe56_Co56_temp, j_Fe56_Co56_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Fe56_to_Co56) = rate;
    rate_eval.enuc_weak += C::Legacy::n_A * Y(Fe56) * (edot_nu + edot_gamma);

    tabular_evaluate(j_n_p_meta, j_n_p_rhoy, j_n_p_temp, j_n_p_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_n_to_p) = rate;
    rate_eval.enuc_weak += C::Legacy::n_A * Y(N) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Ni56_Co56_meta, j_Ni56_Co56_rhoy, j_Ni56_Co56_temp, j_Ni56_Co56_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Ni56_to_Co56) = rate;
    rate_eval.enuc_weak += C::Legacy::n_A * Y(Ni56) * (edot_nu + edot_gamma);

    tabular_evaluate(j_p_n_meta, j_p_n_rhoy, j_p_n_temp, j_p_n_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_p_to_n) = rate;
    rate_eval.enuc_weak += C::Legacy::n_A * Y(H1) * (edot_nu + edot_gamma);

    auto screened_rates = rate_eval.screened_rates;

    ydot_nuc(N) =
        (-screened_rates(k_n_to_p)*Y(N) + screened_rates(k_p_to_n)*Y(H1));

    ydot_nuc(H1) =
        (screened_rates(k_n_to_p)*Y(N) + -screened_rates(k_p_to_n)*Y(H1));

    ydot_nuc(P_nse) = 0.0_rt;

    ydot_nuc(He4) = 0.0_rt;

    ydot_nuc(C12) = 0.0_rt;

    ydot_nuc(N13) = 0.0_rt;

    ydot_nuc(N14) = 0.0_rt;

    ydot_nuc(O16) = 0.0_rt;

    ydot_nuc(F18) = 0.0_rt;

    ydot_nuc(Ne20) = 0.0_rt;

    ydot_nuc(Ne21) = 0.0_rt;

    ydot_nuc(Na22) = 0.0_rt;

    ydot_nuc(Na23) = 0.0_rt;

    ydot_nuc(Mg24) = 0.0_rt;

    ydot_nuc(Al27) = 0.0_rt;

    ydot_nuc(Si28) = 0.0_rt;

    ydot_nuc(P31) = 0.0_rt;

    ydot_nuc(S32) = 0.0_rt;

    ydot_nuc(Ar36) = 0.0_rt;

    ydot_nuc(Ca40) = 0.0_rt;

    ydot_nuc(Ti44) = 0.0_rt;

    ydot_nuc(Cr48) = 0.0_rt;

    ydot_nuc(Mn51) = 0.0_rt;

    ydot_nuc(Fe52) = 0.0_rt;

    ydot_nuc(Fe54) = 0.0_rt;

    ydot_nuc(Fe56) =
        (screened_rates(k_Co56_to_Fe56)*Y(Co56) + -screened_rates(k_Fe56_to_Co56)*Y(Fe56));

    ydot_nuc(Co55) = 0.0_rt;

    ydot_nuc(Co56) =
        (-screened_rates(k_Co56_to_Fe56)*Y(Co56) + screened_rates(k_Fe56_to_Co56)*Y(Fe56)) +
        (screened_rates(k_Ni56_to_Co56)*Y(Ni56) + -screened_rates(k_Co56_to_Ni56)*Y(Co56));

    ydot_nuc(Co57) = 0.0_rt;

    ydot_nuc(Ni56) =
        (-screened_rates(k_Ni56_to_Co56)*Y(Ni56) + screened_rates(k_Co56_to_Ni56)*Y(Co56));

    ydot_nuc(Ni58) = 0.0_rt;

    enuc_weak = rate_eval.enuc_weak;
}
#endif


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rhs_nuc(const burn_t& state,
             amrex::Array1D<amrex::Real, 1, neqs>& ydot_nuc,
             const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
             const amrex::Array1D<amrex::Real, 1, NumRates>& screened_rates) {

    using namespace Rates;

    ydot_nuc(N) =
        (-screened_rates(k_n_Co55_to_Co56)*Y(Co55)*Y(N)*state.rho + screened_rates(k_Co56_to_n_Co55_derived)*Y(Co56)) +
        (-screened_rates(k_n_Co56_to_Co57)*Y(Co56)*Y(N)*state.rho + screened_rates(k_Co57_to_n_Co56_derived)*Y(Co57)) +
        (-screened_rates(k_n_Co56_to_p_nse_Fe56)*Y(Co56)*Y(N)*state.rho + screened_rates(k_p_nse_Fe56_to_n_Co56_derived)*Y(Fe56)*Y(P_nse)*state.rho) +
        (-screened_rates(k_n_Ni56_to_p_nse_Co56)*Y(Ni56)*Y(N)*state.rho + screened_rates(k_p_nse_Co56_to_n_Ni56_derived)*Y(Co56)*Y(P_nse)*state.rho) +
        (-2.0*screened_rates(k_Fe52_n_n_to_Fe54_approx)*Y(Fe52)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + 2.0*screened_rates(k_Fe54_to_Fe52_n_n_approx)*Y(Fe54)) +
        (-2.0*screened_rates(k_Fe54_n_n_to_Fe56_approx)*Y(Fe54)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + 2.0*screened_rates(k_Fe56_to_Fe54_n_n_approx)*Y(Fe56)) +
        (-2.0*screened_rates(k_Ni56_n_n_to_Ni58_approx)*Y(Ni56)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + 2.0*screened_rates(k_Ni58_to_Ni56_n_n_approx)*Y(Ni58)) +
        (-screened_rates(k_n_to_p)*Y(N) + screened_rates(k_p_to_n)*Y(H1));

    ydot_nuc(H1) =
        (-screened_rates(k_p_C12_to_N13)*Y(C12)*Y(H1)*state.rho + screened_rates(k_N13_to_p_C12_derived)*Y(N13)) +
        (-screened_rates(k_p_Ne21_to_Na22)*Y(Ne21)*Y(H1)*state.rho + screened_rates(k_Na22_to_p_Ne21_derived)*Y(Na22)) +
        (-screened_rates(k_p_Na23_to_Mg24)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_Mg24_to_p_Na23_derived)*Y(Mg24)) +
        (-screened_rates(k_p_Al27_to_Si28)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_Si28_to_p_Al27_derived)*Y(Si28)) +
        (-screened_rates(k_p_P31_to_S32)*Y(P31)*Y(H1)*state.rho + screened_rates(k_S32_to_p_P31_derived)*Y(S32)) +
        0.5*screened_rates(k_C12_C12_to_p_Na23)*amrex::Math::powi<2>(Y(C12))*state.rho +
        (screened_rates(k_He4_N13_to_p_O16)*Y(He4)*Y(N13)*state.rho + -screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*Y(H1)*state.rho) +
        screened_rates(k_C12_O16_to_p_Al27)*Y(C12)*Y(O16)*state.rho +
        0.5*screened_rates(k_O16_O16_to_p_P31)*amrex::Math::powi<2>(Y(O16))*state.rho +
        (screened_rates(k_He4_F18_to_p_Ne21)*Y(F18)*Y(He4)*state.rho + -screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(Ne21)*Y(H1)*state.rho) +
        (-screened_rates(k_p_Na23_to_He4_Ne20)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(He4)*Y(Ne20)*state.rho) +
        (-screened_rates(k_p_Al27_to_He4_Mg24)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(He4)*Y(Mg24)*state.rho) +
        (-screened_rates(k_p_P31_to_He4_Si28)*Y(P31)*Y(H1)*state.rho + screened_rates(k_He4_Si28_to_p_P31_derived)*Y(He4)*Y(Si28)*state.rho) +
        (screened_rates(k_n_to_p)*Y(N) + -screened_rates(k_p_to_n)*Y(H1));

    ydot_nuc(P_nse) =
        (-screened_rates(k_p_nse_Mn51_to_Fe52)*Y(Mn51)*Y(P_nse)*state.rho + screened_rates(k_Fe52_to_p_nse_Mn51_derived)*Y(Fe52)) +
        (-screened_rates(k_p_nse_Co55_to_Ni56)*Y(Co55)*Y(P_nse)*state.rho + screened_rates(k_Ni56_to_p_nse_Co55_derived)*Y(Ni56)) +
        (screened_rates(k_He4_Cr48_to_p_nse_Mn51)*Y(Cr48)*Y(He4)*state.rho + -screened_rates(k_p_nse_Mn51_to_He4_Cr48_derived)*Y(Mn51)*Y(P_nse)*state.rho) +
        (screened_rates(k_He4_Fe52_to_p_nse_Co55)*Y(Fe52)*Y(He4)*state.rho + -screened_rates(k_p_nse_Co55_to_He4_Fe52_derived)*Y(Co55)*Y(P_nse)*state.rho) +
        (-screened_rates(k_p_nse_Fe54_to_Co55)*Y(Fe54)*Y(P_nse)*state.rho + screened_rates(k_Co55_to_p_nse_Fe54_derived)*Y(Co55)) +
        (-screened_rates(k_p_nse_Fe56_to_Co57)*Y(Fe56)*Y(P_nse)*state.rho + screened_rates(k_Co57_to_p_nse_Fe56_derived)*Y(Co57)) +
        (-screened_rates(k_p_nse_Co57_to_Ni58)*Y(Co57)*Y(P_nse)*state.rho + screened_rates(k_Ni58_to_p_nse_Co57_derived)*Y(Ni58)) +
        (screened_rates(k_He4_Mn51_to_p_nse_Fe54)*Y(He4)*Y(Mn51)*state.rho + -screened_rates(k_p_nse_Fe54_to_He4_Mn51_derived)*Y(Fe54)*Y(P_nse)*state.rho) +
        (screened_rates(k_He4_Co55_to_p_nse_Ni58)*Y(Co55)*Y(He4)*state.rho + -screened_rates(k_p_nse_Ni58_to_He4_Co55_derived)*Y(Ni58)*Y(P_nse)*state.rho) +
        (screened_rates(k_n_Co56_to_p_nse_Fe56)*Y(Co56)*Y(N)*state.rho + -screened_rates(k_p_nse_Fe56_to_n_Co56_derived)*Y(Fe56)*Y(P_nse)*state.rho) +
        (-screened_rates(k_p_nse_Co57_to_He4_Fe54)*Y(Co57)*Y(P_nse)*state.rho + screened_rates(k_He4_Fe54_to_p_nse_Co57_derived)*Y(Fe54)*Y(He4)*state.rho) +
        (screened_rates(k_n_Ni56_to_p_nse_Co56)*Y(Ni56)*Y(N)*state.rho + -screened_rates(k_p_nse_Co56_to_n_Ni56_derived)*Y(Co56)*Y(P_nse)*state.rho);

    ydot_nuc(He4) =
        (-screened_rates(k_He4_C12_to_O16)*Y(C12)*Y(He4)*state.rho + screened_rates(k_O16_to_He4_C12_derived)*Y(O16)) +
        (-screened_rates(k_He4_N14_to_F18)*Y(He4)*Y(N14)*state.rho + screened_rates(k_F18_to_He4_N14_derived)*Y(F18)) +
        (-screened_rates(k_He4_O16_to_Ne20)*Y(He4)*Y(O16)*state.rho + screened_rates(k_Ne20_to_He4_O16_derived)*Y(Ne20)) +
        (-screened_rates(k_He4_F18_to_Na22)*Y(F18)*Y(He4)*state.rho + screened_rates(k_Na22_to_He4_F18_derived)*Y(Na22)) +
        (-screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*Y(Ne20)*state.rho + screened_rates(k_Mg24_to_He4_Ne20_derived)*Y(Mg24)) +
        (-screened_rates(k_He4_Mg24_to_Si28)*Y(He4)*Y(Mg24)*state.rho + screened_rates(k_Si28_to_He4_Mg24_derived)*Y(Si28)) +
        (-screened_rates(k_He4_Si28_to_S32)*Y(He4)*Y(Si28)*state.rho + screened_rates(k_S32_to_He4_Si28_derived)*Y(S32)) +
        (-screened_rates(k_He4_Cr48_to_Fe52)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_Fe52_to_He4_Cr48_derived)*Y(Fe52)) +
        (-screened_rates(k_He4_Mn51_to_Co55)*Y(He4)*Y(Mn51)*state.rho + screened_rates(k_Co55_to_He4_Mn51_derived)*Y(Co55)) +
        (-screened_rates(k_He4_Fe52_to_Ni56)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_Ni56_to_He4_Fe52_derived)*Y(Ni56)) +
        0.5*screened_rates(k_C12_C12_to_He4_Ne20)*amrex::Math::powi<2>(Y(C12))*state.rho +
        (-screened_rates(k_He4_N13_to_p_O16)*Y(He4)*Y(N13)*state.rho + screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*Y(H1)*state.rho) +
        screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*Y(O16)*state.rho +
        0.5*screened_rates(k_O16_O16_to_He4_Si28)*amrex::Math::powi<2>(Y(O16))*state.rho +
        (-screened_rates(k_He4_F18_to_p_Ne21)*Y(F18)*Y(He4)*state.rho + screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(Ne21)*Y(H1)*state.rho) +
        (screened_rates(k_p_Na23_to_He4_Ne20)*Y(Na23)*Y(H1)*state.rho + -screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(He4)*Y(Ne20)*state.rho) +
        (screened_rates(k_p_Al27_to_He4_Mg24)*Y(Al27)*Y(H1)*state.rho + -screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(He4)*Y(Mg24)*state.rho) +
        (screened_rates(k_p_P31_to_He4_Si28)*Y(P31)*Y(H1)*state.rho + -screened_rates(k_He4_Si28_to_p_P31_derived)*Y(He4)*Y(Si28)*state.rho) +
        (-screened_rates(k_He4_Cr48_to_p_nse_Mn51)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_p_nse_Mn51_to_He4_Cr48_derived)*Y(Mn51)*Y(P_nse)*state.rho) +
        (-screened_rates(k_He4_Fe52_to_p_nse_Co55)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_p_nse_Co55_to_He4_Fe52_derived)*Y(Co55)*Y(P_nse)*state.rho) +
        (-0.5*screened_rates(k_He4_He4_He4_to_C12)*amrex::Math::powi<3>(Y(He4))*amrex::Math::powi<2>(state.rho) + 3.0*screened_rates(k_C12_to_He4_He4_He4_derived)*Y(C12)) +
        (-screened_rates(k_He4_Fe54_to_Ni58)*Y(Fe54)*Y(He4)*state.rho + screened_rates(k_Ni58_to_He4_Fe54_derived)*Y(Ni58)) +
        (-screened_rates(k_He4_Mn51_to_p_nse_Fe54)*Y(He4)*Y(Mn51)*state.rho + screened_rates(k_p_nse_Fe54_to_He4_Mn51_derived)*Y(Fe54)*Y(P_nse)*state.rho) +
        (-screened_rates(k_He4_Co55_to_p_nse_Ni58)*Y(Co55)*Y(He4)*state.rho + screened_rates(k_p_nse_Ni58_to_He4_Co55_derived)*Y(Ni58)*Y(P_nse)*state.rho) +
        (screened_rates(k_p_nse_Co57_to_He4_Fe54)*Y(Co57)*Y(P_nse)*state.rho + -screened_rates(k_He4_Fe54_to_p_nse_Co57_derived)*Y(Fe54)*Y(He4)*state.rho) +
        (-screened_rates(k_S32_He4_to_Ar36_approx)*Y(He4)*Y(S32)*state.rho + screened_rates(k_Ar36_to_S32_He4_approx)*Y(Ar36)) +
        (-screened_rates(k_Ar36_He4_to_Ca40_approx)*Y(Ar36)*Y(He4)*state.rho + screened_rates(k_Ca40_to_Ar36_He4_approx)*Y(Ca40)) +
        (-screened_rates(k_Ca40_He4_to_Ti44_approx)*Y(Ca40)*Y(He4)*state.rho + screened_rates(k_Ti44_to_Ca40_He4_approx)*Y(Ti44)) +
        (-screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(He4)*Y(Ti44)*state.rho + screened_rates(k_Cr48_to_Ti44_He4_approx)*Y(Cr48));

    ydot_nuc(C12) =
        (-screened_rates(k_p_C12_to_N13)*Y(C12)*Y(H1)*state.rho + screened_rates(k_N13_to_p_C12_derived)*Y(N13)) +
        (-screened_rates(k_He4_C12_to_O16)*Y(C12)*Y(He4)*state.rho + screened_rates(k_O16_to_He4_C12_derived)*Y(O16)) +
        -screened_rates(k_C12_C12_to_p_Na23)*amrex::Math::powi<2>(Y(C12))*state.rho +
        -screened_rates(k_C12_C12_to_He4_Ne20)*amrex::Math::powi<2>(Y(C12))*state.rho +
        -screened_rates(k_C12_O16_to_p_Al27)*Y(C12)*Y(O16)*state.rho +
        -screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*Y(O16)*state.rho +
        (0.16666666666666667*screened_rates(k_He4_He4_He4_to_C12)*amrex::Math::powi<3>(Y(He4))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_C12_to_He4_He4_He4_derived)*Y(C12)) +
        -screened_rates(k_C12_C12_to_Mg24_modified)*amrex::Math::powi<2>(Y(C12))*state.rho +
        -screened_rates(k_C12_O16_to_Si28_modified)*Y(C12)*Y(O16)*state.rho;

    ydot_nuc(N13) =
        (screened_rates(k_p_C12_to_N13)*Y(C12)*Y(H1)*state.rho + -screened_rates(k_N13_to_p_C12_derived)*Y(N13)) +
        (-screened_rates(k_He4_N13_to_p_O16)*Y(He4)*Y(N13)*state.rho + screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*Y(H1)*state.rho);

    ydot_nuc(N14) =
        (-screened_rates(k_He4_N14_to_F18)*Y(He4)*Y(N14)*state.rho + screened_rates(k_F18_to_He4_N14_derived)*Y(F18));

    ydot_nuc(O16) =
        (screened_rates(k_He4_C12_to_O16)*Y(C12)*Y(He4)*state.rho + -screened_rates(k_O16_to_He4_C12_derived)*Y(O16)) +
        (-screened_rates(k_He4_O16_to_Ne20)*Y(He4)*Y(O16)*state.rho + screened_rates(k_Ne20_to_He4_O16_derived)*Y(Ne20)) +
        (screened_rates(k_He4_N13_to_p_O16)*Y(He4)*Y(N13)*state.rho + -screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*Y(H1)*state.rho) +
        -screened_rates(k_C12_O16_to_p_Al27)*Y(C12)*Y(O16)*state.rho +
        -screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*Y(O16)*state.rho +
        -screened_rates(k_O16_O16_to_p_P31)*amrex::Math::powi<2>(Y(O16))*state.rho +
        -screened_rates(k_O16_O16_to_He4_Si28)*amrex::Math::powi<2>(Y(O16))*state.rho +
        -screened_rates(k_O16_O16_to_S32_modified)*amrex::Math::powi<2>(Y(O16))*state.rho +
        -screened_rates(k_C12_O16_to_Si28_modified)*Y(C12)*Y(O16)*state.rho;

    ydot_nuc(F18) =
        (screened_rates(k_He4_N14_to_F18)*Y(He4)*Y(N14)*state.rho + -screened_rates(k_F18_to_He4_N14_derived)*Y(F18)) +
        (-screened_rates(k_He4_F18_to_Na22)*Y(F18)*Y(He4)*state.rho + screened_rates(k_Na22_to_He4_F18_derived)*Y(Na22)) +
        (-screened_rates(k_He4_F18_to_p_Ne21)*Y(F18)*Y(He4)*state.rho + screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(Ne21)*Y(H1)*state.rho);

    ydot_nuc(Ne20) =
        (screened_rates(k_He4_O16_to_Ne20)*Y(He4)*Y(O16)*state.rho + -screened_rates(k_Ne20_to_He4_O16_derived)*Y(Ne20)) +
        (-screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*Y(Ne20)*state.rho + screened_rates(k_Mg24_to_He4_Ne20_derived)*Y(Mg24)) +
        0.5*screened_rates(k_C12_C12_to_He4_Ne20)*amrex::Math::powi<2>(Y(C12))*state.rho +
        (screened_rates(k_p_Na23_to_He4_Ne20)*Y(Na23)*Y(H1)*state.rho + -screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(He4)*Y(Ne20)*state.rho);

    ydot_nuc(Ne21) =
        (-screened_rates(k_p_Ne21_to_Na22)*Y(Ne21)*Y(H1)*state.rho + screened_rates(k_Na22_to_p_Ne21_derived)*Y(Na22)) +
        (screened_rates(k_He4_F18_to_p_Ne21)*Y(F18)*Y(He4)*state.rho + -screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(Ne21)*Y(H1)*state.rho);

    ydot_nuc(Na22) =
        (screened_rates(k_He4_F18_to_Na22)*Y(F18)*Y(He4)*state.rho + -screened_rates(k_Na22_to_He4_F18_derived)*Y(Na22)) +
        (screened_rates(k_p_Ne21_to_Na22)*Y(Ne21)*Y(H1)*state.rho + -screened_rates(k_Na22_to_p_Ne21_derived)*Y(Na22));

    ydot_nuc(Na23) =
        (-screened_rates(k_p_Na23_to_Mg24)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_Mg24_to_p_Na23_derived)*Y(Mg24)) +
        0.5*screened_rates(k_C12_C12_to_p_Na23)*amrex::Math::powi<2>(Y(C12))*state.rho +
        (-screened_rates(k_p_Na23_to_He4_Ne20)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(He4)*Y(Ne20)*state.rho);

    ydot_nuc(Mg24) =
        (screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*Y(Ne20)*state.rho + -screened_rates(k_Mg24_to_He4_Ne20_derived)*Y(Mg24)) +
        (screened_rates(k_p_Na23_to_Mg24)*Y(Na23)*Y(H1)*state.rho + -screened_rates(k_Mg24_to_p_Na23_derived)*Y(Mg24)) +
        (-screened_rates(k_He4_Mg24_to_Si28)*Y(He4)*Y(Mg24)*state.rho + screened_rates(k_Si28_to_He4_Mg24_derived)*Y(Si28)) +
        screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*Y(O16)*state.rho +
        (screened_rates(k_p_Al27_to_He4_Mg24)*Y(Al27)*Y(H1)*state.rho + -screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(He4)*Y(Mg24)*state.rho) +
        0.5*screened_rates(k_C12_C12_to_Mg24_modified)*amrex::Math::powi<2>(Y(C12))*state.rho;

    ydot_nuc(Al27) =
        (-screened_rates(k_p_Al27_to_Si28)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_Si28_to_p_Al27_derived)*Y(Si28)) +
        screened_rates(k_C12_O16_to_p_Al27)*Y(C12)*Y(O16)*state.rho +
        (-screened_rates(k_p_Al27_to_He4_Mg24)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(He4)*Y(Mg24)*state.rho);

    ydot_nuc(Si28) =
        (screened_rates(k_He4_Mg24_to_Si28)*Y(He4)*Y(Mg24)*state.rho + -screened_rates(k_Si28_to_He4_Mg24_derived)*Y(Si28)) +
        (screened_rates(k_p_Al27_to_Si28)*Y(Al27)*Y(H1)*state.rho + -screened_rates(k_Si28_to_p_Al27_derived)*Y(Si28)) +
        (-screened_rates(k_He4_Si28_to_S32)*Y(He4)*Y(Si28)*state.rho + screened_rates(k_S32_to_He4_Si28_derived)*Y(S32)) +
        0.5*screened_rates(k_O16_O16_to_He4_Si28)*amrex::Math::powi<2>(Y(O16))*state.rho +
        (screened_rates(k_p_P31_to_He4_Si28)*Y(P31)*Y(H1)*state.rho + -screened_rates(k_He4_Si28_to_p_P31_derived)*Y(He4)*Y(Si28)*state.rho) +
        screened_rates(k_C12_O16_to_Si28_modified)*Y(C12)*Y(O16)*state.rho;

    ydot_nuc(P31) =
        (-screened_rates(k_p_P31_to_S32)*Y(P31)*Y(H1)*state.rho + screened_rates(k_S32_to_p_P31_derived)*Y(S32)) +
        0.5*screened_rates(k_O16_O16_to_p_P31)*amrex::Math::powi<2>(Y(O16))*state.rho +
        (-screened_rates(k_p_P31_to_He4_Si28)*Y(P31)*Y(H1)*state.rho + screened_rates(k_He4_Si28_to_p_P31_derived)*Y(He4)*Y(Si28)*state.rho);

    ydot_nuc(S32) =
        (screened_rates(k_He4_Si28_to_S32)*Y(He4)*Y(Si28)*state.rho + -screened_rates(k_S32_to_He4_Si28_derived)*Y(S32)) +
        (screened_rates(k_p_P31_to_S32)*Y(P31)*Y(H1)*state.rho + -screened_rates(k_S32_to_p_P31_derived)*Y(S32)) +
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
        (-screened_rates(k_He4_Cr48_to_Fe52)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_Fe52_to_He4_Cr48_derived)*Y(Fe52)) +
        (-screened_rates(k_He4_Cr48_to_p_nse_Mn51)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_p_nse_Mn51_to_He4_Cr48_derived)*Y(Mn51)*Y(P_nse)*state.rho) +
        (screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(He4)*Y(Ti44)*state.rho + -screened_rates(k_Cr48_to_Ti44_He4_approx)*Y(Cr48));

    ydot_nuc(Mn51) =
        (-screened_rates(k_p_nse_Mn51_to_Fe52)*Y(Mn51)*Y(P_nse)*state.rho + screened_rates(k_Fe52_to_p_nse_Mn51_derived)*Y(Fe52)) +
        (-screened_rates(k_He4_Mn51_to_Co55)*Y(He4)*Y(Mn51)*state.rho + screened_rates(k_Co55_to_He4_Mn51_derived)*Y(Co55)) +
        (screened_rates(k_He4_Cr48_to_p_nse_Mn51)*Y(Cr48)*Y(He4)*state.rho + -screened_rates(k_p_nse_Mn51_to_He4_Cr48_derived)*Y(Mn51)*Y(P_nse)*state.rho) +
        (-screened_rates(k_He4_Mn51_to_p_nse_Fe54)*Y(He4)*Y(Mn51)*state.rho + screened_rates(k_p_nse_Fe54_to_He4_Mn51_derived)*Y(Fe54)*Y(P_nse)*state.rho);

    ydot_nuc(Fe52) =
        (screened_rates(k_He4_Cr48_to_Fe52)*Y(Cr48)*Y(He4)*state.rho + -screened_rates(k_Fe52_to_He4_Cr48_derived)*Y(Fe52)) +
        (screened_rates(k_p_nse_Mn51_to_Fe52)*Y(Mn51)*Y(P_nse)*state.rho + -screened_rates(k_Fe52_to_p_nse_Mn51_derived)*Y(Fe52)) +
        (-screened_rates(k_He4_Fe52_to_Ni56)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_Ni56_to_He4_Fe52_derived)*Y(Ni56)) +
        (-screened_rates(k_He4_Fe52_to_p_nse_Co55)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_p_nse_Co55_to_He4_Fe52_derived)*Y(Co55)*Y(P_nse)*state.rho) +
        (-screened_rates(k_Fe52_n_n_to_Fe54_approx)*Y(Fe52)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + screened_rates(k_Fe54_to_Fe52_n_n_approx)*Y(Fe54));

    ydot_nuc(Fe54) =
        (-screened_rates(k_p_nse_Fe54_to_Co55)*Y(Fe54)*Y(P_nse)*state.rho + screened_rates(k_Co55_to_p_nse_Fe54_derived)*Y(Co55)) +
        (-screened_rates(k_He4_Fe54_to_Ni58)*Y(Fe54)*Y(He4)*state.rho + screened_rates(k_Ni58_to_He4_Fe54_derived)*Y(Ni58)) +
        (screened_rates(k_He4_Mn51_to_p_nse_Fe54)*Y(He4)*Y(Mn51)*state.rho + -screened_rates(k_p_nse_Fe54_to_He4_Mn51_derived)*Y(Fe54)*Y(P_nse)*state.rho) +
        (screened_rates(k_p_nse_Co57_to_He4_Fe54)*Y(Co57)*Y(P_nse)*state.rho + -screened_rates(k_He4_Fe54_to_p_nse_Co57_derived)*Y(Fe54)*Y(He4)*state.rho) +
        (screened_rates(k_Fe52_n_n_to_Fe54_approx)*Y(Fe52)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_Fe54_to_Fe52_n_n_approx)*Y(Fe54)) +
        (-screened_rates(k_Fe54_n_n_to_Fe56_approx)*Y(Fe54)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + screened_rates(k_Fe56_to_Fe54_n_n_approx)*Y(Fe56));

    ydot_nuc(Fe56) =
        (-screened_rates(k_p_nse_Fe56_to_Co57)*Y(Fe56)*Y(P_nse)*state.rho + screened_rates(k_Co57_to_p_nse_Fe56_derived)*Y(Co57)) +
        (screened_rates(k_n_Co56_to_p_nse_Fe56)*Y(Co56)*Y(N)*state.rho + -screened_rates(k_p_nse_Fe56_to_n_Co56_derived)*Y(Fe56)*Y(P_nse)*state.rho) +
        (screened_rates(k_Fe54_n_n_to_Fe56_approx)*Y(Fe54)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_Fe56_to_Fe54_n_n_approx)*Y(Fe56)) +
        (screened_rates(k_Co56_to_Fe56)*Y(Co56) + -screened_rates(k_Fe56_to_Co56)*Y(Fe56));

    ydot_nuc(Co55) =
        (screened_rates(k_He4_Mn51_to_Co55)*Y(He4)*Y(Mn51)*state.rho + -screened_rates(k_Co55_to_He4_Mn51_derived)*Y(Co55)) +
        (-screened_rates(k_p_nse_Co55_to_Ni56)*Y(Co55)*Y(P_nse)*state.rho + screened_rates(k_Ni56_to_p_nse_Co55_derived)*Y(Ni56)) +
        (screened_rates(k_He4_Fe52_to_p_nse_Co55)*Y(Fe52)*Y(He4)*state.rho + -screened_rates(k_p_nse_Co55_to_He4_Fe52_derived)*Y(Co55)*Y(P_nse)*state.rho) +
        (screened_rates(k_p_nse_Fe54_to_Co55)*Y(Fe54)*Y(P_nse)*state.rho + -screened_rates(k_Co55_to_p_nse_Fe54_derived)*Y(Co55)) +
        (-screened_rates(k_n_Co55_to_Co56)*Y(Co55)*Y(N)*state.rho + screened_rates(k_Co56_to_n_Co55_derived)*Y(Co56)) +
        (-screened_rates(k_He4_Co55_to_p_nse_Ni58)*Y(Co55)*Y(He4)*state.rho + screened_rates(k_p_nse_Ni58_to_He4_Co55_derived)*Y(Ni58)*Y(P_nse)*state.rho);

    ydot_nuc(Co56) =
        (screened_rates(k_n_Co55_to_Co56)*Y(Co55)*Y(N)*state.rho + -screened_rates(k_Co56_to_n_Co55_derived)*Y(Co56)) +
        (-screened_rates(k_n_Co56_to_Co57)*Y(Co56)*Y(N)*state.rho + screened_rates(k_Co57_to_n_Co56_derived)*Y(Co57)) +
        (-screened_rates(k_n_Co56_to_p_nse_Fe56)*Y(Co56)*Y(N)*state.rho + screened_rates(k_p_nse_Fe56_to_n_Co56_derived)*Y(Fe56)*Y(P_nse)*state.rho) +
        (screened_rates(k_n_Ni56_to_p_nse_Co56)*Y(Ni56)*Y(N)*state.rho + -screened_rates(k_p_nse_Co56_to_n_Ni56_derived)*Y(Co56)*Y(P_nse)*state.rho) +
        (-screened_rates(k_Co56_to_Fe56)*Y(Co56) + screened_rates(k_Fe56_to_Co56)*Y(Fe56)) +
        (screened_rates(k_Ni56_to_Co56)*Y(Ni56) + -screened_rates(k_Co56_to_Ni56)*Y(Co56));

    ydot_nuc(Co57) =
        (screened_rates(k_p_nse_Fe56_to_Co57)*Y(Fe56)*Y(P_nse)*state.rho + -screened_rates(k_Co57_to_p_nse_Fe56_derived)*Y(Co57)) +
        (screened_rates(k_n_Co56_to_Co57)*Y(Co56)*Y(N)*state.rho + -screened_rates(k_Co57_to_n_Co56_derived)*Y(Co57)) +
        (-screened_rates(k_p_nse_Co57_to_Ni58)*Y(Co57)*Y(P_nse)*state.rho + screened_rates(k_Ni58_to_p_nse_Co57_derived)*Y(Ni58)) +
        (-screened_rates(k_p_nse_Co57_to_He4_Fe54)*Y(Co57)*Y(P_nse)*state.rho + screened_rates(k_He4_Fe54_to_p_nse_Co57_derived)*Y(Fe54)*Y(He4)*state.rho);

    ydot_nuc(Ni56) =
        (screened_rates(k_He4_Fe52_to_Ni56)*Y(Fe52)*Y(He4)*state.rho + -screened_rates(k_Ni56_to_He4_Fe52_derived)*Y(Ni56)) +
        (screened_rates(k_p_nse_Co55_to_Ni56)*Y(Co55)*Y(P_nse)*state.rho + -screened_rates(k_Ni56_to_p_nse_Co55_derived)*Y(Ni56)) +
        (-screened_rates(k_n_Ni56_to_p_nse_Co56)*Y(Ni56)*Y(N)*state.rho + screened_rates(k_p_nse_Co56_to_n_Ni56_derived)*Y(Co56)*Y(P_nse)*state.rho) +
        (-screened_rates(k_Ni56_n_n_to_Ni58_approx)*Y(Ni56)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + screened_rates(k_Ni58_to_Ni56_n_n_approx)*Y(Ni58)) +
        (-screened_rates(k_Ni56_to_Co56)*Y(Ni56) + screened_rates(k_Co56_to_Ni56)*Y(Co56));

    ydot_nuc(Ni58) =
        (screened_rates(k_He4_Fe54_to_Ni58)*Y(Fe54)*Y(He4)*state.rho + -screened_rates(k_Ni58_to_He4_Fe54_derived)*Y(Ni58)) +
        (screened_rates(k_p_nse_Co57_to_Ni58)*Y(Co57)*Y(P_nse)*state.rho + -screened_rates(k_Ni58_to_p_nse_Co57_derived)*Y(Ni58)) +
        (screened_rates(k_He4_Co55_to_p_nse_Ni58)*Y(Co55)*Y(He4)*state.rho + -screened_rates(k_p_nse_Ni58_to_He4_Co55_derived)*Y(Ni58)*Y(P_nse)*state.rho) +
        (screened_rates(k_Ni56_n_n_to_Ni58_approx)*Y(Ni56)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_Ni58_to_Ni56_n_n_approx)*Y(Ni58));

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

    amrex::Real sneut{};
#ifdef NEUTRINOS
    constexpr int do_derivatives{0};
    amrex::Real dsneutdt{}, dsneutdd{}, dsnuda{}, dsnudz{};
    sneut5<do_derivatives>(state.T, state.rho, state.abar, state.zbar, sneut, dsneutdt, dsneutdd, dsnuda, dsnudz);
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

    scratch = -4.0*screened_rates(k_Fe52_n_n_to_Fe54_approx)*Y(Fe52)*Y(N)*amrex::Math::powi<2>(state.rho) - 4.0*screened_rates(k_Fe54_n_n_to_Fe56_approx)*Y(Fe54)*Y(N)*amrex::Math::powi<2>(state.rho) - 4.0*screened_rates(k_Ni56_n_n_to_Ni58_approx)*Y(Ni56)*Y(N)*amrex::Math::powi<2>(state.rho) - screened_rates(k_n_Co55_to_Co56)*Y(Co55)*state.rho - screened_rates(k_n_Co56_to_Co57)*Y(Co56)*state.rho - screened_rates(k_n_Co56_to_p_nse_Fe56)*Y(Co56)*state.rho - screened_rates(k_n_Ni56_to_p_nse_Co56)*Y(Ni56)*state.rho - screened_rates(k_n_to_p);
    jac.set(N, N, scratch);

    scratch = screened_rates(k_p_to_n);
    jac.set(N, H1, scratch);

    scratch = screened_rates(k_p_nse_Co56_to_n_Ni56_derived)*Y(Co56)*state.rho + screened_rates(k_p_nse_Fe56_to_n_Co56_derived)*Y(Fe56)*state.rho;
    jac.set(N, P_nse, scratch);

    scratch = -2.0*screened_rates(k_Fe52_n_n_to_Fe54_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(N, Fe52, scratch);

    scratch = -2.0*screened_rates(k_Fe54_n_n_to_Fe56_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + 2.0*screened_rates(k_Fe54_to_Fe52_n_n_approx);
    jac.set(N, Fe54, scratch);

    scratch = 2.0*screened_rates(k_Fe56_to_Fe54_n_n_approx) + screened_rates(k_p_nse_Fe56_to_n_Co56_derived)*Y(P_nse)*state.rho;
    jac.set(N, Fe56, scratch);

    scratch = -screened_rates(k_n_Co55_to_Co56)*Y(N)*state.rho;
    jac.set(N, Co55, scratch);

    scratch = screened_rates(k_Co56_to_n_Co55_derived) - screened_rates(k_n_Co56_to_Co57)*Y(N)*state.rho - screened_rates(k_n_Co56_to_p_nse_Fe56)*Y(N)*state.rho + screened_rates(k_p_nse_Co56_to_n_Ni56_derived)*Y(P_nse)*state.rho;
    jac.set(N, Co56, scratch);

    scratch = screened_rates(k_Co57_to_n_Co56_derived);
    jac.set(N, Co57, scratch);

    scratch = -2.0*screened_rates(k_Ni56_n_n_to_Ni58_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) - screened_rates(k_n_Ni56_to_p_nse_Co56)*Y(N)*state.rho;
    jac.set(N, Ni56, scratch);

    scratch = 2.0*screened_rates(k_Ni58_to_Ni56_n_n_approx);
    jac.set(N, Ni58, scratch);

    scratch = screened_rates(k_n_to_p);
    jac.set(H1, N, scratch);

    scratch = -screened_rates(k_p_Al27_to_He4_Mg24)*Y(Al27)*state.rho - screened_rates(k_p_Al27_to_Si28)*Y(Al27)*state.rho - screened_rates(k_p_C12_to_N13)*Y(C12)*state.rho - screened_rates(k_p_Na23_to_He4_Ne20)*Y(Na23)*state.rho - screened_rates(k_p_Na23_to_Mg24)*Y(Na23)*state.rho - screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(Ne21)*state.rho - screened_rates(k_p_Ne21_to_Na22)*Y(Ne21)*state.rho - screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*state.rho - screened_rates(k_p_P31_to_He4_Si28)*Y(P31)*state.rho - screened_rates(k_p_P31_to_S32)*Y(P31)*state.rho - screened_rates(k_p_to_n);
    jac.set(H1, H1, scratch);

    scratch = screened_rates(k_He4_F18_to_p_Ne21)*Y(F18)*state.rho + screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(Mg24)*state.rho + screened_rates(k_He4_N13_to_p_O16)*Y(N13)*state.rho + screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(Ne20)*state.rho + screened_rates(k_He4_Si28_to_p_P31_derived)*Y(Si28)*state.rho;
    jac.set(H1, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_p_Na23)*Y(C12)*state.rho + screened_rates(k_C12_O16_to_p_Al27)*Y(O16)*state.rho - screened_rates(k_p_C12_to_N13)*Y(H1)*state.rho;
    jac.set(H1, C12, scratch);

    scratch = screened_rates(k_He4_N13_to_p_O16)*Y(He4)*state.rho + screened_rates(k_N13_to_p_C12_derived);
    jac.set(H1, N13, scratch);

    scratch = screened_rates(k_C12_O16_to_p_Al27)*Y(C12)*state.rho + 1.0*screened_rates(k_O16_O16_to_p_P31)*Y(O16)*state.rho - screened_rates(k_p_O16_to_He4_N13_derived)*Y(H1)*state.rho;
    jac.set(H1, O16, scratch);

    scratch = screened_rates(k_He4_F18_to_p_Ne21)*Y(He4)*state.rho;
    jac.set(H1, F18, scratch);

    scratch = screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(He4)*state.rho;
    jac.set(H1, Ne20, scratch);

    scratch = -screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(H1)*state.rho - screened_rates(k_p_Ne21_to_Na22)*Y(H1)*state.rho;
    jac.set(H1, Ne21, scratch);

    scratch = screened_rates(k_Na22_to_p_Ne21_derived);
    jac.set(H1, Na22, scratch);

    scratch = -screened_rates(k_p_Na23_to_He4_Ne20)*Y(H1)*state.rho - screened_rates(k_p_Na23_to_Mg24)*Y(H1)*state.rho;
    jac.set(H1, Na23, scratch);

    scratch = screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(He4)*state.rho + screened_rates(k_Mg24_to_p_Na23_derived);
    jac.set(H1, Mg24, scratch);

    scratch = -screened_rates(k_p_Al27_to_He4_Mg24)*Y(H1)*state.rho - screened_rates(k_p_Al27_to_Si28)*Y(H1)*state.rho;
    jac.set(H1, Al27, scratch);

    scratch = screened_rates(k_He4_Si28_to_p_P31_derived)*Y(He4)*state.rho + screened_rates(k_Si28_to_p_Al27_derived);
    jac.set(H1, Si28, scratch);

    scratch = -screened_rates(k_p_P31_to_He4_Si28)*Y(H1)*state.rho - screened_rates(k_p_P31_to_S32)*Y(H1)*state.rho;
    jac.set(H1, P31, scratch);

    scratch = screened_rates(k_S32_to_p_P31_derived);
    jac.set(H1, S32, scratch);

    scratch = screened_rates(k_n_Co56_to_p_nse_Fe56)*Y(Co56)*state.rho + screened_rates(k_n_Ni56_to_p_nse_Co56)*Y(Ni56)*state.rho;
    jac.set(P_nse, N, scratch);

    scratch = -screened_rates(k_p_nse_Co55_to_He4_Fe52_derived)*Y(Co55)*state.rho - screened_rates(k_p_nse_Co55_to_Ni56)*Y(Co55)*state.rho - screened_rates(k_p_nse_Co56_to_n_Ni56_derived)*Y(Co56)*state.rho - screened_rates(k_p_nse_Co57_to_He4_Fe54)*Y(Co57)*state.rho - screened_rates(k_p_nse_Co57_to_Ni58)*Y(Co57)*state.rho - screened_rates(k_p_nse_Fe54_to_Co55)*Y(Fe54)*state.rho - screened_rates(k_p_nse_Fe54_to_He4_Mn51_derived)*Y(Fe54)*state.rho - screened_rates(k_p_nse_Fe56_to_Co57)*Y(Fe56)*state.rho - screened_rates(k_p_nse_Fe56_to_n_Co56_derived)*Y(Fe56)*state.rho - screened_rates(k_p_nse_Mn51_to_Fe52)*Y(Mn51)*state.rho - screened_rates(k_p_nse_Mn51_to_He4_Cr48_derived)*Y(Mn51)*state.rho - screened_rates(k_p_nse_Ni58_to_He4_Co55_derived)*Y(Ni58)*state.rho;
    jac.set(P_nse, P_nse, scratch);

    scratch = screened_rates(k_He4_Co55_to_p_nse_Ni58)*Y(Co55)*state.rho + screened_rates(k_He4_Cr48_to_p_nse_Mn51)*Y(Cr48)*state.rho + screened_rates(k_He4_Fe52_to_p_nse_Co55)*Y(Fe52)*state.rho + screened_rates(k_He4_Fe54_to_p_nse_Co57_derived)*Y(Fe54)*state.rho + screened_rates(k_He4_Mn51_to_p_nse_Fe54)*Y(Mn51)*state.rho;
    jac.set(P_nse, He4, scratch);

    scratch = screened_rates(k_He4_Cr48_to_p_nse_Mn51)*Y(He4)*state.rho;
    jac.set(P_nse, Cr48, scratch);

    scratch = screened_rates(k_He4_Mn51_to_p_nse_Fe54)*Y(He4)*state.rho - screened_rates(k_p_nse_Mn51_to_Fe52)*Y(P_nse)*state.rho - screened_rates(k_p_nse_Mn51_to_He4_Cr48_derived)*Y(P_nse)*state.rho;
    jac.set(P_nse, Mn51, scratch);

    scratch = screened_rates(k_Fe52_to_p_nse_Mn51_derived) + screened_rates(k_He4_Fe52_to_p_nse_Co55)*Y(He4)*state.rho;
    jac.set(P_nse, Fe52, scratch);

    scratch = screened_rates(k_He4_Fe54_to_p_nse_Co57_derived)*Y(He4)*state.rho - screened_rates(k_p_nse_Fe54_to_Co55)*Y(P_nse)*state.rho - screened_rates(k_p_nse_Fe54_to_He4_Mn51_derived)*Y(P_nse)*state.rho;
    jac.set(P_nse, Fe54, scratch);

    scratch = -screened_rates(k_p_nse_Fe56_to_Co57)*Y(P_nse)*state.rho - screened_rates(k_p_nse_Fe56_to_n_Co56_derived)*Y(P_nse)*state.rho;
    jac.set(P_nse, Fe56, scratch);

    scratch = screened_rates(k_Co55_to_p_nse_Fe54_derived) + screened_rates(k_He4_Co55_to_p_nse_Ni58)*Y(He4)*state.rho - screened_rates(k_p_nse_Co55_to_He4_Fe52_derived)*Y(P_nse)*state.rho - screened_rates(k_p_nse_Co55_to_Ni56)*Y(P_nse)*state.rho;
    jac.set(P_nse, Co55, scratch);

    scratch = screened_rates(k_n_Co56_to_p_nse_Fe56)*Y(N)*state.rho - screened_rates(k_p_nse_Co56_to_n_Ni56_derived)*Y(P_nse)*state.rho;
    jac.set(P_nse, Co56, scratch);

    scratch = screened_rates(k_Co57_to_p_nse_Fe56_derived) - screened_rates(k_p_nse_Co57_to_He4_Fe54)*Y(P_nse)*state.rho - screened_rates(k_p_nse_Co57_to_Ni58)*Y(P_nse)*state.rho;
    jac.set(P_nse, Co57, scratch);

    scratch = screened_rates(k_Ni56_to_p_nse_Co55_derived) + screened_rates(k_n_Ni56_to_p_nse_Co56)*Y(N)*state.rho;
    jac.set(P_nse, Ni56, scratch);

    scratch = screened_rates(k_Ni58_to_p_nse_Co57_derived) - screened_rates(k_p_nse_Ni58_to_He4_Co55_derived)*Y(P_nse)*state.rho;
    jac.set(P_nse, Ni58, scratch);

    scratch = screened_rates(k_p_Al27_to_He4_Mg24)*Y(Al27)*state.rho + screened_rates(k_p_Na23_to_He4_Ne20)*Y(Na23)*state.rho + screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(Ne21)*state.rho + screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*state.rho + screened_rates(k_p_P31_to_He4_Si28)*Y(P31)*state.rho;
    jac.set(He4, H1, scratch);

    scratch = screened_rates(k_p_nse_Co55_to_He4_Fe52_derived)*Y(Co55)*state.rho + screened_rates(k_p_nse_Co57_to_He4_Fe54)*Y(Co57)*state.rho + screened_rates(k_p_nse_Fe54_to_He4_Mn51_derived)*Y(Fe54)*state.rho + screened_rates(k_p_nse_Mn51_to_He4_Cr48_derived)*Y(Mn51)*state.rho + screened_rates(k_p_nse_Ni58_to_He4_Co55_derived)*Y(Ni58)*state.rho;
    jac.set(He4, P_nse, scratch);

    scratch = -screened_rates(k_Ar36_He4_to_Ca40_approx)*Y(Ar36)*state.rho - screened_rates(k_Ca40_He4_to_Ti44_approx)*Y(Ca40)*state.rho - screened_rates(k_He4_C12_to_O16)*Y(C12)*state.rho - screened_rates(k_He4_Co55_to_p_nse_Ni58)*Y(Co55)*state.rho - screened_rates(k_He4_Cr48_to_Fe52)*Y(Cr48)*state.rho - screened_rates(k_He4_Cr48_to_p_nse_Mn51)*Y(Cr48)*state.rho - screened_rates(k_He4_F18_to_Na22)*Y(F18)*state.rho - screened_rates(k_He4_F18_to_p_Ne21)*Y(F18)*state.rho - screened_rates(k_He4_Fe52_to_Ni56)*Y(Fe52)*state.rho - screened_rates(k_He4_Fe52_to_p_nse_Co55)*Y(Fe52)*state.rho - screened_rates(k_He4_Fe54_to_Ni58)*Y(Fe54)*state.rho - screened_rates(k_He4_Fe54_to_p_nse_Co57_derived)*Y(Fe54)*state.rho - 1.5*screened_rates(k_He4_He4_He4_to_C12)*amrex::Math::powi<2>(Y(He4))*amrex::Math::powi<2>(state.rho) - screened_rates(k_He4_Mg24_to_Si28)*Y(Mg24)*state.rho - screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(Mg24)*state.rho - screened_rates(k_He4_Mn51_to_Co55)*Y(Mn51)*state.rho - screened_rates(k_He4_Mn51_to_p_nse_Fe54)*Y(Mn51)*state.rho - screened_rates(k_He4_N13_to_p_O16)*Y(N13)*state.rho - screened_rates(k_He4_N14_to_F18)*Y(N14)*state.rho - screened_rates(k_He4_Ne20_to_Mg24)*Y(Ne20)*state.rho - screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(Ne20)*state.rho - screened_rates(k_He4_O16_to_Ne20)*Y(O16)*state.rho - screened_rates(k_He4_Si28_to_S32)*Y(Si28)*state.rho - screened_rates(k_He4_Si28_to_p_P31_derived)*Y(Si28)*state.rho - screened_rates(k_S32_He4_to_Ar36_approx)*Y(S32)*state.rho - screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(Ti44)*state.rho;
    jac.set(He4, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_He4_Ne20)*Y(C12)*state.rho + screened_rates(k_C12_O16_to_He4_Mg24)*Y(O16)*state.rho + 3.0*screened_rates(k_C12_to_He4_He4_He4_derived) - screened_rates(k_He4_C12_to_O16)*Y(He4)*state.rho;
    jac.set(He4, C12, scratch);

    scratch = -screened_rates(k_He4_N13_to_p_O16)*Y(He4)*state.rho;
    jac.set(He4, N13, scratch);

    scratch = -screened_rates(k_He4_N14_to_F18)*Y(He4)*state.rho;
    jac.set(He4, N14, scratch);

    scratch = screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*state.rho - screened_rates(k_He4_O16_to_Ne20)*Y(He4)*state.rho + 1.0*screened_rates(k_O16_O16_to_He4_Si28)*Y(O16)*state.rho + screened_rates(k_O16_to_He4_C12_derived) + screened_rates(k_p_O16_to_He4_N13_derived)*Y(H1)*state.rho;
    jac.set(He4, O16, scratch);

    scratch = screened_rates(k_F18_to_He4_N14_derived) - screened_rates(k_He4_F18_to_Na22)*Y(He4)*state.rho - screened_rates(k_He4_F18_to_p_Ne21)*Y(He4)*state.rho;
    jac.set(He4, F18, scratch);

    scratch = -screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*state.rho - screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(He4)*state.rho + screened_rates(k_Ne20_to_He4_O16_derived);
    jac.set(He4, Ne20, scratch);

    scratch = screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(H1)*state.rho;
    jac.set(He4, Ne21, scratch);

    scratch = screened_rates(k_Na22_to_He4_F18_derived);
    jac.set(He4, Na22, scratch);

    scratch = screened_rates(k_p_Na23_to_He4_Ne20)*Y(H1)*state.rho;
    jac.set(He4, Na23, scratch);

    scratch = -screened_rates(k_He4_Mg24_to_Si28)*Y(He4)*state.rho - screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(He4)*state.rho + screened_rates(k_Mg24_to_He4_Ne20_derived);
    jac.set(He4, Mg24, scratch);

    scratch = screened_rates(k_p_Al27_to_He4_Mg24)*Y(H1)*state.rho;
    jac.set(He4, Al27, scratch);

    scratch = -screened_rates(k_He4_Si28_to_S32)*Y(He4)*state.rho - screened_rates(k_He4_Si28_to_p_P31_derived)*Y(He4)*state.rho + screened_rates(k_Si28_to_He4_Mg24_derived);
    jac.set(He4, Si28, scratch);

    scratch = screened_rates(k_p_P31_to_He4_Si28)*Y(H1)*state.rho;
    jac.set(He4, P31, scratch);

    scratch = -screened_rates(k_S32_He4_to_Ar36_approx)*Y(He4)*state.rho + screened_rates(k_S32_to_He4_Si28_derived);
    jac.set(He4, S32, scratch);

    scratch = -screened_rates(k_Ar36_He4_to_Ca40_approx)*Y(He4)*state.rho + screened_rates(k_Ar36_to_S32_He4_approx);
    jac.set(He4, Ar36, scratch);

    scratch = -screened_rates(k_Ca40_He4_to_Ti44_approx)*Y(He4)*state.rho + screened_rates(k_Ca40_to_Ar36_He4_approx);
    jac.set(He4, Ca40, scratch);

    scratch = -screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(He4)*state.rho + screened_rates(k_Ti44_to_Ca40_He4_approx);
    jac.set(He4, Ti44, scratch);

    scratch = screened_rates(k_Cr48_to_Ti44_He4_approx) - screened_rates(k_He4_Cr48_to_Fe52)*Y(He4)*state.rho - screened_rates(k_He4_Cr48_to_p_nse_Mn51)*Y(He4)*state.rho;
    jac.set(He4, Cr48, scratch);

    scratch = -screened_rates(k_He4_Mn51_to_Co55)*Y(He4)*state.rho - screened_rates(k_He4_Mn51_to_p_nse_Fe54)*Y(He4)*state.rho + screened_rates(k_p_nse_Mn51_to_He4_Cr48_derived)*Y(P_nse)*state.rho;
    jac.set(He4, Mn51, scratch);

    scratch = screened_rates(k_Fe52_to_He4_Cr48_derived) - screened_rates(k_He4_Fe52_to_Ni56)*Y(He4)*state.rho - screened_rates(k_He4_Fe52_to_p_nse_Co55)*Y(He4)*state.rho;
    jac.set(He4, Fe52, scratch);

    scratch = -screened_rates(k_He4_Fe54_to_Ni58)*Y(He4)*state.rho - screened_rates(k_He4_Fe54_to_p_nse_Co57_derived)*Y(He4)*state.rho + screened_rates(k_p_nse_Fe54_to_He4_Mn51_derived)*Y(P_nse)*state.rho;
    jac.set(He4, Fe54, scratch);

    scratch = screened_rates(k_Co55_to_He4_Mn51_derived) - screened_rates(k_He4_Co55_to_p_nse_Ni58)*Y(He4)*state.rho + screened_rates(k_p_nse_Co55_to_He4_Fe52_derived)*Y(P_nse)*state.rho;
    jac.set(He4, Co55, scratch);

    scratch = screened_rates(k_p_nse_Co57_to_He4_Fe54)*Y(P_nse)*state.rho;
    jac.set(He4, Co57, scratch);

    scratch = screened_rates(k_Ni56_to_He4_Fe52_derived);
    jac.set(He4, Ni56, scratch);

    scratch = screened_rates(k_Ni58_to_He4_Fe54_derived) + screened_rates(k_p_nse_Ni58_to_He4_Co55_derived)*Y(P_nse)*state.rho;
    jac.set(He4, Ni58, scratch);

    scratch = -screened_rates(k_p_C12_to_N13)*Y(C12)*state.rho;
    jac.set(C12, H1, scratch);

    scratch = -screened_rates(k_He4_C12_to_O16)*Y(C12)*state.rho + 0.5*screened_rates(k_He4_He4_He4_to_C12)*amrex::Math::powi<2>(Y(He4))*amrex::Math::powi<2>(state.rho);
    jac.set(C12, He4, scratch);

    scratch = -2.0*screened_rates(k_C12_C12_to_He4_Ne20)*Y(C12)*state.rho - 2.0*screened_rates(k_C12_C12_to_Mg24_modified)*Y(C12)*state.rho - 2.0*screened_rates(k_C12_C12_to_p_Na23)*Y(C12)*state.rho - screened_rates(k_C12_O16_to_He4_Mg24)*Y(O16)*state.rho - screened_rates(k_C12_O16_to_Si28_modified)*Y(O16)*state.rho - screened_rates(k_C12_O16_to_p_Al27)*Y(O16)*state.rho - screened_rates(k_C12_to_He4_He4_He4_derived) - screened_rates(k_He4_C12_to_O16)*Y(He4)*state.rho - screened_rates(k_p_C12_to_N13)*Y(H1)*state.rho;
    jac.set(C12, C12, scratch);

    scratch = screened_rates(k_N13_to_p_C12_derived);
    jac.set(C12, N13, scratch);

    scratch = -screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*state.rho - screened_rates(k_C12_O16_to_Si28_modified)*Y(C12)*state.rho - screened_rates(k_C12_O16_to_p_Al27)*Y(C12)*state.rho + screened_rates(k_O16_to_He4_C12_derived);
    jac.set(C12, O16, scratch);

    scratch = screened_rates(k_p_C12_to_N13)*Y(C12)*state.rho + screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*state.rho;
    jac.set(N13, H1, scratch);

    scratch = -screened_rates(k_He4_N13_to_p_O16)*Y(N13)*state.rho;
    jac.set(N13, He4, scratch);

    scratch = screened_rates(k_p_C12_to_N13)*Y(H1)*state.rho;
    jac.set(N13, C12, scratch);

    scratch = -screened_rates(k_He4_N13_to_p_O16)*Y(He4)*state.rho - screened_rates(k_N13_to_p_C12_derived);
    jac.set(N13, N13, scratch);

    scratch = screened_rates(k_p_O16_to_He4_N13_derived)*Y(H1)*state.rho;
    jac.set(N13, O16, scratch);

    scratch = -screened_rates(k_He4_N14_to_F18)*Y(N14)*state.rho;
    jac.set(N14, He4, scratch);

    scratch = -screened_rates(k_He4_N14_to_F18)*Y(He4)*state.rho;
    jac.set(N14, N14, scratch);

    scratch = screened_rates(k_F18_to_He4_N14_derived);
    jac.set(N14, F18, scratch);

    scratch = -screened_rates(k_p_O16_to_He4_N13_derived)*Y(O16)*state.rho;
    jac.set(O16, H1, scratch);

    scratch = screened_rates(k_He4_C12_to_O16)*Y(C12)*state.rho + screened_rates(k_He4_N13_to_p_O16)*Y(N13)*state.rho - screened_rates(k_He4_O16_to_Ne20)*Y(O16)*state.rho;
    jac.set(O16, He4, scratch);

    scratch = -screened_rates(k_C12_O16_to_He4_Mg24)*Y(O16)*state.rho - screened_rates(k_C12_O16_to_Si28_modified)*Y(O16)*state.rho - screened_rates(k_C12_O16_to_p_Al27)*Y(O16)*state.rho + screened_rates(k_He4_C12_to_O16)*Y(He4)*state.rho;
    jac.set(O16, C12, scratch);

    scratch = screened_rates(k_He4_N13_to_p_O16)*Y(He4)*state.rho;
    jac.set(O16, N13, scratch);

    scratch = -screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*state.rho - screened_rates(k_C12_O16_to_Si28_modified)*Y(C12)*state.rho - screened_rates(k_C12_O16_to_p_Al27)*Y(C12)*state.rho - screened_rates(k_He4_O16_to_Ne20)*Y(He4)*state.rho - 2.0*screened_rates(k_O16_O16_to_He4_Si28)*Y(O16)*state.rho - 2.0*screened_rates(k_O16_O16_to_S32_modified)*Y(O16)*state.rho - 2.0*screened_rates(k_O16_O16_to_p_P31)*Y(O16)*state.rho - screened_rates(k_O16_to_He4_C12_derived) - screened_rates(k_p_O16_to_He4_N13_derived)*Y(H1)*state.rho;
    jac.set(O16, O16, scratch);

    scratch = screened_rates(k_Ne20_to_He4_O16_derived);
    jac.set(O16, Ne20, scratch);

    scratch = screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(Ne21)*state.rho;
    jac.set(F18, H1, scratch);

    scratch = -screened_rates(k_He4_F18_to_Na22)*Y(F18)*state.rho - screened_rates(k_He4_F18_to_p_Ne21)*Y(F18)*state.rho + screened_rates(k_He4_N14_to_F18)*Y(N14)*state.rho;
    jac.set(F18, He4, scratch);

    scratch = screened_rates(k_He4_N14_to_F18)*Y(He4)*state.rho;
    jac.set(F18, N14, scratch);

    scratch = -screened_rates(k_F18_to_He4_N14_derived) - screened_rates(k_He4_F18_to_Na22)*Y(He4)*state.rho - screened_rates(k_He4_F18_to_p_Ne21)*Y(He4)*state.rho;
    jac.set(F18, F18, scratch);

    scratch = screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(H1)*state.rho;
    jac.set(F18, Ne21, scratch);

    scratch = screened_rates(k_Na22_to_He4_F18_derived);
    jac.set(F18, Na22, scratch);

    scratch = screened_rates(k_p_Na23_to_He4_Ne20)*Y(Na23)*state.rho;
    jac.set(Ne20, H1, scratch);

    scratch = -screened_rates(k_He4_Ne20_to_Mg24)*Y(Ne20)*state.rho - screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(Ne20)*state.rho + screened_rates(k_He4_O16_to_Ne20)*Y(O16)*state.rho;
    jac.set(Ne20, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_He4_Ne20)*Y(C12)*state.rho;
    jac.set(Ne20, C12, scratch);

    scratch = screened_rates(k_He4_O16_to_Ne20)*Y(He4)*state.rho;
    jac.set(Ne20, O16, scratch);

    scratch = -screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*state.rho - screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(He4)*state.rho - screened_rates(k_Ne20_to_He4_O16_derived);
    jac.set(Ne20, Ne20, scratch);

    scratch = screened_rates(k_p_Na23_to_He4_Ne20)*Y(H1)*state.rho;
    jac.set(Ne20, Na23, scratch);

    scratch = screened_rates(k_Mg24_to_He4_Ne20_derived);
    jac.set(Ne20, Mg24, scratch);

    scratch = -screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(Ne21)*state.rho - screened_rates(k_p_Ne21_to_Na22)*Y(Ne21)*state.rho;
    jac.set(Ne21, H1, scratch);

    scratch = screened_rates(k_He4_F18_to_p_Ne21)*Y(F18)*state.rho;
    jac.set(Ne21, He4, scratch);

    scratch = screened_rates(k_He4_F18_to_p_Ne21)*Y(He4)*state.rho;
    jac.set(Ne21, F18, scratch);

    scratch = -screened_rates(k_p_Ne21_to_He4_F18_derived)*Y(H1)*state.rho - screened_rates(k_p_Ne21_to_Na22)*Y(H1)*state.rho;
    jac.set(Ne21, Ne21, scratch);

    scratch = screened_rates(k_Na22_to_p_Ne21_derived);
    jac.set(Ne21, Na22, scratch);

    scratch = screened_rates(k_p_Ne21_to_Na22)*Y(Ne21)*state.rho;
    jac.set(Na22, H1, scratch);

    scratch = screened_rates(k_He4_F18_to_Na22)*Y(F18)*state.rho;
    jac.set(Na22, He4, scratch);

    scratch = screened_rates(k_He4_F18_to_Na22)*Y(He4)*state.rho;
    jac.set(Na22, F18, scratch);

    scratch = screened_rates(k_p_Ne21_to_Na22)*Y(H1)*state.rho;
    jac.set(Na22, Ne21, scratch);

    scratch = -screened_rates(k_Na22_to_He4_F18_derived) - screened_rates(k_Na22_to_p_Ne21_derived);
    jac.set(Na22, Na22, scratch);

    scratch = -screened_rates(k_p_Na23_to_He4_Ne20)*Y(Na23)*state.rho - screened_rates(k_p_Na23_to_Mg24)*Y(Na23)*state.rho;
    jac.set(Na23, H1, scratch);

    scratch = screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(Ne20)*state.rho;
    jac.set(Na23, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_p_Na23)*Y(C12)*state.rho;
    jac.set(Na23, C12, scratch);

    scratch = screened_rates(k_He4_Ne20_to_p_Na23_derived)*Y(He4)*state.rho;
    jac.set(Na23, Ne20, scratch);

    scratch = -screened_rates(k_p_Na23_to_He4_Ne20)*Y(H1)*state.rho - screened_rates(k_p_Na23_to_Mg24)*Y(H1)*state.rho;
    jac.set(Na23, Na23, scratch);

    scratch = screened_rates(k_Mg24_to_p_Na23_derived);
    jac.set(Na23, Mg24, scratch);

    scratch = screened_rates(k_p_Al27_to_He4_Mg24)*Y(Al27)*state.rho + screened_rates(k_p_Na23_to_Mg24)*Y(Na23)*state.rho;
    jac.set(Mg24, H1, scratch);

    scratch = -screened_rates(k_He4_Mg24_to_Si28)*Y(Mg24)*state.rho - screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(Mg24)*state.rho + screened_rates(k_He4_Ne20_to_Mg24)*Y(Ne20)*state.rho;
    jac.set(Mg24, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_Mg24_modified)*Y(C12)*state.rho + screened_rates(k_C12_O16_to_He4_Mg24)*Y(O16)*state.rho;
    jac.set(Mg24, C12, scratch);

    scratch = screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*state.rho;
    jac.set(Mg24, O16, scratch);

    scratch = screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*state.rho;
    jac.set(Mg24, Ne20, scratch);

    scratch = screened_rates(k_p_Na23_to_Mg24)*Y(H1)*state.rho;
    jac.set(Mg24, Na23, scratch);

    scratch = -screened_rates(k_He4_Mg24_to_Si28)*Y(He4)*state.rho - screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(He4)*state.rho - screened_rates(k_Mg24_to_He4_Ne20_derived) - screened_rates(k_Mg24_to_p_Na23_derived);
    jac.set(Mg24, Mg24, scratch);

    scratch = screened_rates(k_p_Al27_to_He4_Mg24)*Y(H1)*state.rho;
    jac.set(Mg24, Al27, scratch);

    scratch = screened_rates(k_Si28_to_He4_Mg24_derived);
    jac.set(Mg24, Si28, scratch);

    scratch = -screened_rates(k_p_Al27_to_He4_Mg24)*Y(Al27)*state.rho - screened_rates(k_p_Al27_to_Si28)*Y(Al27)*state.rho;
    jac.set(Al27, H1, scratch);

    scratch = screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(Mg24)*state.rho;
    jac.set(Al27, He4, scratch);

    scratch = screened_rates(k_C12_O16_to_p_Al27)*Y(O16)*state.rho;
    jac.set(Al27, C12, scratch);

    scratch = screened_rates(k_C12_O16_to_p_Al27)*Y(C12)*state.rho;
    jac.set(Al27, O16, scratch);

    scratch = screened_rates(k_He4_Mg24_to_p_Al27_derived)*Y(He4)*state.rho;
    jac.set(Al27, Mg24, scratch);

    scratch = -screened_rates(k_p_Al27_to_He4_Mg24)*Y(H1)*state.rho - screened_rates(k_p_Al27_to_Si28)*Y(H1)*state.rho;
    jac.set(Al27, Al27, scratch);

    scratch = screened_rates(k_Si28_to_p_Al27_derived);
    jac.set(Al27, Si28, scratch);

    scratch = screened_rates(k_p_Al27_to_Si28)*Y(Al27)*state.rho + screened_rates(k_p_P31_to_He4_Si28)*Y(P31)*state.rho;
    jac.set(Si28, H1, scratch);

    scratch = screened_rates(k_He4_Mg24_to_Si28)*Y(Mg24)*state.rho - screened_rates(k_He4_Si28_to_S32)*Y(Si28)*state.rho - screened_rates(k_He4_Si28_to_p_P31_derived)*Y(Si28)*state.rho;
    jac.set(Si28, He4, scratch);

    scratch = screened_rates(k_C12_O16_to_Si28_modified)*Y(O16)*state.rho;
    jac.set(Si28, C12, scratch);

    scratch = screened_rates(k_C12_O16_to_Si28_modified)*Y(C12)*state.rho + 1.0*screened_rates(k_O16_O16_to_He4_Si28)*Y(O16)*state.rho;
    jac.set(Si28, O16, scratch);

    scratch = screened_rates(k_He4_Mg24_to_Si28)*Y(He4)*state.rho;
    jac.set(Si28, Mg24, scratch);

    scratch = screened_rates(k_p_Al27_to_Si28)*Y(H1)*state.rho;
    jac.set(Si28, Al27, scratch);

    scratch = -screened_rates(k_He4_Si28_to_S32)*Y(He4)*state.rho - screened_rates(k_He4_Si28_to_p_P31_derived)*Y(He4)*state.rho - screened_rates(k_Si28_to_He4_Mg24_derived) - screened_rates(k_Si28_to_p_Al27_derived);
    jac.set(Si28, Si28, scratch);

    scratch = screened_rates(k_p_P31_to_He4_Si28)*Y(H1)*state.rho;
    jac.set(Si28, P31, scratch);

    scratch = screened_rates(k_S32_to_He4_Si28_derived);
    jac.set(Si28, S32, scratch);

    scratch = -screened_rates(k_p_P31_to_He4_Si28)*Y(P31)*state.rho - screened_rates(k_p_P31_to_S32)*Y(P31)*state.rho;
    jac.set(P31, H1, scratch);

    scratch = screened_rates(k_He4_Si28_to_p_P31_derived)*Y(Si28)*state.rho;
    jac.set(P31, He4, scratch);

    scratch = 1.0*screened_rates(k_O16_O16_to_p_P31)*Y(O16)*state.rho;
    jac.set(P31, O16, scratch);

    scratch = screened_rates(k_He4_Si28_to_p_P31_derived)*Y(He4)*state.rho;
    jac.set(P31, Si28, scratch);

    scratch = -screened_rates(k_p_P31_to_He4_Si28)*Y(H1)*state.rho - screened_rates(k_p_P31_to_S32)*Y(H1)*state.rho;
    jac.set(P31, P31, scratch);

    scratch = screened_rates(k_S32_to_p_P31_derived);
    jac.set(P31, S32, scratch);

    scratch = screened_rates(k_p_P31_to_S32)*Y(P31)*state.rho;
    jac.set(S32, H1, scratch);

    scratch = screened_rates(k_He4_Si28_to_S32)*Y(Si28)*state.rho - screened_rates(k_S32_He4_to_Ar36_approx)*Y(S32)*state.rho;
    jac.set(S32, He4, scratch);

    scratch = 1.0*screened_rates(k_O16_O16_to_S32_modified)*Y(O16)*state.rho;
    jac.set(S32, O16, scratch);

    scratch = screened_rates(k_He4_Si28_to_S32)*Y(He4)*state.rho;
    jac.set(S32, Si28, scratch);

    scratch = screened_rates(k_p_P31_to_S32)*Y(H1)*state.rho;
    jac.set(S32, P31, scratch);

    scratch = -screened_rates(k_S32_He4_to_Ar36_approx)*Y(He4)*state.rho - screened_rates(k_S32_to_He4_Si28_derived) - screened_rates(k_S32_to_p_P31_derived);
    jac.set(S32, S32, scratch);

    scratch = screened_rates(k_Ar36_to_S32_He4_approx);
    jac.set(S32, Ar36, scratch);

    scratch = -screened_rates(k_Ar36_He4_to_Ca40_approx)*Y(Ar36)*state.rho + screened_rates(k_S32_He4_to_Ar36_approx)*Y(S32)*state.rho;
    jac.set(Ar36, He4, scratch);

    scratch = screened_rates(k_S32_He4_to_Ar36_approx)*Y(He4)*state.rho;
    jac.set(Ar36, S32, scratch);

    scratch = -screened_rates(k_Ar36_He4_to_Ca40_approx)*Y(He4)*state.rho - screened_rates(k_Ar36_to_S32_He4_approx);
    jac.set(Ar36, Ar36, scratch);

    scratch = screened_rates(k_Ca40_to_Ar36_He4_approx);
    jac.set(Ar36, Ca40, scratch);

    scratch = screened_rates(k_Ar36_He4_to_Ca40_approx)*Y(Ar36)*state.rho - screened_rates(k_Ca40_He4_to_Ti44_approx)*Y(Ca40)*state.rho;
    jac.set(Ca40, He4, scratch);

    scratch = screened_rates(k_Ar36_He4_to_Ca40_approx)*Y(He4)*state.rho;
    jac.set(Ca40, Ar36, scratch);

    scratch = -screened_rates(k_Ca40_He4_to_Ti44_approx)*Y(He4)*state.rho - screened_rates(k_Ca40_to_Ar36_He4_approx);
    jac.set(Ca40, Ca40, scratch);

    scratch = screened_rates(k_Ti44_to_Ca40_He4_approx);
    jac.set(Ca40, Ti44, scratch);

    scratch = screened_rates(k_Ca40_He4_to_Ti44_approx)*Y(Ca40)*state.rho - screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(Ti44)*state.rho;
    jac.set(Ti44, He4, scratch);

    scratch = screened_rates(k_Ca40_He4_to_Ti44_approx)*Y(He4)*state.rho;
    jac.set(Ti44, Ca40, scratch);

    scratch = -screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(He4)*state.rho - screened_rates(k_Ti44_to_Ca40_He4_approx);
    jac.set(Ti44, Ti44, scratch);

    scratch = screened_rates(k_Cr48_to_Ti44_He4_approx);
    jac.set(Ti44, Cr48, scratch);

    scratch = screened_rates(k_p_nse_Mn51_to_He4_Cr48_derived)*Y(Mn51)*state.rho;
    jac.set(Cr48, P_nse, scratch);

    scratch = -screened_rates(k_He4_Cr48_to_Fe52)*Y(Cr48)*state.rho - screened_rates(k_He4_Cr48_to_p_nse_Mn51)*Y(Cr48)*state.rho + screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(Ti44)*state.rho;
    jac.set(Cr48, He4, scratch);

    scratch = screened_rates(k_Ti44_He4_to_Cr48_approx)*Y(He4)*state.rho;
    jac.set(Cr48, Ti44, scratch);

    scratch = -screened_rates(k_Cr48_to_Ti44_He4_approx) - screened_rates(k_He4_Cr48_to_Fe52)*Y(He4)*state.rho - screened_rates(k_He4_Cr48_to_p_nse_Mn51)*Y(He4)*state.rho;
    jac.set(Cr48, Cr48, scratch);

    scratch = screened_rates(k_p_nse_Mn51_to_He4_Cr48_derived)*Y(P_nse)*state.rho;
    jac.set(Cr48, Mn51, scratch);

    scratch = screened_rates(k_Fe52_to_He4_Cr48_derived);
    jac.set(Cr48, Fe52, scratch);

    scratch = screened_rates(k_p_nse_Fe54_to_He4_Mn51_derived)*Y(Fe54)*state.rho - screened_rates(k_p_nse_Mn51_to_Fe52)*Y(Mn51)*state.rho - screened_rates(k_p_nse_Mn51_to_He4_Cr48_derived)*Y(Mn51)*state.rho;
    jac.set(Mn51, P_nse, scratch);

    scratch = screened_rates(k_He4_Cr48_to_p_nse_Mn51)*Y(Cr48)*state.rho - screened_rates(k_He4_Mn51_to_Co55)*Y(Mn51)*state.rho - screened_rates(k_He4_Mn51_to_p_nse_Fe54)*Y(Mn51)*state.rho;
    jac.set(Mn51, He4, scratch);

    scratch = screened_rates(k_He4_Cr48_to_p_nse_Mn51)*Y(He4)*state.rho;
    jac.set(Mn51, Cr48, scratch);

    scratch = -screened_rates(k_He4_Mn51_to_Co55)*Y(He4)*state.rho - screened_rates(k_He4_Mn51_to_p_nse_Fe54)*Y(He4)*state.rho - screened_rates(k_p_nse_Mn51_to_Fe52)*Y(P_nse)*state.rho - screened_rates(k_p_nse_Mn51_to_He4_Cr48_derived)*Y(P_nse)*state.rho;
    jac.set(Mn51, Mn51, scratch);

    scratch = screened_rates(k_Fe52_to_p_nse_Mn51_derived);
    jac.set(Mn51, Fe52, scratch);

    scratch = screened_rates(k_p_nse_Fe54_to_He4_Mn51_derived)*Y(P_nse)*state.rho;
    jac.set(Mn51, Fe54, scratch);

    scratch = screened_rates(k_Co55_to_He4_Mn51_derived);
    jac.set(Mn51, Co55, scratch);

    scratch = -2.0*screened_rates(k_Fe52_n_n_to_Fe54_approx)*Y(Fe52)*Y(N)*amrex::Math::powi<2>(state.rho);
    jac.set(Fe52, N, scratch);

    scratch = screened_rates(k_p_nse_Co55_to_He4_Fe52_derived)*Y(Co55)*state.rho + screened_rates(k_p_nse_Mn51_to_Fe52)*Y(Mn51)*state.rho;
    jac.set(Fe52, P_nse, scratch);

    scratch = screened_rates(k_He4_Cr48_to_Fe52)*Y(Cr48)*state.rho - screened_rates(k_He4_Fe52_to_Ni56)*Y(Fe52)*state.rho - screened_rates(k_He4_Fe52_to_p_nse_Co55)*Y(Fe52)*state.rho;
    jac.set(Fe52, He4, scratch);

    scratch = screened_rates(k_He4_Cr48_to_Fe52)*Y(He4)*state.rho;
    jac.set(Fe52, Cr48, scratch);

    scratch = screened_rates(k_p_nse_Mn51_to_Fe52)*Y(P_nse)*state.rho;
    jac.set(Fe52, Mn51, scratch);

    scratch = -screened_rates(k_Fe52_n_n_to_Fe54_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) - screened_rates(k_Fe52_to_He4_Cr48_derived) - screened_rates(k_Fe52_to_p_nse_Mn51_derived) - screened_rates(k_He4_Fe52_to_Ni56)*Y(He4)*state.rho - screened_rates(k_He4_Fe52_to_p_nse_Co55)*Y(He4)*state.rho;
    jac.set(Fe52, Fe52, scratch);

    scratch = screened_rates(k_Fe54_to_Fe52_n_n_approx);
    jac.set(Fe52, Fe54, scratch);

    scratch = screened_rates(k_p_nse_Co55_to_He4_Fe52_derived)*Y(P_nse)*state.rho;
    jac.set(Fe52, Co55, scratch);

    scratch = screened_rates(k_Ni56_to_He4_Fe52_derived);
    jac.set(Fe52, Ni56, scratch);

    scratch = 2.0*screened_rates(k_Fe52_n_n_to_Fe54_approx)*Y(Fe52)*Y(N)*amrex::Math::powi<2>(state.rho) - 2.0*screened_rates(k_Fe54_n_n_to_Fe56_approx)*Y(Fe54)*Y(N)*amrex::Math::powi<2>(state.rho);
    jac.set(Fe54, N, scratch);

    scratch = screened_rates(k_p_nse_Co57_to_He4_Fe54)*Y(Co57)*state.rho - screened_rates(k_p_nse_Fe54_to_Co55)*Y(Fe54)*state.rho - screened_rates(k_p_nse_Fe54_to_He4_Mn51_derived)*Y(Fe54)*state.rho;
    jac.set(Fe54, P_nse, scratch);

    scratch = -screened_rates(k_He4_Fe54_to_Ni58)*Y(Fe54)*state.rho - screened_rates(k_He4_Fe54_to_p_nse_Co57_derived)*Y(Fe54)*state.rho + screened_rates(k_He4_Mn51_to_p_nse_Fe54)*Y(Mn51)*state.rho;
    jac.set(Fe54, He4, scratch);

    scratch = screened_rates(k_He4_Mn51_to_p_nse_Fe54)*Y(He4)*state.rho;
    jac.set(Fe54, Mn51, scratch);

    scratch = screened_rates(k_Fe52_n_n_to_Fe54_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(Fe54, Fe52, scratch);

    scratch = -screened_rates(k_Fe54_n_n_to_Fe56_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) - screened_rates(k_Fe54_to_Fe52_n_n_approx) - screened_rates(k_He4_Fe54_to_Ni58)*Y(He4)*state.rho - screened_rates(k_He4_Fe54_to_p_nse_Co57_derived)*Y(He4)*state.rho - screened_rates(k_p_nse_Fe54_to_Co55)*Y(P_nse)*state.rho - screened_rates(k_p_nse_Fe54_to_He4_Mn51_derived)*Y(P_nse)*state.rho;
    jac.set(Fe54, Fe54, scratch);

    scratch = screened_rates(k_Fe56_to_Fe54_n_n_approx);
    jac.set(Fe54, Fe56, scratch);

    scratch = screened_rates(k_Co55_to_p_nse_Fe54_derived);
    jac.set(Fe54, Co55, scratch);

    scratch = screened_rates(k_p_nse_Co57_to_He4_Fe54)*Y(P_nse)*state.rho;
    jac.set(Fe54, Co57, scratch);

    scratch = screened_rates(k_Ni58_to_He4_Fe54_derived);
    jac.set(Fe54, Ni58, scratch);

    scratch = 2.0*screened_rates(k_Fe54_n_n_to_Fe56_approx)*Y(Fe54)*Y(N)*amrex::Math::powi<2>(state.rho) + screened_rates(k_n_Co56_to_p_nse_Fe56)*Y(Co56)*state.rho;
    jac.set(Fe56, N, scratch);

    scratch = -screened_rates(k_p_nse_Fe56_to_Co57)*Y(Fe56)*state.rho - screened_rates(k_p_nse_Fe56_to_n_Co56_derived)*Y(Fe56)*state.rho;
    jac.set(Fe56, P_nse, scratch);

    scratch = screened_rates(k_Fe54_n_n_to_Fe56_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(Fe56, Fe54, scratch);

    scratch = -screened_rates(k_Fe56_to_Co56) - screened_rates(k_Fe56_to_Fe54_n_n_approx) - screened_rates(k_p_nse_Fe56_to_Co57)*Y(P_nse)*state.rho - screened_rates(k_p_nse_Fe56_to_n_Co56_derived)*Y(P_nse)*state.rho;
    jac.set(Fe56, Fe56, scratch);

    scratch = screened_rates(k_Co56_to_Fe56) + screened_rates(k_n_Co56_to_p_nse_Fe56)*Y(N)*state.rho;
    jac.set(Fe56, Co56, scratch);

    scratch = screened_rates(k_Co57_to_p_nse_Fe56_derived);
    jac.set(Fe56, Co57, scratch);

    scratch = -screened_rates(k_n_Co55_to_Co56)*Y(Co55)*state.rho;
    jac.set(Co55, N, scratch);

    scratch = -screened_rates(k_p_nse_Co55_to_He4_Fe52_derived)*Y(Co55)*state.rho - screened_rates(k_p_nse_Co55_to_Ni56)*Y(Co55)*state.rho + screened_rates(k_p_nse_Fe54_to_Co55)*Y(Fe54)*state.rho + screened_rates(k_p_nse_Ni58_to_He4_Co55_derived)*Y(Ni58)*state.rho;
    jac.set(Co55, P_nse, scratch);

    scratch = -screened_rates(k_He4_Co55_to_p_nse_Ni58)*Y(Co55)*state.rho + screened_rates(k_He4_Fe52_to_p_nse_Co55)*Y(Fe52)*state.rho + screened_rates(k_He4_Mn51_to_Co55)*Y(Mn51)*state.rho;
    jac.set(Co55, He4, scratch);

    scratch = screened_rates(k_He4_Mn51_to_Co55)*Y(He4)*state.rho;
    jac.set(Co55, Mn51, scratch);

    scratch = screened_rates(k_He4_Fe52_to_p_nse_Co55)*Y(He4)*state.rho;
    jac.set(Co55, Fe52, scratch);

    scratch = screened_rates(k_p_nse_Fe54_to_Co55)*Y(P_nse)*state.rho;
    jac.set(Co55, Fe54, scratch);

    scratch = -screened_rates(k_Co55_to_He4_Mn51_derived) - screened_rates(k_Co55_to_p_nse_Fe54_derived) - screened_rates(k_He4_Co55_to_p_nse_Ni58)*Y(He4)*state.rho - screened_rates(k_n_Co55_to_Co56)*Y(N)*state.rho - screened_rates(k_p_nse_Co55_to_He4_Fe52_derived)*Y(P_nse)*state.rho - screened_rates(k_p_nse_Co55_to_Ni56)*Y(P_nse)*state.rho;
    jac.set(Co55, Co55, scratch);

    scratch = screened_rates(k_Co56_to_n_Co55_derived);
    jac.set(Co55, Co56, scratch);

    scratch = screened_rates(k_Ni56_to_p_nse_Co55_derived);
    jac.set(Co55, Ni56, scratch);

    scratch = screened_rates(k_p_nse_Ni58_to_He4_Co55_derived)*Y(P_nse)*state.rho;
    jac.set(Co55, Ni58, scratch);

    scratch = screened_rates(k_n_Co55_to_Co56)*Y(Co55)*state.rho - screened_rates(k_n_Co56_to_Co57)*Y(Co56)*state.rho - screened_rates(k_n_Co56_to_p_nse_Fe56)*Y(Co56)*state.rho + screened_rates(k_n_Ni56_to_p_nse_Co56)*Y(Ni56)*state.rho;
    jac.set(Co56, N, scratch);

    scratch = -screened_rates(k_p_nse_Co56_to_n_Ni56_derived)*Y(Co56)*state.rho + screened_rates(k_p_nse_Fe56_to_n_Co56_derived)*Y(Fe56)*state.rho;
    jac.set(Co56, P_nse, scratch);

    scratch = screened_rates(k_Fe56_to_Co56) + screened_rates(k_p_nse_Fe56_to_n_Co56_derived)*Y(P_nse)*state.rho;
    jac.set(Co56, Fe56, scratch);

    scratch = screened_rates(k_n_Co55_to_Co56)*Y(N)*state.rho;
    jac.set(Co56, Co55, scratch);

    scratch = -screened_rates(k_Co56_to_Fe56) - screened_rates(k_Co56_to_Ni56) - screened_rates(k_Co56_to_n_Co55_derived) - screened_rates(k_n_Co56_to_Co57)*Y(N)*state.rho - screened_rates(k_n_Co56_to_p_nse_Fe56)*Y(N)*state.rho - screened_rates(k_p_nse_Co56_to_n_Ni56_derived)*Y(P_nse)*state.rho;
    jac.set(Co56, Co56, scratch);

    scratch = screened_rates(k_Co57_to_n_Co56_derived);
    jac.set(Co56, Co57, scratch);

    scratch = screened_rates(k_Ni56_to_Co56) + screened_rates(k_n_Ni56_to_p_nse_Co56)*Y(N)*state.rho;
    jac.set(Co56, Ni56, scratch);

    scratch = screened_rates(k_n_Co56_to_Co57)*Y(Co56)*state.rho;
    jac.set(Co57, N, scratch);

    scratch = -screened_rates(k_p_nse_Co57_to_He4_Fe54)*Y(Co57)*state.rho - screened_rates(k_p_nse_Co57_to_Ni58)*Y(Co57)*state.rho + screened_rates(k_p_nse_Fe56_to_Co57)*Y(Fe56)*state.rho;
    jac.set(Co57, P_nse, scratch);

    scratch = screened_rates(k_He4_Fe54_to_p_nse_Co57_derived)*Y(Fe54)*state.rho;
    jac.set(Co57, He4, scratch);

    scratch = screened_rates(k_He4_Fe54_to_p_nse_Co57_derived)*Y(He4)*state.rho;
    jac.set(Co57, Fe54, scratch);

    scratch = screened_rates(k_p_nse_Fe56_to_Co57)*Y(P_nse)*state.rho;
    jac.set(Co57, Fe56, scratch);

    scratch = screened_rates(k_n_Co56_to_Co57)*Y(N)*state.rho;
    jac.set(Co57, Co56, scratch);

    scratch = -screened_rates(k_Co57_to_n_Co56_derived) - screened_rates(k_Co57_to_p_nse_Fe56_derived) - screened_rates(k_p_nse_Co57_to_He4_Fe54)*Y(P_nse)*state.rho - screened_rates(k_p_nse_Co57_to_Ni58)*Y(P_nse)*state.rho;
    jac.set(Co57, Co57, scratch);

    scratch = screened_rates(k_Ni58_to_p_nse_Co57_derived);
    jac.set(Co57, Ni58, scratch);

    scratch = -2.0*screened_rates(k_Ni56_n_n_to_Ni58_approx)*Y(Ni56)*Y(N)*amrex::Math::powi<2>(state.rho) - screened_rates(k_n_Ni56_to_p_nse_Co56)*Y(Ni56)*state.rho;
    jac.set(Ni56, N, scratch);

    scratch = screened_rates(k_p_nse_Co55_to_Ni56)*Y(Co55)*state.rho + screened_rates(k_p_nse_Co56_to_n_Ni56_derived)*Y(Co56)*state.rho;
    jac.set(Ni56, P_nse, scratch);

    scratch = screened_rates(k_He4_Fe52_to_Ni56)*Y(Fe52)*state.rho;
    jac.set(Ni56, He4, scratch);

    scratch = screened_rates(k_He4_Fe52_to_Ni56)*Y(He4)*state.rho;
    jac.set(Ni56, Fe52, scratch);

    scratch = screened_rates(k_p_nse_Co55_to_Ni56)*Y(P_nse)*state.rho;
    jac.set(Ni56, Co55, scratch);

    scratch = screened_rates(k_Co56_to_Ni56) + screened_rates(k_p_nse_Co56_to_n_Ni56_derived)*Y(P_nse)*state.rho;
    jac.set(Ni56, Co56, scratch);

    scratch = -screened_rates(k_Ni56_n_n_to_Ni58_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) - screened_rates(k_Ni56_to_Co56) - screened_rates(k_Ni56_to_He4_Fe52_derived) - screened_rates(k_Ni56_to_p_nse_Co55_derived) - screened_rates(k_n_Ni56_to_p_nse_Co56)*Y(N)*state.rho;
    jac.set(Ni56, Ni56, scratch);

    scratch = screened_rates(k_Ni58_to_Ni56_n_n_approx);
    jac.set(Ni56, Ni58, scratch);

    scratch = 2.0*screened_rates(k_Ni56_n_n_to_Ni58_approx)*Y(Ni56)*Y(N)*amrex::Math::powi<2>(state.rho);
    jac.set(Ni58, N, scratch);

    scratch = screened_rates(k_p_nse_Co57_to_Ni58)*Y(Co57)*state.rho - screened_rates(k_p_nse_Ni58_to_He4_Co55_derived)*Y(Ni58)*state.rho;
    jac.set(Ni58, P_nse, scratch);

    scratch = screened_rates(k_He4_Co55_to_p_nse_Ni58)*Y(Co55)*state.rho + screened_rates(k_He4_Fe54_to_Ni58)*Y(Fe54)*state.rho;
    jac.set(Ni58, He4, scratch);

    scratch = screened_rates(k_He4_Fe54_to_Ni58)*Y(He4)*state.rho;
    jac.set(Ni58, Fe54, scratch);

    scratch = screened_rates(k_He4_Co55_to_p_nse_Ni58)*Y(He4)*state.rho;
    jac.set(Ni58, Co55, scratch);

    scratch = screened_rates(k_p_nse_Co57_to_Ni58)*Y(P_nse)*state.rho;
    jac.set(Ni58, Co57, scratch);

    scratch = screened_rates(k_Ni56_n_n_to_Ni58_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(Ni58, Ni56, scratch);

    scratch = -screened_rates(k_Ni58_to_He4_Fe54_derived) - screened_rates(k_Ni58_to_Ni56_n_n_approx) - screened_rates(k_Ni58_to_p_nse_Co57_derived) - screened_rates(k_p_nse_Ni58_to_He4_Co55_derived)*Y(P_nse)*state.rho;
    jac.set(Ni58, Ni58, scratch);


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

    amrex::Real dsneutdt{};
#ifdef NEUTRINOS
    amrex::Real sneut, dsneutdd, dsnuda, dsnudz;
    constexpr int do_derivatives{1};
    sneut5<do_derivatives>(state.T, state.rho, state.abar, state.zbar, sneut, dsneutdt, dsneutdd, dsnuda, dsnudz);

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
