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
void ener_gener_rate(T const& dydt, Real& enuc)
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

    Array1D<Real, 1, NumSpec> Y;
    for (int n = 1; n <= NumSpec; ++n) {
        Y(n) = state.xn[n-1] * aion_inv[n-1];
    }

    [[maybe_unused]] Real rhoy = state.rho * state.y_e;

    // Calculate Reaclib rates

    plasma_state_t pstate{};
    fill_plasma_state(pstate, state.T, state.rho, Y);

    tf_t tfactors = evaluate_tfactors(state.T);

    fill_reaclib_rates<do_T_derivatives, T>(tfactors, rate_eval);



    // Evaluate screening factors

    Real ratraw, dratraw_dT;
    Real scor, dscor_dt;
    [[maybe_unused]] Real scor2, dscor2_dt;


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 6.0_rt, 12.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_c12_to_n13);
    rate_eval.screened_rates(k_p_c12_to_n13) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_c12_to_n13);
        rate_eval.dscreened_rates_dT(k_p_c12_to_n13) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 6.0_rt, 12.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_c12_to_o16);
    rate_eval.screened_rates(k_he4_c12_to_o16) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_c12_to_o16);
        rate_eval.dscreened_rates_dT(k_he4_c12_to_o16) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 7.0_rt, 14.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_n14_to_f18);
    rate_eval.screened_rates(k_he4_n14_to_f18) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_n14_to_f18);
        rate_eval.dscreened_rates_dT(k_he4_n14_to_f18) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 8.0_rt, 16.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_o16_to_ne20);
    rate_eval.screened_rates(k_he4_o16_to_ne20) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_o16_to_ne20);
        rate_eval.dscreened_rates_dT(k_he4_o16_to_ne20) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 9.0_rt, 18.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_f18_to_na22);
    rate_eval.screened_rates(k_he4_f18_to_na22) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_f18_to_na22);
        rate_eval.dscreened_rates_dT(k_he4_f18_to_na22) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_he4_f18_to_p_ne21);
    rate_eval.screened_rates(k_he4_f18_to_p_ne21) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_f18_to_p_ne21);
        rate_eval.dscreened_rates_dT(k_he4_f18_to_p_ne21) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 10.0_rt, 20.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_ne20_to_mg24);
    rate_eval.screened_rates(k_he4_ne20_to_mg24) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_ne20_to_mg24);
        rate_eval.dscreened_rates_dT(k_he4_ne20_to_mg24) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_he4_ne20_to_p_na23_derived);
    rate_eval.screened_rates(k_he4_ne20_to_p_na23_derived) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_ne20_to_p_na23_derived);
        rate_eval.dscreened_rates_dT(k_he4_ne20_to_p_na23_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 10.0_rt, 21.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_ne21_to_na22);
    rate_eval.screened_rates(k_p_ne21_to_na22) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_ne21_to_na22);
        rate_eval.dscreened_rates_dT(k_p_ne21_to_na22) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_ne21_to_he4_f18_derived);
    rate_eval.screened_rates(k_p_ne21_to_he4_f18_derived) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_ne21_to_he4_f18_derived);
        rate_eval.dscreened_rates_dT(k_p_ne21_to_he4_f18_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 11.0_rt, 23.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_na23_to_mg24);
    rate_eval.screened_rates(k_p_na23_to_mg24) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_na23_to_mg24);
        rate_eval.dscreened_rates_dT(k_p_na23_to_mg24) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_na23_to_he4_ne20);
    rate_eval.screened_rates(k_p_na23_to_he4_ne20) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_na23_to_he4_ne20);
        rate_eval.dscreened_rates_dT(k_p_na23_to_he4_ne20) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 12.0_rt, 24.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_mg24_to_si28);
    rate_eval.screened_rates(k_he4_mg24_to_si28) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_mg24_to_si28);
        rate_eval.dscreened_rates_dT(k_he4_mg24_to_si28) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_he4_mg24_to_p_al27_derived);
    rate_eval.screened_rates(k_he4_mg24_to_p_al27_derived) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_mg24_to_p_al27_derived);
        rate_eval.dscreened_rates_dT(k_he4_mg24_to_p_al27_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 13.0_rt, 27.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_al27_to_si28);
    rate_eval.screened_rates(k_p_al27_to_si28) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_al27_to_si28);
        rate_eval.dscreened_rates_dT(k_p_al27_to_si28) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_al27_to_he4_mg24);
    rate_eval.screened_rates(k_p_al27_to_he4_mg24) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_al27_to_he4_mg24);
        rate_eval.dscreened_rates_dT(k_p_al27_to_he4_mg24) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 14.0_rt, 28.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_si28_to_s32);
    rate_eval.screened_rates(k_he4_si28_to_s32) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_si28_to_s32);
        rate_eval.dscreened_rates_dT(k_he4_si28_to_s32) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_he4_si28_to_p_p31_derived);
    rate_eval.screened_rates(k_he4_si28_to_p_p31_derived) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_si28_to_p_p31_derived);
        rate_eval.dscreened_rates_dT(k_he4_si28_to_p_p31_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 15.0_rt, 31.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_p31_to_s32);
    rate_eval.screened_rates(k_p_p31_to_s32) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_p31_to_s32);
        rate_eval.dscreened_rates_dT(k_p_p31_to_s32) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_p31_to_he4_si28);
    rate_eval.screened_rates(k_p_p31_to_he4_si28) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_p31_to_he4_si28);
        rate_eval.dscreened_rates_dT(k_p_p31_to_he4_si28) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 24.0_rt, 48.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_cr48_to_fe52);
    rate_eval.screened_rates(k_he4_cr48_to_fe52) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_cr48_to_fe52);
        rate_eval.dscreened_rates_dT(k_he4_cr48_to_fe52) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_he4_cr48_to_p_mn51);
    rate_eval.screened_rates(k_he4_cr48_to_p_mn51) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_cr48_to_p_mn51);
        rate_eval.dscreened_rates_dT(k_he4_cr48_to_p_mn51) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 25.0_rt, 51.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_mn51_to_fe52);
    rate_eval.screened_rates(k_p_mn51_to_fe52) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_mn51_to_fe52);
        rate_eval.dscreened_rates_dT(k_p_mn51_to_fe52) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_mn51_to_he4_cr48_derived);
    rate_eval.screened_rates(k_p_mn51_to_he4_cr48_derived) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_mn51_to_he4_cr48_derived);
        rate_eval.dscreened_rates_dT(k_p_mn51_to_he4_cr48_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 25.0_rt, 51.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_mn51_to_co55);
    rate_eval.screened_rates(k_he4_mn51_to_co55) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_mn51_to_co55);
        rate_eval.dscreened_rates_dT(k_he4_mn51_to_co55) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_he4_mn51_to_p_fe54);
    rate_eval.screened_rates(k_he4_mn51_to_p_fe54) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_mn51_to_p_fe54);
        rate_eval.dscreened_rates_dT(k_he4_mn51_to_p_fe54) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 26.0_rt, 52.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_fe52_to_ni56);
    rate_eval.screened_rates(k_he4_fe52_to_ni56) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_fe52_to_ni56);
        rate_eval.dscreened_rates_dT(k_he4_fe52_to_ni56) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_he4_fe52_to_p_co55);
    rate_eval.screened_rates(k_he4_fe52_to_p_co55) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_fe52_to_p_co55);
        rate_eval.dscreened_rates_dT(k_he4_fe52_to_p_co55) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 27.0_rt, 55.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_co55_to_ni56);
    rate_eval.screened_rates(k_p_co55_to_ni56) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_co55_to_ni56);
        rate_eval.dscreened_rates_dT(k_p_co55_to_ni56) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_co55_to_he4_fe52_derived);
    rate_eval.screened_rates(k_p_co55_to_he4_fe52_derived) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_co55_to_he4_fe52_derived);
        rate_eval.dscreened_rates_dT(k_p_co55_to_he4_fe52_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 27.0_rt, 55.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_co55_to_cu59);
    rate_eval.screened_rates(k_he4_co55_to_cu59) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_co55_to_cu59);
        rate_eval.dscreened_rates_dT(k_he4_co55_to_cu59) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_he4_co55_to_p_ni58);
    rate_eval.screened_rates(k_he4_co55_to_p_ni58) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_co55_to_p_ni58);
        rate_eval.dscreened_rates_dT(k_he4_co55_to_p_ni58) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 28.0_rt, 56.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_ni56_to_zn60);
    rate_eval.screened_rates(k_he4_ni56_to_zn60) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_ni56_to_zn60);
        rate_eval.dscreened_rates_dT(k_he4_ni56_to_zn60) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_he4_ni56_to_p_cu59_derived);
    rate_eval.screened_rates(k_he4_ni56_to_p_cu59_derived) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_ni56_to_p_cu59_derived);
        rate_eval.dscreened_rates_dT(k_he4_ni56_to_p_cu59_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 29.0_rt, 59.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_cu59_to_zn60);
    rate_eval.screened_rates(k_p_cu59_to_zn60) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_cu59_to_zn60);
        rate_eval.dscreened_rates_dT(k_p_cu59_to_zn60) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_cu59_to_he4_ni56);
    rate_eval.screened_rates(k_p_cu59_to_he4_ni56) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_cu59_to_he4_ni56);
        rate_eval.dscreened_rates_dT(k_p_cu59_to_he4_ni56) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(6.0_rt, 12.0_rt, 6.0_rt, 12.0_rt);


        static_assert(scn_fac.z1 == 6.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_c12_c12_to_p_na23);
    rate_eval.screened_rates(k_c12_c12_to_p_na23) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_c12_c12_to_p_na23);
        rate_eval.dscreened_rates_dT(k_c12_c12_to_p_na23) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_c12_c12_to_he4_ne20);
    rate_eval.screened_rates(k_c12_c12_to_he4_ne20) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_c12_c12_to_he4_ne20);
        rate_eval.dscreened_rates_dT(k_c12_c12_to_he4_ne20) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_c12_c12_to_mg24_modified);
    rate_eval.screened_rates(k_c12_c12_to_mg24_modified) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_c12_c12_to_mg24_modified);
        rate_eval.dscreened_rates_dT(k_c12_c12_to_mg24_modified) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 7.0_rt, 13.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_n13_to_p_o16);
    rate_eval.screened_rates(k_he4_n13_to_p_o16) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_n13_to_p_o16);
        rate_eval.dscreened_rates_dT(k_he4_n13_to_p_o16) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(6.0_rt, 12.0_rt, 8.0_rt, 16.0_rt);


        static_assert(scn_fac.z1 == 6.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_c12_o16_to_p_al27);
    rate_eval.screened_rates(k_c12_o16_to_p_al27) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_c12_o16_to_p_al27);
        rate_eval.dscreened_rates_dT(k_c12_o16_to_p_al27) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_c12_o16_to_he4_mg24);
    rate_eval.screened_rates(k_c12_o16_to_he4_mg24) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_c12_o16_to_he4_mg24);
        rate_eval.dscreened_rates_dT(k_c12_o16_to_he4_mg24) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_c12_o16_to_si28_modified);
    rate_eval.screened_rates(k_c12_o16_to_si28_modified) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_c12_o16_to_si28_modified);
        rate_eval.dscreened_rates_dT(k_c12_o16_to_si28_modified) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(8.0_rt, 16.0_rt, 8.0_rt, 16.0_rt);


        static_assert(scn_fac.z1 == 8.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_o16_o16_to_p_p31);
    rate_eval.screened_rates(k_o16_o16_to_p_p31) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_o16_o16_to_p_p31);
        rate_eval.dscreened_rates_dT(k_o16_o16_to_p_p31) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_o16_o16_to_he4_si28);
    rate_eval.screened_rates(k_o16_o16_to_he4_si28) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_o16_o16_to_he4_si28);
        rate_eval.dscreened_rates_dT(k_o16_o16_to_he4_si28) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_o16_o16_to_s32_modified);
    rate_eval.screened_rates(k_o16_o16_to_s32_modified) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_o16_o16_to_s32_modified);
        rate_eval.dscreened_rates_dT(k_o16_o16_to_s32_modified) = ratraw * dscor_dt + dratraw_dT * scor;
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


    ratraw = rate_eval.screened_rates(k_he4_he4_he4_to_c12);
    rate_eval.screened_rates(k_he4_he4_he4_to_c12) *= scor * scor2;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_he4_he4_to_c12);
        rate_eval.dscreened_rates_dT(k_he4_he4_he4_to_c12) = ratraw * (scor * dscor2_dt + dscor_dt * scor2) + dratraw_dT * scor * scor2;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 25.0_rt, 55.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_mn55_to_fe56);
    rate_eval.screened_rates(k_p_mn55_to_fe56) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_mn55_to_fe56);
        rate_eval.dscreened_rates_dT(k_p_mn55_to_fe56) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_mn55_to_n_fe55_derived);
    rate_eval.screened_rates(k_p_mn55_to_n_fe55_derived) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_mn55_to_n_fe55_derived);
        rate_eval.dscreened_rates_dT(k_p_mn55_to_n_fe55_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 26.0_rt, 53.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_fe53_to_ni57);
    rate_eval.screened_rates(k_he4_fe53_to_ni57) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_fe53_to_ni57);
        rate_eval.dscreened_rates_dT(k_he4_fe53_to_ni57) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_he4_fe53_to_p_co56);
    rate_eval.screened_rates(k_he4_fe53_to_p_co56) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_fe53_to_p_co56);
        rate_eval.dscreened_rates_dT(k_he4_fe53_to_p_co56) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_he4_fe53_to_n_ni56_derived);
    rate_eval.screened_rates(k_he4_fe53_to_n_ni56_derived) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_fe53_to_n_ni56_derived);
        rate_eval.dscreened_rates_dT(k_he4_fe53_to_n_ni56_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 26.0_rt, 54.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_fe54_to_co55);
    rate_eval.screened_rates(k_p_fe54_to_co55) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_fe54_to_co55);
        rate_eval.dscreened_rates_dT(k_p_fe54_to_co55) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_fe54_to_he4_mn51_derived);
    rate_eval.screened_rates(k_p_fe54_to_he4_mn51_derived) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_fe54_to_he4_mn51_derived);
        rate_eval.dscreened_rates_dT(k_p_fe54_to_he4_mn51_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 26.0_rt, 54.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_fe54_to_ni58);
    rate_eval.screened_rates(k_he4_fe54_to_ni58) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_fe54_to_ni58);
        rate_eval.dscreened_rates_dT(k_he4_fe54_to_ni58) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_he4_fe54_to_n_ni57_derived);
    rate_eval.screened_rates(k_he4_fe54_to_n_ni57_derived) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_fe54_to_n_ni57_derived);
        rate_eval.dscreened_rates_dT(k_he4_fe54_to_n_ni57_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_he4_fe54_to_p_co57_derived);
    rate_eval.screened_rates(k_he4_fe54_to_p_co57_derived) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_fe54_to_p_co57_derived);
        rate_eval.dscreened_rates_dT(k_he4_fe54_to_p_co57_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 26.0_rt, 55.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_fe55_to_co56);
    rate_eval.screened_rates(k_p_fe55_to_co56) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_fe55_to_co56);
        rate_eval.dscreened_rates_dT(k_p_fe55_to_co56) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_fe55_to_n_co55_derived);
    rate_eval.screened_rates(k_p_fe55_to_n_co55_derived) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_fe55_to_n_co55_derived);
        rate_eval.dscreened_rates_dT(k_p_fe55_to_n_co55_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 26.0_rt, 56.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_fe56_to_co57);
    rate_eval.screened_rates(k_p_fe56_to_co57) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_fe56_to_co57);
        rate_eval.dscreened_rates_dT(k_p_fe56_to_co57) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_fe56_to_n_co56_derived);
    rate_eval.screened_rates(k_p_fe56_to_n_co56_derived) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_fe56_to_n_co56_derived);
        rate_eval.dscreened_rates_dT(k_p_fe56_to_n_co56_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 27.0_rt, 56.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_co56_to_ni57);
    rate_eval.screened_rates(k_p_co56_to_ni57) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_co56_to_ni57);
        rate_eval.dscreened_rates_dT(k_p_co56_to_ni57) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_co56_to_n_ni56_derived);
    rate_eval.screened_rates(k_p_co56_to_n_ni56_derived) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_co56_to_n_ni56_derived);
        rate_eval.dscreened_rates_dT(k_p_co56_to_n_ni56_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_co56_to_he4_fe53_derived);
    rate_eval.screened_rates(k_p_co56_to_he4_fe53_derived) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_co56_to_he4_fe53_derived);
        rate_eval.dscreened_rates_dT(k_p_co56_to_he4_fe53_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 27.0_rt, 57.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_co57_to_ni58);
    rate_eval.screened_rates(k_p_co57_to_ni58) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_co57_to_ni58);
        rate_eval.dscreened_rates_dT(k_p_co57_to_ni58) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_co57_to_he4_fe54);
    rate_eval.screened_rates(k_p_co57_to_he4_fe54) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_co57_to_he4_fe54);
        rate_eval.dscreened_rates_dT(k_p_co57_to_he4_fe54) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_co57_to_n_ni57_derived);
    rate_eval.screened_rates(k_p_co57_to_n_ni57_derived) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_co57_to_n_ni57_derived);
        rate_eval.dscreened_rates_dT(k_p_co57_to_n_ni57_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 28.0_rt, 58.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_ni58_to_cu59);
    rate_eval.screened_rates(k_p_ni58_to_cu59) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_ni58_to_cu59);
        rate_eval.dscreened_rates_dT(k_p_ni58_to_cu59) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_ni58_to_he4_co55_derived);
    rate_eval.screened_rates(k_p_ni58_to_he4_co55_derived) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_ni58_to_he4_co55_derived);
        rate_eval.dscreened_rates_dT(k_p_ni58_to_he4_co55_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 8.0_rt, 16.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_o16_to_he4_n13_derived);
    rate_eval.screened_rates(k_p_o16_to_he4_n13_derived) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_o16_to_he4_n13_derived);
        rate_eval.dscreened_rates_dT(k_p_o16_to_he4_n13_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 26.0_rt, 55.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_fe55_to_n_ni58_derived);
    rate_eval.screened_rates(k_he4_fe55_to_n_ni58_derived) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_fe55_to_n_ni58_derived);
        rate_eval.dscreened_rates_dT(k_he4_fe55_to_n_ni58_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 27.0_rt, 56.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_co56_to_n_cu59_derived);
    rate_eval.screened_rates(k_he4_co56_to_n_cu59_derived) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_co56_to_n_cu59_derived);
        rate_eval.dscreened_rates_dT(k_he4_co56_to_n_cu59_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 28.0_rt, 57.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_ni57_to_n_zn60_derived);
    rate_eval.screened_rates(k_he4_ni57_to_n_zn60_derived) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_ni57_to_n_zn60_derived);
        rate_eval.dscreened_rates_dT(k_he4_ni57_to_n_zn60_derived) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 16.0_rt, 32.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_s32_to_ar36_removed);
    rate_eval.screened_rates(k_he4_s32_to_ar36_removed) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_s32_to_ar36_removed);
        rate_eval.dscreened_rates_dT(k_he4_s32_to_ar36_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_he4_s32_to_p_cl35_derived_removed);
    rate_eval.screened_rates(k_he4_s32_to_p_cl35_derived_removed) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_s32_to_p_cl35_derived_removed);
        rate_eval.dscreened_rates_dT(k_he4_s32_to_p_cl35_derived_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 17.0_rt, 35.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_cl35_to_ar36_removed);
    rate_eval.screened_rates(k_p_cl35_to_ar36_removed) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_cl35_to_ar36_removed);
        rate_eval.dscreened_rates_dT(k_p_cl35_to_ar36_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_cl35_to_he4_s32_removed);
    rate_eval.screened_rates(k_p_cl35_to_he4_s32_removed) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_cl35_to_he4_s32_removed);
        rate_eval.dscreened_rates_dT(k_p_cl35_to_he4_s32_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 18.0_rt, 36.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_ar36_to_ca40_removed);
    rate_eval.screened_rates(k_he4_ar36_to_ca40_removed) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_ar36_to_ca40_removed);
        rate_eval.dscreened_rates_dT(k_he4_ar36_to_ca40_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_he4_ar36_to_p_k39_derived_removed);
    rate_eval.screened_rates(k_he4_ar36_to_p_k39_derived_removed) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_ar36_to_p_k39_derived_removed);
        rate_eval.dscreened_rates_dT(k_he4_ar36_to_p_k39_derived_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 19.0_rt, 39.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_k39_to_ca40_removed);
    rate_eval.screened_rates(k_p_k39_to_ca40_removed) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_k39_to_ca40_removed);
        rate_eval.dscreened_rates_dT(k_p_k39_to_ca40_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_k39_to_he4_ar36_removed);
    rate_eval.screened_rates(k_p_k39_to_he4_ar36_removed) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_k39_to_he4_ar36_removed);
        rate_eval.dscreened_rates_dT(k_p_k39_to_he4_ar36_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 20.0_rt, 40.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_ca40_to_ti44_removed);
    rate_eval.screened_rates(k_he4_ca40_to_ti44_removed) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_ca40_to_ti44_removed);
        rate_eval.dscreened_rates_dT(k_he4_ca40_to_ti44_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_he4_ca40_to_p_sc43_derived_removed);
    rate_eval.screened_rates(k_he4_ca40_to_p_sc43_derived_removed) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_ca40_to_p_sc43_derived_removed);
        rate_eval.dscreened_rates_dT(k_he4_ca40_to_p_sc43_derived_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 21.0_rt, 43.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_sc43_to_ti44_removed);
    rate_eval.screened_rates(k_p_sc43_to_ti44_removed) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_sc43_to_ti44_removed);
        rate_eval.dscreened_rates_dT(k_p_sc43_to_ti44_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_sc43_to_he4_ca40_removed);
    rate_eval.screened_rates(k_p_sc43_to_he4_ca40_removed) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_sc43_to_he4_ca40_removed);
        rate_eval.dscreened_rates_dT(k_p_sc43_to_he4_ca40_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 22.0_rt, 44.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_ti44_to_cr48_removed);
    rate_eval.screened_rates(k_he4_ti44_to_cr48_removed) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_ti44_to_cr48_removed);
        rate_eval.dscreened_rates_dT(k_he4_ti44_to_cr48_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_he4_ti44_to_p_v47_removed);
    rate_eval.screened_rates(k_he4_ti44_to_p_v47_removed) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_ti44_to_p_v47_removed);
        rate_eval.dscreened_rates_dT(k_he4_ti44_to_p_v47_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 23.0_rt, 47.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_v47_to_cr48_removed);
    rate_eval.screened_rates(k_p_v47_to_cr48_removed) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_v47_to_cr48_removed);
        rate_eval.dscreened_rates_dT(k_p_v47_to_cr48_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_v47_to_he4_ti44_derived_removed);
    rate_eval.screened_rates(k_p_v47_to_he4_ti44_derived_removed) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_v47_to_he4_ti44_derived_removed);
        rate_eval.dscreened_rates_dT(k_p_v47_to_he4_ti44_derived_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    // Fill approximate rates

    fill_approx_rates<do_T_derivatives, T>(tfactors, rate_eval);

    // Calculate tabular rates

    [[maybe_unused]] Real rate, drate_dt, edot_nu, edot_gamma;

    tabular_evaluate(j_co55_fe55_meta, j_co55_fe55_rhoy, j_co55_fe55_temp, j_co55_fe55_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_co55_to_fe55) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_co55_to_fe55) = drate_dt;
    }
    rate_eval.add_energy_rate(k_co55_to_fe55) = edot_nu + edot_gamma;

    tabular_evaluate(j_co56_fe56_meta, j_co56_fe56_rhoy, j_co56_fe56_temp, j_co56_fe56_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_co56_to_fe56) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_co56_to_fe56) = drate_dt;
    }
    rate_eval.add_energy_rate(k_co56_to_fe56) = edot_nu + edot_gamma;

    tabular_evaluate(j_co56_ni56_meta, j_co56_ni56_rhoy, j_co56_ni56_temp, j_co56_ni56_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_co56_to_ni56) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_co56_to_ni56) = drate_dt;
    }
    rate_eval.add_energy_rate(k_co56_to_ni56) = edot_nu + edot_gamma;

    tabular_evaluate(j_co57_ni57_meta, j_co57_ni57_rhoy, j_co57_ni57_temp, j_co57_ni57_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_co57_to_ni57) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_co57_to_ni57) = drate_dt;
    }
    rate_eval.add_energy_rate(k_co57_to_ni57) = edot_nu + edot_gamma;

    tabular_evaluate(j_fe55_co55_meta, j_fe55_co55_rhoy, j_fe55_co55_temp, j_fe55_co55_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_fe55_to_co55) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_fe55_to_co55) = drate_dt;
    }
    rate_eval.add_energy_rate(k_fe55_to_co55) = edot_nu + edot_gamma;

    tabular_evaluate(j_fe55_mn55_meta, j_fe55_mn55_rhoy, j_fe55_mn55_temp, j_fe55_mn55_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_fe55_to_mn55) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_fe55_to_mn55) = drate_dt;
    }
    rate_eval.add_energy_rate(k_fe55_to_mn55) = edot_nu + edot_gamma;

    tabular_evaluate(j_fe56_co56_meta, j_fe56_co56_rhoy, j_fe56_co56_temp, j_fe56_co56_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_fe56_to_co56) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_fe56_to_co56) = drate_dt;
    }
    rate_eval.add_energy_rate(k_fe56_to_co56) = edot_nu + edot_gamma;

    tabular_evaluate(j_mn55_fe55_meta, j_mn55_fe55_rhoy, j_mn55_fe55_temp, j_mn55_fe55_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_mn55_to_fe55) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_mn55_to_fe55) = drate_dt;
    }
    rate_eval.add_energy_rate(k_mn55_to_fe55) = edot_nu + edot_gamma;

    tabular_evaluate(j_ni56_co56_meta, j_ni56_co56_rhoy, j_ni56_co56_temp, j_ni56_co56_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_ni56_to_co56) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_ni56_to_co56) = drate_dt;
    }
    rate_eval.add_energy_rate(k_ni56_to_co56) = edot_nu + edot_gamma;

    tabular_evaluate(j_ni57_co57_meta, j_ni57_co57_rhoy, j_ni57_co57_temp, j_ni57_co57_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_ni57_to_co57) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_ni57_to_co57) = drate_dt;
    }
    rate_eval.add_energy_rate(k_ni57_to_co57) = edot_nu + edot_gamma;


}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rhs_nuc(const burn_t& state,
             Array1D<Real, 1, neqs>& ydot_nuc,
             const Array1D<Real, 1, NumSpec>& Y,
             const Array1D<Real, 1, NumRates>& screened_rates) {

    using namespace Rates;

    ydot_nuc(N) =
        -screened_rates(k_n_to_p_weak_wc12)*Y(N) +
        (-screened_rates(k_n_fe52_to_fe53)*Y(Fe52)*Y(N)*state.rho + screened_rates(k_fe53_to_n_fe52_derived)*Y(Fe53)) +
        (-screened_rates(k_n_fe53_to_fe54)*Y(Fe53)*Y(N)*state.rho + screened_rates(k_fe54_to_n_fe53_derived)*Y(Fe54)) +
        (-screened_rates(k_n_fe54_to_fe55)*Y(Fe54)*Y(N)*state.rho + screened_rates(k_fe55_to_n_fe54_derived)*Y(Fe55)) +
        (-screened_rates(k_n_fe55_to_fe56)*Y(Fe55)*Y(N)*state.rho + screened_rates(k_fe56_to_n_fe55_derived)*Y(Fe56)) +
        (-screened_rates(k_n_co55_to_co56)*Y(Co55)*Y(N)*state.rho + screened_rates(k_co56_to_n_co55_derived)*Y(Co56)) +
        (-screened_rates(k_n_co56_to_co57)*Y(Co56)*Y(N)*state.rho + screened_rates(k_co57_to_n_co56_derived)*Y(Co57)) +
        (-screened_rates(k_n_ni56_to_ni57)*Y(N)*Y(Ni56)*state.rho + screened_rates(k_ni57_to_n_ni56_derived)*Y(Ni57)) +
        (-screened_rates(k_n_ni57_to_ni58)*Y(N)*Y(Ni57)*state.rho + screened_rates(k_ni58_to_n_ni57_derived)*Y(Ni58)) +
        (-screened_rates(k_n_fe55_to_p_mn55)*Y(Fe55)*Y(N)*state.rho + screened_rates(k_p_mn55_to_n_fe55_derived)*Y(Mn55)*Y(H1)*state.rho) +
        (-screened_rates(k_n_co55_to_p_fe55)*Y(Co55)*Y(N)*state.rho + screened_rates(k_p_fe55_to_n_co55_derived)*Y(Fe55)*Y(H1)*state.rho) +
        (-screened_rates(k_n_co56_to_p_fe56)*Y(Co56)*Y(N)*state.rho + screened_rates(k_p_fe56_to_n_co56_derived)*Y(Fe56)*Y(H1)*state.rho) +
        (-screened_rates(k_n_ni56_to_p_co56)*Y(N)*Y(Ni56)*state.rho + screened_rates(k_p_co56_to_n_ni56_derived)*Y(Co56)*Y(H1)*state.rho) +
        (-screened_rates(k_n_ni56_to_he4_fe53)*Y(N)*Y(Ni56)*state.rho + screened_rates(k_he4_fe53_to_n_ni56_derived)*Y(Fe53)*Y(He4)*state.rho) +
        (-screened_rates(k_n_ni57_to_p_co57)*Y(N)*Y(Ni57)*state.rho + screened_rates(k_p_co57_to_n_ni57_derived)*Y(Co57)*Y(H1)*state.rho) +
        (-screened_rates(k_n_ni57_to_he4_fe54)*Y(N)*Y(Ni57)*state.rho + screened_rates(k_he4_fe54_to_n_ni57_derived)*Y(Fe54)*Y(He4)*state.rho) +
        (-screened_rates(k_n_ni58_to_he4_fe55)*Y(N)*Y(Ni58)*state.rho + screened_rates(k_he4_fe55_to_n_ni58_derived)*Y(Fe55)*Y(He4)*state.rho) +
        (-screened_rates(k_n_cu59_to_he4_co56)*Y(Cu59)*Y(N)*state.rho + screened_rates(k_he4_co56_to_n_cu59_derived)*Y(Co56)*Y(He4)*state.rho) +
        (-screened_rates(k_n_zn60_to_he4_ni57)*Y(N)*Y(Zn60)*state.rho + screened_rates(k_he4_ni57_to_n_zn60_derived)*Y(He4)*Y(Ni57)*state.rho);

    ydot_nuc(H1) =
        (-screened_rates(k_p_c12_to_n13)*Y(C12)*Y(H1)*state.rho + screened_rates(k_n13_to_p_c12_derived)*Y(N13)) +
        (-screened_rates(k_p_ne21_to_na22)*Y(Ne21)*Y(H1)*state.rho + screened_rates(k_na22_to_p_ne21_derived)*Y(Na22)) +
        (-screened_rates(k_p_na23_to_mg24)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_mg24_to_p_na23_derived)*Y(Mg24)) +
        (-screened_rates(k_p_al27_to_si28)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_si28_to_p_al27_derived)*Y(Si28)) +
        (-screened_rates(k_p_p31_to_s32)*Y(P31)*Y(H1)*state.rho + screened_rates(k_s32_to_p_p31_derived)*Y(S32)) +
        (-screened_rates(k_p_mn51_to_fe52)*Y(Mn51)*Y(H1)*state.rho + screened_rates(k_fe52_to_p_mn51_derived)*Y(Fe52)) +
        (-screened_rates(k_p_co55_to_ni56)*Y(Co55)*Y(H1)*state.rho + screened_rates(k_ni56_to_p_co55_derived)*Y(Ni56)) +
        (-screened_rates(k_p_cu59_to_zn60)*Y(Cu59)*Y(H1)*state.rho + screened_rates(k_zn60_to_p_cu59_derived)*Y(Zn60)) +
        0.5*screened_rates(k_c12_c12_to_p_na23)*std::pow(Y(C12), 2)*state.rho +
        (screened_rates(k_he4_n13_to_p_o16)*Y(He4)*Y(N13)*state.rho + -screened_rates(k_p_o16_to_he4_n13_derived)*Y(O16)*Y(H1)*state.rho) +
        screened_rates(k_c12_o16_to_p_al27)*Y(C12)*Y(O16)*state.rho +
        0.5*screened_rates(k_o16_o16_to_p_p31)*std::pow(Y(O16), 2)*state.rho +
        (screened_rates(k_he4_f18_to_p_ne21)*Y(F18)*Y(He4)*state.rho + -screened_rates(k_p_ne21_to_he4_f18_derived)*Y(Ne21)*Y(H1)*state.rho) +
        (-screened_rates(k_p_na23_to_he4_ne20)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_he4_ne20_to_p_na23_derived)*Y(He4)*Y(Ne20)*state.rho) +
        (-screened_rates(k_p_al27_to_he4_mg24)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_he4_mg24_to_p_al27_derived)*Y(He4)*Y(Mg24)*state.rho) +
        (-screened_rates(k_p_p31_to_he4_si28)*Y(P31)*Y(H1)*state.rho + screened_rates(k_he4_si28_to_p_p31_derived)*Y(He4)*Y(Si28)*state.rho) +
        (screened_rates(k_he4_cr48_to_p_mn51)*Y(Cr48)*Y(He4)*state.rho + -screened_rates(k_p_mn51_to_he4_cr48_derived)*Y(Mn51)*Y(H1)*state.rho) +
        (screened_rates(k_he4_fe52_to_p_co55)*Y(Fe52)*Y(He4)*state.rho + -screened_rates(k_p_co55_to_he4_fe52_derived)*Y(Co55)*Y(H1)*state.rho) +
        (-screened_rates(k_p_cu59_to_he4_ni56)*Y(Cu59)*Y(H1)*state.rho + screened_rates(k_he4_ni56_to_p_cu59_derived)*Y(He4)*Y(Ni56)*state.rho) +
        screened_rates(k_n_to_p_weak_wc12)*Y(N) +
        (-screened_rates(k_p_mn55_to_fe56)*Y(Mn55)*Y(H1)*state.rho + screened_rates(k_fe56_to_p_mn55_derived)*Y(Fe56)) +
        (-screened_rates(k_p_fe54_to_co55)*Y(Fe54)*Y(H1)*state.rho + screened_rates(k_co55_to_p_fe54_derived)*Y(Co55)) +
        (-screened_rates(k_p_fe55_to_co56)*Y(Fe55)*Y(H1)*state.rho + screened_rates(k_co56_to_p_fe55_derived)*Y(Co56)) +
        (-screened_rates(k_p_fe56_to_co57)*Y(Fe56)*Y(H1)*state.rho + screened_rates(k_co57_to_p_fe56_derived)*Y(Co57)) +
        (-screened_rates(k_p_co56_to_ni57)*Y(Co56)*Y(H1)*state.rho + screened_rates(k_ni57_to_p_co56_derived)*Y(Ni57)) +
        (-screened_rates(k_p_co57_to_ni58)*Y(Co57)*Y(H1)*state.rho + screened_rates(k_ni58_to_p_co57_derived)*Y(Ni58)) +
        (-screened_rates(k_p_ni58_to_cu59)*Y(Ni58)*Y(H1)*state.rho + screened_rates(k_cu59_to_p_ni58_derived)*Y(Cu59)) +
        (screened_rates(k_he4_mn51_to_p_fe54)*Y(He4)*Y(Mn51)*state.rho + -screened_rates(k_p_fe54_to_he4_mn51_derived)*Y(Fe54)*Y(H1)*state.rho) +
        (screened_rates(k_he4_fe53_to_p_co56)*Y(Fe53)*Y(He4)*state.rho + -screened_rates(k_p_co56_to_he4_fe53_derived)*Y(Co56)*Y(H1)*state.rho) +
        (screened_rates(k_n_fe55_to_p_mn55)*Y(Fe55)*Y(N)*state.rho + -screened_rates(k_p_mn55_to_n_fe55_derived)*Y(Mn55)*Y(H1)*state.rho) +
        (screened_rates(k_n_co55_to_p_fe55)*Y(Co55)*Y(N)*state.rho + -screened_rates(k_p_fe55_to_n_co55_derived)*Y(Fe55)*Y(H1)*state.rho) +
        (screened_rates(k_he4_co55_to_p_ni58)*Y(Co55)*Y(He4)*state.rho + -screened_rates(k_p_ni58_to_he4_co55_derived)*Y(Ni58)*Y(H1)*state.rho) +
        (screened_rates(k_n_co56_to_p_fe56)*Y(Co56)*Y(N)*state.rho + -screened_rates(k_p_fe56_to_n_co56_derived)*Y(Fe56)*Y(H1)*state.rho) +
        (-screened_rates(k_p_co57_to_he4_fe54)*Y(Co57)*Y(H1)*state.rho + screened_rates(k_he4_fe54_to_p_co57_derived)*Y(Fe54)*Y(He4)*state.rho) +
        (screened_rates(k_n_ni56_to_p_co56)*Y(N)*Y(Ni56)*state.rho + -screened_rates(k_p_co56_to_n_ni56_derived)*Y(Co56)*Y(H1)*state.rho) +
        (screened_rates(k_n_ni57_to_p_co57)*Y(N)*Y(Ni57)*state.rho + -screened_rates(k_p_co57_to_n_ni57_derived)*Y(Co57)*Y(H1)*state.rho);

    ydot_nuc(He4) =
        (-screened_rates(k_he4_c12_to_o16)*Y(C12)*Y(He4)*state.rho + screened_rates(k_o16_to_he4_c12_derived)*Y(O16)) +
        (-screened_rates(k_he4_n14_to_f18)*Y(He4)*Y(N14)*state.rho + screened_rates(k_f18_to_he4_n14_derived)*Y(F18)) +
        (-screened_rates(k_he4_o16_to_ne20)*Y(He4)*Y(O16)*state.rho + screened_rates(k_ne20_to_he4_o16_derived)*Y(Ne20)) +
        (-screened_rates(k_he4_f18_to_na22)*Y(F18)*Y(He4)*state.rho + screened_rates(k_na22_to_he4_f18_derived)*Y(Na22)) +
        (-screened_rates(k_he4_ne20_to_mg24)*Y(He4)*Y(Ne20)*state.rho + screened_rates(k_mg24_to_he4_ne20_derived)*Y(Mg24)) +
        (-screened_rates(k_he4_mg24_to_si28)*Y(He4)*Y(Mg24)*state.rho + screened_rates(k_si28_to_he4_mg24_derived)*Y(Si28)) +
        (-screened_rates(k_he4_si28_to_s32)*Y(He4)*Y(Si28)*state.rho + screened_rates(k_s32_to_he4_si28_derived)*Y(S32)) +
        (-screened_rates(k_he4_cr48_to_fe52)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_fe52_to_he4_cr48_derived)*Y(Fe52)) +
        (-screened_rates(k_he4_mn51_to_co55)*Y(He4)*Y(Mn51)*state.rho + screened_rates(k_co55_to_he4_mn51_derived)*Y(Co55)) +
        (-screened_rates(k_he4_fe52_to_ni56)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_ni56_to_he4_fe52_derived)*Y(Ni56)) +
        (-screened_rates(k_he4_co55_to_cu59)*Y(Co55)*Y(He4)*state.rho + screened_rates(k_cu59_to_he4_co55_derived)*Y(Cu59)) +
        (-screened_rates(k_he4_ni56_to_zn60)*Y(He4)*Y(Ni56)*state.rho + screened_rates(k_zn60_to_he4_ni56_derived)*Y(Zn60)) +
        0.5*screened_rates(k_c12_c12_to_he4_ne20)*std::pow(Y(C12), 2)*state.rho +
        (-screened_rates(k_he4_n13_to_p_o16)*Y(He4)*Y(N13)*state.rho + screened_rates(k_p_o16_to_he4_n13_derived)*Y(O16)*Y(H1)*state.rho) +
        screened_rates(k_c12_o16_to_he4_mg24)*Y(C12)*Y(O16)*state.rho +
        0.5*screened_rates(k_o16_o16_to_he4_si28)*std::pow(Y(O16), 2)*state.rho +
        (-screened_rates(k_he4_f18_to_p_ne21)*Y(F18)*Y(He4)*state.rho + screened_rates(k_p_ne21_to_he4_f18_derived)*Y(Ne21)*Y(H1)*state.rho) +
        (screened_rates(k_p_na23_to_he4_ne20)*Y(Na23)*Y(H1)*state.rho + -screened_rates(k_he4_ne20_to_p_na23_derived)*Y(He4)*Y(Ne20)*state.rho) +
        (screened_rates(k_p_al27_to_he4_mg24)*Y(Al27)*Y(H1)*state.rho + -screened_rates(k_he4_mg24_to_p_al27_derived)*Y(He4)*Y(Mg24)*state.rho) +
        (screened_rates(k_p_p31_to_he4_si28)*Y(P31)*Y(H1)*state.rho + -screened_rates(k_he4_si28_to_p_p31_derived)*Y(He4)*Y(Si28)*state.rho) +
        (-screened_rates(k_he4_cr48_to_p_mn51)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_p_mn51_to_he4_cr48_derived)*Y(Mn51)*Y(H1)*state.rho) +
        (-screened_rates(k_he4_fe52_to_p_co55)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_p_co55_to_he4_fe52_derived)*Y(Co55)*Y(H1)*state.rho) +
        (screened_rates(k_p_cu59_to_he4_ni56)*Y(Cu59)*Y(H1)*state.rho + -screened_rates(k_he4_ni56_to_p_cu59_derived)*Y(He4)*Y(Ni56)*state.rho) +
        (-0.5*screened_rates(k_he4_he4_he4_to_c12)*std::pow(Y(He4), 3)*std::pow(state.rho, 2) + 3.0*screened_rates(k_c12_to_he4_he4_he4_derived)*Y(C12)) +
        (-screened_rates(k_he4_fe53_to_ni57)*Y(Fe53)*Y(He4)*state.rho + screened_rates(k_ni57_to_he4_fe53_derived)*Y(Ni57)) +
        (-screened_rates(k_he4_fe54_to_ni58)*Y(Fe54)*Y(He4)*state.rho + screened_rates(k_ni58_to_he4_fe54_derived)*Y(Ni58)) +
        (-screened_rates(k_he4_mn51_to_p_fe54)*Y(He4)*Y(Mn51)*state.rho + screened_rates(k_p_fe54_to_he4_mn51_derived)*Y(Fe54)*Y(H1)*state.rho) +
        (-screened_rates(k_he4_fe53_to_p_co56)*Y(Fe53)*Y(He4)*state.rho + screened_rates(k_p_co56_to_he4_fe53_derived)*Y(Co56)*Y(H1)*state.rho) +
        (-screened_rates(k_he4_co55_to_p_ni58)*Y(Co55)*Y(He4)*state.rho + screened_rates(k_p_ni58_to_he4_co55_derived)*Y(Ni58)*Y(H1)*state.rho) +
        (screened_rates(k_p_co57_to_he4_fe54)*Y(Co57)*Y(H1)*state.rho + -screened_rates(k_he4_fe54_to_p_co57_derived)*Y(Fe54)*Y(He4)*state.rho) +
        (screened_rates(k_n_ni56_to_he4_fe53)*Y(N)*Y(Ni56)*state.rho + -screened_rates(k_he4_fe53_to_n_ni56_derived)*Y(Fe53)*Y(He4)*state.rho) +
        (screened_rates(k_n_ni57_to_he4_fe54)*Y(N)*Y(Ni57)*state.rho + -screened_rates(k_he4_fe54_to_n_ni57_derived)*Y(Fe54)*Y(He4)*state.rho) +
        (screened_rates(k_n_ni58_to_he4_fe55)*Y(N)*Y(Ni58)*state.rho + -screened_rates(k_he4_fe55_to_n_ni58_derived)*Y(Fe55)*Y(He4)*state.rho) +
        (screened_rates(k_n_cu59_to_he4_co56)*Y(Cu59)*Y(N)*state.rho + -screened_rates(k_he4_co56_to_n_cu59_derived)*Y(Co56)*Y(He4)*state.rho) +
        (screened_rates(k_n_zn60_to_he4_ni57)*Y(N)*Y(Zn60)*state.rho + -screened_rates(k_he4_ni57_to_n_zn60_derived)*Y(He4)*Y(Ni57)*state.rho) +
        (-screened_rates(k_s32_he4_to_ar36_approx)*Y(He4)*Y(S32)*state.rho + screened_rates(k_ar36_to_s32_he4_approx)*Y(Ar36)) +
        (-screened_rates(k_ar36_he4_to_ca40_approx)*Y(Ar36)*Y(He4)*state.rho + screened_rates(k_ca40_to_ar36_he4_approx)*Y(Ca40)) +
        (-screened_rates(k_ca40_he4_to_ti44_approx)*Y(Ca40)*Y(He4)*state.rho + screened_rates(k_ti44_to_ca40_he4_approx)*Y(Ti44)) +
        (-screened_rates(k_ti44_he4_to_cr48_approx)*Y(He4)*Y(Ti44)*state.rho + screened_rates(k_cr48_to_ti44_he4_approx)*Y(Cr48));

    ydot_nuc(C12) =
        (-screened_rates(k_p_c12_to_n13)*Y(C12)*Y(H1)*state.rho + screened_rates(k_n13_to_p_c12_derived)*Y(N13)) +
        (-screened_rates(k_he4_c12_to_o16)*Y(C12)*Y(He4)*state.rho + screened_rates(k_o16_to_he4_c12_derived)*Y(O16)) +
        -screened_rates(k_c12_c12_to_p_na23)*std::pow(Y(C12), 2)*state.rho +
        -screened_rates(k_c12_c12_to_he4_ne20)*std::pow(Y(C12), 2)*state.rho +
        -screened_rates(k_c12_o16_to_p_al27)*Y(C12)*Y(O16)*state.rho +
        -screened_rates(k_c12_o16_to_he4_mg24)*Y(C12)*Y(O16)*state.rho +
        (0.16666666666666667*screened_rates(k_he4_he4_he4_to_c12)*std::pow(Y(He4), 3)*std::pow(state.rho, 2) + -screened_rates(k_c12_to_he4_he4_he4_derived)*Y(C12)) +
        -screened_rates(k_c12_c12_to_mg24_modified)*std::pow(Y(C12), 2)*state.rho +
        -screened_rates(k_c12_o16_to_si28_modified)*Y(C12)*Y(O16)*state.rho;

    ydot_nuc(N13) =
        (screened_rates(k_p_c12_to_n13)*Y(C12)*Y(H1)*state.rho + -screened_rates(k_n13_to_p_c12_derived)*Y(N13)) +
        (-screened_rates(k_he4_n13_to_p_o16)*Y(He4)*Y(N13)*state.rho + screened_rates(k_p_o16_to_he4_n13_derived)*Y(O16)*Y(H1)*state.rho);

    ydot_nuc(N14) =
        (-screened_rates(k_he4_n14_to_f18)*Y(He4)*Y(N14)*state.rho + screened_rates(k_f18_to_he4_n14_derived)*Y(F18));

    ydot_nuc(O16) =
        (screened_rates(k_he4_c12_to_o16)*Y(C12)*Y(He4)*state.rho + -screened_rates(k_o16_to_he4_c12_derived)*Y(O16)) +
        (-screened_rates(k_he4_o16_to_ne20)*Y(He4)*Y(O16)*state.rho + screened_rates(k_ne20_to_he4_o16_derived)*Y(Ne20)) +
        (screened_rates(k_he4_n13_to_p_o16)*Y(He4)*Y(N13)*state.rho + -screened_rates(k_p_o16_to_he4_n13_derived)*Y(O16)*Y(H1)*state.rho) +
        -screened_rates(k_c12_o16_to_p_al27)*Y(C12)*Y(O16)*state.rho +
        -screened_rates(k_c12_o16_to_he4_mg24)*Y(C12)*Y(O16)*state.rho +
        -screened_rates(k_o16_o16_to_p_p31)*std::pow(Y(O16), 2)*state.rho +
        -screened_rates(k_o16_o16_to_he4_si28)*std::pow(Y(O16), 2)*state.rho +
        -screened_rates(k_o16_o16_to_s32_modified)*std::pow(Y(O16), 2)*state.rho +
        -screened_rates(k_c12_o16_to_si28_modified)*Y(C12)*Y(O16)*state.rho;

    ydot_nuc(F18) =
        (screened_rates(k_he4_n14_to_f18)*Y(He4)*Y(N14)*state.rho + -screened_rates(k_f18_to_he4_n14_derived)*Y(F18)) +
        (-screened_rates(k_he4_f18_to_na22)*Y(F18)*Y(He4)*state.rho + screened_rates(k_na22_to_he4_f18_derived)*Y(Na22)) +
        (-screened_rates(k_he4_f18_to_p_ne21)*Y(F18)*Y(He4)*state.rho + screened_rates(k_p_ne21_to_he4_f18_derived)*Y(Ne21)*Y(H1)*state.rho);

    ydot_nuc(Ne20) =
        (screened_rates(k_he4_o16_to_ne20)*Y(He4)*Y(O16)*state.rho + -screened_rates(k_ne20_to_he4_o16_derived)*Y(Ne20)) +
        (-screened_rates(k_he4_ne20_to_mg24)*Y(He4)*Y(Ne20)*state.rho + screened_rates(k_mg24_to_he4_ne20_derived)*Y(Mg24)) +
        0.5*screened_rates(k_c12_c12_to_he4_ne20)*std::pow(Y(C12), 2)*state.rho +
        (screened_rates(k_p_na23_to_he4_ne20)*Y(Na23)*Y(H1)*state.rho + -screened_rates(k_he4_ne20_to_p_na23_derived)*Y(He4)*Y(Ne20)*state.rho);

    ydot_nuc(Ne21) =
        (-screened_rates(k_p_ne21_to_na22)*Y(Ne21)*Y(H1)*state.rho + screened_rates(k_na22_to_p_ne21_derived)*Y(Na22)) +
        (screened_rates(k_he4_f18_to_p_ne21)*Y(F18)*Y(He4)*state.rho + -screened_rates(k_p_ne21_to_he4_f18_derived)*Y(Ne21)*Y(H1)*state.rho);

    ydot_nuc(Na22) =
        (screened_rates(k_he4_f18_to_na22)*Y(F18)*Y(He4)*state.rho + -screened_rates(k_na22_to_he4_f18_derived)*Y(Na22)) +
        (screened_rates(k_p_ne21_to_na22)*Y(Ne21)*Y(H1)*state.rho + -screened_rates(k_na22_to_p_ne21_derived)*Y(Na22));

    ydot_nuc(Na23) =
        (-screened_rates(k_p_na23_to_mg24)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_mg24_to_p_na23_derived)*Y(Mg24)) +
        0.5*screened_rates(k_c12_c12_to_p_na23)*std::pow(Y(C12), 2)*state.rho +
        (-screened_rates(k_p_na23_to_he4_ne20)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_he4_ne20_to_p_na23_derived)*Y(He4)*Y(Ne20)*state.rho);

    ydot_nuc(Mg24) =
        (screened_rates(k_he4_ne20_to_mg24)*Y(He4)*Y(Ne20)*state.rho + -screened_rates(k_mg24_to_he4_ne20_derived)*Y(Mg24)) +
        (screened_rates(k_p_na23_to_mg24)*Y(Na23)*Y(H1)*state.rho + -screened_rates(k_mg24_to_p_na23_derived)*Y(Mg24)) +
        (-screened_rates(k_he4_mg24_to_si28)*Y(He4)*Y(Mg24)*state.rho + screened_rates(k_si28_to_he4_mg24_derived)*Y(Si28)) +
        screened_rates(k_c12_o16_to_he4_mg24)*Y(C12)*Y(O16)*state.rho +
        (screened_rates(k_p_al27_to_he4_mg24)*Y(Al27)*Y(H1)*state.rho + -screened_rates(k_he4_mg24_to_p_al27_derived)*Y(He4)*Y(Mg24)*state.rho) +
        0.5*screened_rates(k_c12_c12_to_mg24_modified)*std::pow(Y(C12), 2)*state.rho;

    ydot_nuc(Al27) =
        (-screened_rates(k_p_al27_to_si28)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_si28_to_p_al27_derived)*Y(Si28)) +
        screened_rates(k_c12_o16_to_p_al27)*Y(C12)*Y(O16)*state.rho +
        (-screened_rates(k_p_al27_to_he4_mg24)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_he4_mg24_to_p_al27_derived)*Y(He4)*Y(Mg24)*state.rho);

    ydot_nuc(Si28) =
        (screened_rates(k_he4_mg24_to_si28)*Y(He4)*Y(Mg24)*state.rho + -screened_rates(k_si28_to_he4_mg24_derived)*Y(Si28)) +
        (screened_rates(k_p_al27_to_si28)*Y(Al27)*Y(H1)*state.rho + -screened_rates(k_si28_to_p_al27_derived)*Y(Si28)) +
        (-screened_rates(k_he4_si28_to_s32)*Y(He4)*Y(Si28)*state.rho + screened_rates(k_s32_to_he4_si28_derived)*Y(S32)) +
        0.5*screened_rates(k_o16_o16_to_he4_si28)*std::pow(Y(O16), 2)*state.rho +
        (screened_rates(k_p_p31_to_he4_si28)*Y(P31)*Y(H1)*state.rho + -screened_rates(k_he4_si28_to_p_p31_derived)*Y(He4)*Y(Si28)*state.rho) +
        screened_rates(k_c12_o16_to_si28_modified)*Y(C12)*Y(O16)*state.rho;

    ydot_nuc(P31) =
        (-screened_rates(k_p_p31_to_s32)*Y(P31)*Y(H1)*state.rho + screened_rates(k_s32_to_p_p31_derived)*Y(S32)) +
        0.5*screened_rates(k_o16_o16_to_p_p31)*std::pow(Y(O16), 2)*state.rho +
        (-screened_rates(k_p_p31_to_he4_si28)*Y(P31)*Y(H1)*state.rho + screened_rates(k_he4_si28_to_p_p31_derived)*Y(He4)*Y(Si28)*state.rho);

    ydot_nuc(S32) =
        (screened_rates(k_he4_si28_to_s32)*Y(He4)*Y(Si28)*state.rho + -screened_rates(k_s32_to_he4_si28_derived)*Y(S32)) +
        (screened_rates(k_p_p31_to_s32)*Y(P31)*Y(H1)*state.rho + -screened_rates(k_s32_to_p_p31_derived)*Y(S32)) +
        0.5*screened_rates(k_o16_o16_to_s32_modified)*std::pow(Y(O16), 2)*state.rho +
        (-screened_rates(k_s32_he4_to_ar36_approx)*Y(He4)*Y(S32)*state.rho + screened_rates(k_ar36_to_s32_he4_approx)*Y(Ar36));

    ydot_nuc(Ar36) =
        (screened_rates(k_s32_he4_to_ar36_approx)*Y(He4)*Y(S32)*state.rho + -screened_rates(k_ar36_to_s32_he4_approx)*Y(Ar36)) +
        (-screened_rates(k_ar36_he4_to_ca40_approx)*Y(Ar36)*Y(He4)*state.rho + screened_rates(k_ca40_to_ar36_he4_approx)*Y(Ca40));

    ydot_nuc(Ca40) =
        (screened_rates(k_ar36_he4_to_ca40_approx)*Y(Ar36)*Y(He4)*state.rho + -screened_rates(k_ca40_to_ar36_he4_approx)*Y(Ca40)) +
        (-screened_rates(k_ca40_he4_to_ti44_approx)*Y(Ca40)*Y(He4)*state.rho + screened_rates(k_ti44_to_ca40_he4_approx)*Y(Ti44));

    ydot_nuc(Ti44) =
        (screened_rates(k_ca40_he4_to_ti44_approx)*Y(Ca40)*Y(He4)*state.rho + -screened_rates(k_ti44_to_ca40_he4_approx)*Y(Ti44)) +
        (-screened_rates(k_ti44_he4_to_cr48_approx)*Y(He4)*Y(Ti44)*state.rho + screened_rates(k_cr48_to_ti44_he4_approx)*Y(Cr48));

    ydot_nuc(Cr48) =
        (-screened_rates(k_he4_cr48_to_fe52)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_fe52_to_he4_cr48_derived)*Y(Fe52)) +
        (-screened_rates(k_he4_cr48_to_p_mn51)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_p_mn51_to_he4_cr48_derived)*Y(Mn51)*Y(H1)*state.rho) +
        (screened_rates(k_ti44_he4_to_cr48_approx)*Y(He4)*Y(Ti44)*state.rho + -screened_rates(k_cr48_to_ti44_he4_approx)*Y(Cr48));

    ydot_nuc(Mn51) =
        (-screened_rates(k_p_mn51_to_fe52)*Y(Mn51)*Y(H1)*state.rho + screened_rates(k_fe52_to_p_mn51_derived)*Y(Fe52)) +
        (-screened_rates(k_he4_mn51_to_co55)*Y(He4)*Y(Mn51)*state.rho + screened_rates(k_co55_to_he4_mn51_derived)*Y(Co55)) +
        (screened_rates(k_he4_cr48_to_p_mn51)*Y(Cr48)*Y(He4)*state.rho + -screened_rates(k_p_mn51_to_he4_cr48_derived)*Y(Mn51)*Y(H1)*state.rho) +
        (-screened_rates(k_he4_mn51_to_p_fe54)*Y(He4)*Y(Mn51)*state.rho + screened_rates(k_p_fe54_to_he4_mn51_derived)*Y(Fe54)*Y(H1)*state.rho);

    ydot_nuc(Mn55) =
        (-screened_rates(k_p_mn55_to_fe56)*Y(Mn55)*Y(H1)*state.rho + screened_rates(k_fe56_to_p_mn55_derived)*Y(Fe56)) +
        (screened_rates(k_n_fe55_to_p_mn55)*Y(Fe55)*Y(N)*state.rho + -screened_rates(k_p_mn55_to_n_fe55_derived)*Y(Mn55)*Y(H1)*state.rho) +
        (screened_rates(k_fe55_to_mn55)*Y(Fe55) + -screened_rates(k_mn55_to_fe55)*Y(Mn55));

    ydot_nuc(Fe52) =
        (screened_rates(k_he4_cr48_to_fe52)*Y(Cr48)*Y(He4)*state.rho + -screened_rates(k_fe52_to_he4_cr48_derived)*Y(Fe52)) +
        (screened_rates(k_p_mn51_to_fe52)*Y(Mn51)*Y(H1)*state.rho + -screened_rates(k_fe52_to_p_mn51_derived)*Y(Fe52)) +
        (-screened_rates(k_he4_fe52_to_ni56)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_ni56_to_he4_fe52_derived)*Y(Ni56)) +
        (-screened_rates(k_he4_fe52_to_p_co55)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_p_co55_to_he4_fe52_derived)*Y(Co55)*Y(H1)*state.rho) +
        (-screened_rates(k_n_fe52_to_fe53)*Y(Fe52)*Y(N)*state.rho + screened_rates(k_fe53_to_n_fe52_derived)*Y(Fe53));

    ydot_nuc(Fe53) =
        (screened_rates(k_n_fe52_to_fe53)*Y(Fe52)*Y(N)*state.rho + -screened_rates(k_fe53_to_n_fe52_derived)*Y(Fe53)) +
        (-screened_rates(k_n_fe53_to_fe54)*Y(Fe53)*Y(N)*state.rho + screened_rates(k_fe54_to_n_fe53_derived)*Y(Fe54)) +
        (-screened_rates(k_he4_fe53_to_ni57)*Y(Fe53)*Y(He4)*state.rho + screened_rates(k_ni57_to_he4_fe53_derived)*Y(Ni57)) +
        (-screened_rates(k_he4_fe53_to_p_co56)*Y(Fe53)*Y(He4)*state.rho + screened_rates(k_p_co56_to_he4_fe53_derived)*Y(Co56)*Y(H1)*state.rho) +
        (screened_rates(k_n_ni56_to_he4_fe53)*Y(N)*Y(Ni56)*state.rho + -screened_rates(k_he4_fe53_to_n_ni56_derived)*Y(Fe53)*Y(He4)*state.rho);

    ydot_nuc(Fe54) =
        (screened_rates(k_n_fe53_to_fe54)*Y(Fe53)*Y(N)*state.rho + -screened_rates(k_fe54_to_n_fe53_derived)*Y(Fe54)) +
        (-screened_rates(k_n_fe54_to_fe55)*Y(Fe54)*Y(N)*state.rho + screened_rates(k_fe55_to_n_fe54_derived)*Y(Fe55)) +
        (-screened_rates(k_p_fe54_to_co55)*Y(Fe54)*Y(H1)*state.rho + screened_rates(k_co55_to_p_fe54_derived)*Y(Co55)) +
        (-screened_rates(k_he4_fe54_to_ni58)*Y(Fe54)*Y(He4)*state.rho + screened_rates(k_ni58_to_he4_fe54_derived)*Y(Ni58)) +
        (screened_rates(k_he4_mn51_to_p_fe54)*Y(He4)*Y(Mn51)*state.rho + -screened_rates(k_p_fe54_to_he4_mn51_derived)*Y(Fe54)*Y(H1)*state.rho) +
        (screened_rates(k_p_co57_to_he4_fe54)*Y(Co57)*Y(H1)*state.rho + -screened_rates(k_he4_fe54_to_p_co57_derived)*Y(Fe54)*Y(He4)*state.rho) +
        (screened_rates(k_n_ni57_to_he4_fe54)*Y(N)*Y(Ni57)*state.rho + -screened_rates(k_he4_fe54_to_n_ni57_derived)*Y(Fe54)*Y(He4)*state.rho);

    ydot_nuc(Fe55) =
        (screened_rates(k_n_fe54_to_fe55)*Y(Fe54)*Y(N)*state.rho + -screened_rates(k_fe55_to_n_fe54_derived)*Y(Fe55)) +
        (-screened_rates(k_n_fe55_to_fe56)*Y(Fe55)*Y(N)*state.rho + screened_rates(k_fe56_to_n_fe55_derived)*Y(Fe56)) +
        (-screened_rates(k_p_fe55_to_co56)*Y(Fe55)*Y(H1)*state.rho + screened_rates(k_co56_to_p_fe55_derived)*Y(Co56)) +
        (-screened_rates(k_n_fe55_to_p_mn55)*Y(Fe55)*Y(N)*state.rho + screened_rates(k_p_mn55_to_n_fe55_derived)*Y(Mn55)*Y(H1)*state.rho) +
        (screened_rates(k_n_co55_to_p_fe55)*Y(Co55)*Y(N)*state.rho + -screened_rates(k_p_fe55_to_n_co55_derived)*Y(Fe55)*Y(H1)*state.rho) +
        (screened_rates(k_n_ni58_to_he4_fe55)*Y(N)*Y(Ni58)*state.rho + -screened_rates(k_he4_fe55_to_n_ni58_derived)*Y(Fe55)*Y(He4)*state.rho) +
        (screened_rates(k_co55_to_fe55)*Y(Co55) + -screened_rates(k_fe55_to_co55)*Y(Fe55)) +
        (-screened_rates(k_fe55_to_mn55)*Y(Fe55) + screened_rates(k_mn55_to_fe55)*Y(Mn55));

    ydot_nuc(Fe56) =
        (screened_rates(k_p_mn55_to_fe56)*Y(Mn55)*Y(H1)*state.rho + -screened_rates(k_fe56_to_p_mn55_derived)*Y(Fe56)) +
        (screened_rates(k_n_fe55_to_fe56)*Y(Fe55)*Y(N)*state.rho + -screened_rates(k_fe56_to_n_fe55_derived)*Y(Fe56)) +
        (-screened_rates(k_p_fe56_to_co57)*Y(Fe56)*Y(H1)*state.rho + screened_rates(k_co57_to_p_fe56_derived)*Y(Co57)) +
        (screened_rates(k_n_co56_to_p_fe56)*Y(Co56)*Y(N)*state.rho + -screened_rates(k_p_fe56_to_n_co56_derived)*Y(Fe56)*Y(H1)*state.rho) +
        (screened_rates(k_co56_to_fe56)*Y(Co56) + -screened_rates(k_fe56_to_co56)*Y(Fe56));

    ydot_nuc(Co55) =
        (screened_rates(k_he4_mn51_to_co55)*Y(He4)*Y(Mn51)*state.rho + -screened_rates(k_co55_to_he4_mn51_derived)*Y(Co55)) +
        (-screened_rates(k_p_co55_to_ni56)*Y(Co55)*Y(H1)*state.rho + screened_rates(k_ni56_to_p_co55_derived)*Y(Ni56)) +
        (-screened_rates(k_he4_co55_to_cu59)*Y(Co55)*Y(He4)*state.rho + screened_rates(k_cu59_to_he4_co55_derived)*Y(Cu59)) +
        (screened_rates(k_he4_fe52_to_p_co55)*Y(Fe52)*Y(He4)*state.rho + -screened_rates(k_p_co55_to_he4_fe52_derived)*Y(Co55)*Y(H1)*state.rho) +
        (screened_rates(k_p_fe54_to_co55)*Y(Fe54)*Y(H1)*state.rho + -screened_rates(k_co55_to_p_fe54_derived)*Y(Co55)) +
        (-screened_rates(k_n_co55_to_co56)*Y(Co55)*Y(N)*state.rho + screened_rates(k_co56_to_n_co55_derived)*Y(Co56)) +
        (-screened_rates(k_n_co55_to_p_fe55)*Y(Co55)*Y(N)*state.rho + screened_rates(k_p_fe55_to_n_co55_derived)*Y(Fe55)*Y(H1)*state.rho) +
        (-screened_rates(k_he4_co55_to_p_ni58)*Y(Co55)*Y(He4)*state.rho + screened_rates(k_p_ni58_to_he4_co55_derived)*Y(Ni58)*Y(H1)*state.rho) +
        (-screened_rates(k_co55_to_fe55)*Y(Co55) + screened_rates(k_fe55_to_co55)*Y(Fe55));

    ydot_nuc(Co56) =
        (screened_rates(k_p_fe55_to_co56)*Y(Fe55)*Y(H1)*state.rho + -screened_rates(k_co56_to_p_fe55_derived)*Y(Co56)) +
        (screened_rates(k_n_co55_to_co56)*Y(Co55)*Y(N)*state.rho + -screened_rates(k_co56_to_n_co55_derived)*Y(Co56)) +
        (-screened_rates(k_n_co56_to_co57)*Y(Co56)*Y(N)*state.rho + screened_rates(k_co57_to_n_co56_derived)*Y(Co57)) +
        (-screened_rates(k_p_co56_to_ni57)*Y(Co56)*Y(H1)*state.rho + screened_rates(k_ni57_to_p_co56_derived)*Y(Ni57)) +
        (screened_rates(k_he4_fe53_to_p_co56)*Y(Fe53)*Y(He4)*state.rho + -screened_rates(k_p_co56_to_he4_fe53_derived)*Y(Co56)*Y(H1)*state.rho) +
        (-screened_rates(k_n_co56_to_p_fe56)*Y(Co56)*Y(N)*state.rho + screened_rates(k_p_fe56_to_n_co56_derived)*Y(Fe56)*Y(H1)*state.rho) +
        (screened_rates(k_n_ni56_to_p_co56)*Y(N)*Y(Ni56)*state.rho + -screened_rates(k_p_co56_to_n_ni56_derived)*Y(Co56)*Y(H1)*state.rho) +
        (screened_rates(k_n_cu59_to_he4_co56)*Y(Cu59)*Y(N)*state.rho + -screened_rates(k_he4_co56_to_n_cu59_derived)*Y(Co56)*Y(He4)*state.rho) +
        (-screened_rates(k_co56_to_fe56)*Y(Co56) + screened_rates(k_fe56_to_co56)*Y(Fe56)) +
        (screened_rates(k_ni56_to_co56)*Y(Ni56) + -screened_rates(k_co56_to_ni56)*Y(Co56));

    ydot_nuc(Co57) =
        (screened_rates(k_p_fe56_to_co57)*Y(Fe56)*Y(H1)*state.rho + -screened_rates(k_co57_to_p_fe56_derived)*Y(Co57)) +
        (screened_rates(k_n_co56_to_co57)*Y(Co56)*Y(N)*state.rho + -screened_rates(k_co57_to_n_co56_derived)*Y(Co57)) +
        (-screened_rates(k_p_co57_to_ni58)*Y(Co57)*Y(H1)*state.rho + screened_rates(k_ni58_to_p_co57_derived)*Y(Ni58)) +
        (-screened_rates(k_p_co57_to_he4_fe54)*Y(Co57)*Y(H1)*state.rho + screened_rates(k_he4_fe54_to_p_co57_derived)*Y(Fe54)*Y(He4)*state.rho) +
        (screened_rates(k_n_ni57_to_p_co57)*Y(N)*Y(Ni57)*state.rho + -screened_rates(k_p_co57_to_n_ni57_derived)*Y(Co57)*Y(H1)*state.rho) +
        (screened_rates(k_ni57_to_co57)*Y(Ni57) + -screened_rates(k_co57_to_ni57)*Y(Co57));

    ydot_nuc(Ni56) =
        (screened_rates(k_he4_fe52_to_ni56)*Y(Fe52)*Y(He4)*state.rho + -screened_rates(k_ni56_to_he4_fe52_derived)*Y(Ni56)) +
        (screened_rates(k_p_co55_to_ni56)*Y(Co55)*Y(H1)*state.rho + -screened_rates(k_ni56_to_p_co55_derived)*Y(Ni56)) +
        (-screened_rates(k_he4_ni56_to_zn60)*Y(He4)*Y(Ni56)*state.rho + screened_rates(k_zn60_to_he4_ni56_derived)*Y(Zn60)) +
        (screened_rates(k_p_cu59_to_he4_ni56)*Y(Cu59)*Y(H1)*state.rho + -screened_rates(k_he4_ni56_to_p_cu59_derived)*Y(He4)*Y(Ni56)*state.rho) +
        (-screened_rates(k_n_ni56_to_ni57)*Y(N)*Y(Ni56)*state.rho + screened_rates(k_ni57_to_n_ni56_derived)*Y(Ni57)) +
        (-screened_rates(k_n_ni56_to_p_co56)*Y(N)*Y(Ni56)*state.rho + screened_rates(k_p_co56_to_n_ni56_derived)*Y(Co56)*Y(H1)*state.rho) +
        (-screened_rates(k_n_ni56_to_he4_fe53)*Y(N)*Y(Ni56)*state.rho + screened_rates(k_he4_fe53_to_n_ni56_derived)*Y(Fe53)*Y(He4)*state.rho) +
        (-screened_rates(k_ni56_to_co56)*Y(Ni56) + screened_rates(k_co56_to_ni56)*Y(Co56));

    ydot_nuc(Ni57) =
        (screened_rates(k_he4_fe53_to_ni57)*Y(Fe53)*Y(He4)*state.rho + -screened_rates(k_ni57_to_he4_fe53_derived)*Y(Ni57)) +
        (screened_rates(k_p_co56_to_ni57)*Y(Co56)*Y(H1)*state.rho + -screened_rates(k_ni57_to_p_co56_derived)*Y(Ni57)) +
        (screened_rates(k_n_ni56_to_ni57)*Y(N)*Y(Ni56)*state.rho + -screened_rates(k_ni57_to_n_ni56_derived)*Y(Ni57)) +
        (-screened_rates(k_n_ni57_to_ni58)*Y(N)*Y(Ni57)*state.rho + screened_rates(k_ni58_to_n_ni57_derived)*Y(Ni58)) +
        (-screened_rates(k_n_ni57_to_p_co57)*Y(N)*Y(Ni57)*state.rho + screened_rates(k_p_co57_to_n_ni57_derived)*Y(Co57)*Y(H1)*state.rho) +
        (-screened_rates(k_n_ni57_to_he4_fe54)*Y(N)*Y(Ni57)*state.rho + screened_rates(k_he4_fe54_to_n_ni57_derived)*Y(Fe54)*Y(He4)*state.rho) +
        (screened_rates(k_n_zn60_to_he4_ni57)*Y(N)*Y(Zn60)*state.rho + -screened_rates(k_he4_ni57_to_n_zn60_derived)*Y(He4)*Y(Ni57)*state.rho) +
        (-screened_rates(k_ni57_to_co57)*Y(Ni57) + screened_rates(k_co57_to_ni57)*Y(Co57));

    ydot_nuc(Ni58) =
        (screened_rates(k_he4_fe54_to_ni58)*Y(Fe54)*Y(He4)*state.rho + -screened_rates(k_ni58_to_he4_fe54_derived)*Y(Ni58)) +
        (screened_rates(k_p_co57_to_ni58)*Y(Co57)*Y(H1)*state.rho + -screened_rates(k_ni58_to_p_co57_derived)*Y(Ni58)) +
        (screened_rates(k_n_ni57_to_ni58)*Y(N)*Y(Ni57)*state.rho + -screened_rates(k_ni58_to_n_ni57_derived)*Y(Ni58)) +
        (-screened_rates(k_p_ni58_to_cu59)*Y(Ni58)*Y(H1)*state.rho + screened_rates(k_cu59_to_p_ni58_derived)*Y(Cu59)) +
        (screened_rates(k_he4_co55_to_p_ni58)*Y(Co55)*Y(He4)*state.rho + -screened_rates(k_p_ni58_to_he4_co55_derived)*Y(Ni58)*Y(H1)*state.rho) +
        (-screened_rates(k_n_ni58_to_he4_fe55)*Y(N)*Y(Ni58)*state.rho + screened_rates(k_he4_fe55_to_n_ni58_derived)*Y(Fe55)*Y(He4)*state.rho);

    ydot_nuc(Cu59) =
        (screened_rates(k_he4_co55_to_cu59)*Y(Co55)*Y(He4)*state.rho + -screened_rates(k_cu59_to_he4_co55_derived)*Y(Cu59)) +
        (-screened_rates(k_p_cu59_to_zn60)*Y(Cu59)*Y(H1)*state.rho + screened_rates(k_zn60_to_p_cu59_derived)*Y(Zn60)) +
        (-screened_rates(k_p_cu59_to_he4_ni56)*Y(Cu59)*Y(H1)*state.rho + screened_rates(k_he4_ni56_to_p_cu59_derived)*Y(He4)*Y(Ni56)*state.rho) +
        (screened_rates(k_p_ni58_to_cu59)*Y(Ni58)*Y(H1)*state.rho + -screened_rates(k_cu59_to_p_ni58_derived)*Y(Cu59)) +
        (-screened_rates(k_n_cu59_to_he4_co56)*Y(Cu59)*Y(N)*state.rho + screened_rates(k_he4_co56_to_n_cu59_derived)*Y(Co56)*Y(He4)*state.rho);

    ydot_nuc(Zn60) =
        (screened_rates(k_he4_ni56_to_zn60)*Y(He4)*Y(Ni56)*state.rho + -screened_rates(k_zn60_to_he4_ni56_derived)*Y(Zn60)) +
        (screened_rates(k_p_cu59_to_zn60)*Y(Cu59)*Y(H1)*state.rho + -screened_rates(k_zn60_to_p_cu59_derived)*Y(Zn60)) +
        (-screened_rates(k_n_zn60_to_he4_ni57)*Y(N)*Y(Zn60)*state.rho + screened_rates(k_he4_ni57_to_n_zn60_derived)*Y(He4)*Y(Ni57)*state.rho);

}


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_rhs (burn_t& state, Array1D<Real, 1, neqs>& ydot)
{
    for (int i = 1; i <= neqs; ++i) {
        ydot(i) = 0.0_rt;
    }


    // Set molar abundances
    Array1D<Real, 1, NumSpec> Y;
    for (int i = 1; i <= NumSpec; ++i) {
        Y(i) = state.xn[i-1] * aion_inv[i-1];
    }

    // build the rates

    rate_t rate_eval;

    constexpr int do_T_derivatives = 0;
    evaluate_rates<do_T_derivatives, rate_t>(state, rate_eval);

    rhs_nuc(state, ydot, Y, rate_eval.screened_rates);

    // ion binding energy contributions

    Real enuc;
    ener_gener_rate(ydot, enuc);

    // include reaction neutrino losses (non-thermal) and gamma heating
    enuc += C::Legacy::n_A * Y(Co55) * rate_eval.add_energy_rate(k_co55_to_fe55);
    enuc += C::Legacy::n_A * Y(Co56) * rate_eval.add_energy_rate(k_co56_to_fe56);
    enuc += C::Legacy::n_A * Y(Co56) * rate_eval.add_energy_rate(k_co56_to_ni56);
    enuc += C::Legacy::n_A * Y(Co57) * rate_eval.add_energy_rate(k_co57_to_ni57);
    enuc += C::Legacy::n_A * Y(Fe55) * rate_eval.add_energy_rate(k_fe55_to_co55);
    enuc += C::Legacy::n_A * Y(Fe55) * rate_eval.add_energy_rate(k_fe55_to_mn55);
    enuc += C::Legacy::n_A * Y(Fe56) * rate_eval.add_energy_rate(k_fe56_to_co56);
    enuc += C::Legacy::n_A * Y(Mn55) * rate_eval.add_energy_rate(k_mn55_to_fe55);
    enuc += C::Legacy::n_A * Y(Ni56) * rate_eval.add_energy_rate(k_ni56_to_co56);
    enuc += C::Legacy::n_A * Y(Ni57) * rate_eval.add_energy_rate(k_ni57_to_co57);

    // Get the thermal neutrino losses

    Real sneut, dsneutdt, dsneutdd, dsnuda, dsnudz;
    constexpr int do_derivatives{0};
    sneut5<do_derivatives>(state.T, state.rho, state.abar, state.zbar, sneut, dsneutdt, dsneutdd, dsnuda, dsnudz);

    // Append the energy equation (this is erg/g/s)

    ydot(net_ienuc) = enuc - sneut;

}


template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void jac_nuc(const burn_t& state,
             MatrixType& jac,
             const Array1D<Real, 1, NumSpec>& Y,
             const Array1D<Real, 1, NumRates>& screened_rates)
{

    Real scratch;

    scratch = -screened_rates(k_n_co55_to_co56)*Y(Co55)*state.rho - screened_rates(k_n_co55_to_p_fe55)*Y(Co55)*state.rho - screened_rates(k_n_co56_to_co57)*Y(Co56)*state.rho - screened_rates(k_n_co56_to_p_fe56)*Y(Co56)*state.rho - screened_rates(k_n_cu59_to_he4_co56)*Y(Cu59)*state.rho - screened_rates(k_n_fe52_to_fe53)*Y(Fe52)*state.rho - screened_rates(k_n_fe53_to_fe54)*Y(Fe53)*state.rho - screened_rates(k_n_fe54_to_fe55)*Y(Fe54)*state.rho - screened_rates(k_n_fe55_to_fe56)*Y(Fe55)*state.rho - screened_rates(k_n_fe55_to_p_mn55)*Y(Fe55)*state.rho - screened_rates(k_n_ni56_to_he4_fe53)*Y(Ni56)*state.rho - screened_rates(k_n_ni56_to_ni57)*Y(Ni56)*state.rho - screened_rates(k_n_ni56_to_p_co56)*Y(Ni56)*state.rho - screened_rates(k_n_ni57_to_he4_fe54)*Y(Ni57)*state.rho - screened_rates(k_n_ni57_to_ni58)*Y(Ni57)*state.rho - screened_rates(k_n_ni57_to_p_co57)*Y(Ni57)*state.rho - screened_rates(k_n_ni58_to_he4_fe55)*Y(Ni58)*state.rho - screened_rates(k_n_to_p_weak_wc12) - screened_rates(k_n_zn60_to_he4_ni57)*Y(Zn60)*state.rho;
    jac.set(N, N, scratch);

    scratch = screened_rates(k_p_co56_to_n_ni56_derived)*Y(Co56)*state.rho + screened_rates(k_p_co57_to_n_ni57_derived)*Y(Co57)*state.rho + screened_rates(k_p_fe55_to_n_co55_derived)*Y(Fe55)*state.rho + screened_rates(k_p_fe56_to_n_co56_derived)*Y(Fe56)*state.rho + screened_rates(k_p_mn55_to_n_fe55_derived)*Y(Mn55)*state.rho;
    jac.set(N, H1, scratch);

    scratch = screened_rates(k_he4_co56_to_n_cu59_derived)*Y(Co56)*state.rho + screened_rates(k_he4_fe53_to_n_ni56_derived)*Y(Fe53)*state.rho + screened_rates(k_he4_fe54_to_n_ni57_derived)*Y(Fe54)*state.rho + screened_rates(k_he4_fe55_to_n_ni58_derived)*Y(Fe55)*state.rho + screened_rates(k_he4_ni57_to_n_zn60_derived)*Y(Ni57)*state.rho;
    jac.set(N, He4, scratch);

    scratch = screened_rates(k_p_mn55_to_n_fe55_derived)*Y(H1)*state.rho;
    jac.set(N, Mn55, scratch);

    scratch = -screened_rates(k_n_fe52_to_fe53)*Y(N)*state.rho;
    jac.set(N, Fe52, scratch);

    scratch = screened_rates(k_fe53_to_n_fe52_derived) + screened_rates(k_he4_fe53_to_n_ni56_derived)*Y(He4)*state.rho - screened_rates(k_n_fe53_to_fe54)*Y(N)*state.rho;
    jac.set(N, Fe53, scratch);

    scratch = screened_rates(k_fe54_to_n_fe53_derived) + screened_rates(k_he4_fe54_to_n_ni57_derived)*Y(He4)*state.rho - screened_rates(k_n_fe54_to_fe55)*Y(N)*state.rho;
    jac.set(N, Fe54, scratch);

    scratch = screened_rates(k_fe55_to_n_fe54_derived) + screened_rates(k_he4_fe55_to_n_ni58_derived)*Y(He4)*state.rho - screened_rates(k_n_fe55_to_fe56)*Y(N)*state.rho - screened_rates(k_n_fe55_to_p_mn55)*Y(N)*state.rho + screened_rates(k_p_fe55_to_n_co55_derived)*Y(H1)*state.rho;
    jac.set(N, Fe55, scratch);

    scratch = screened_rates(k_fe56_to_n_fe55_derived) + screened_rates(k_p_fe56_to_n_co56_derived)*Y(H1)*state.rho;
    jac.set(N, Fe56, scratch);

    scratch = -screened_rates(k_n_co55_to_co56)*Y(N)*state.rho - screened_rates(k_n_co55_to_p_fe55)*Y(N)*state.rho;
    jac.set(N, Co55, scratch);

    scratch = screened_rates(k_co56_to_n_co55_derived) + screened_rates(k_he4_co56_to_n_cu59_derived)*Y(He4)*state.rho - screened_rates(k_n_co56_to_co57)*Y(N)*state.rho - screened_rates(k_n_co56_to_p_fe56)*Y(N)*state.rho + screened_rates(k_p_co56_to_n_ni56_derived)*Y(H1)*state.rho;
    jac.set(N, Co56, scratch);

    scratch = screened_rates(k_co57_to_n_co56_derived) + screened_rates(k_p_co57_to_n_ni57_derived)*Y(H1)*state.rho;
    jac.set(N, Co57, scratch);

    scratch = -screened_rates(k_n_ni56_to_he4_fe53)*Y(N)*state.rho - screened_rates(k_n_ni56_to_ni57)*Y(N)*state.rho - screened_rates(k_n_ni56_to_p_co56)*Y(N)*state.rho;
    jac.set(N, Ni56, scratch);

    scratch = screened_rates(k_he4_ni57_to_n_zn60_derived)*Y(He4)*state.rho - screened_rates(k_n_ni57_to_he4_fe54)*Y(N)*state.rho - screened_rates(k_n_ni57_to_ni58)*Y(N)*state.rho - screened_rates(k_n_ni57_to_p_co57)*Y(N)*state.rho + screened_rates(k_ni57_to_n_ni56_derived);
    jac.set(N, Ni57, scratch);

    scratch = -screened_rates(k_n_ni58_to_he4_fe55)*Y(N)*state.rho + screened_rates(k_ni58_to_n_ni57_derived);
    jac.set(N, Ni58, scratch);

    scratch = -screened_rates(k_n_cu59_to_he4_co56)*Y(N)*state.rho;
    jac.set(N, Cu59, scratch);

    scratch = -screened_rates(k_n_zn60_to_he4_ni57)*Y(N)*state.rho;
    jac.set(N, Zn60, scratch);

    scratch = screened_rates(k_n_co55_to_p_fe55)*Y(Co55)*state.rho + screened_rates(k_n_co56_to_p_fe56)*Y(Co56)*state.rho + screened_rates(k_n_fe55_to_p_mn55)*Y(Fe55)*state.rho + screened_rates(k_n_ni56_to_p_co56)*Y(Ni56)*state.rho + screened_rates(k_n_ni57_to_p_co57)*Y(Ni57)*state.rho + screened_rates(k_n_to_p_weak_wc12);
    jac.set(H1, N, scratch);

    scratch = -screened_rates(k_p_al27_to_he4_mg24)*Y(Al27)*state.rho - screened_rates(k_p_al27_to_si28)*Y(Al27)*state.rho - screened_rates(k_p_c12_to_n13)*Y(C12)*state.rho - screened_rates(k_p_co55_to_he4_fe52_derived)*Y(Co55)*state.rho - screened_rates(k_p_co55_to_ni56)*Y(Co55)*state.rho - screened_rates(k_p_co56_to_he4_fe53_derived)*Y(Co56)*state.rho - screened_rates(k_p_co56_to_n_ni56_derived)*Y(Co56)*state.rho - screened_rates(k_p_co56_to_ni57)*Y(Co56)*state.rho - screened_rates(k_p_co57_to_he4_fe54)*Y(Co57)*state.rho - screened_rates(k_p_co57_to_n_ni57_derived)*Y(Co57)*state.rho - screened_rates(k_p_co57_to_ni58)*Y(Co57)*state.rho - screened_rates(k_p_cu59_to_he4_ni56)*Y(Cu59)*state.rho - screened_rates(k_p_cu59_to_zn60)*Y(Cu59)*state.rho - screened_rates(k_p_fe54_to_co55)*Y(Fe54)*state.rho - screened_rates(k_p_fe54_to_he4_mn51_derived)*Y(Fe54)*state.rho - screened_rates(k_p_fe55_to_co56)*Y(Fe55)*state.rho - screened_rates(k_p_fe55_to_n_co55_derived)*Y(Fe55)*state.rho - screened_rates(k_p_fe56_to_co57)*Y(Fe56)*state.rho - screened_rates(k_p_fe56_to_n_co56_derived)*Y(Fe56)*state.rho - screened_rates(k_p_mn51_to_fe52)*Y(Mn51)*state.rho - screened_rates(k_p_mn51_to_he4_cr48_derived)*Y(Mn51)*state.rho - screened_rates(k_p_mn55_to_fe56)*Y(Mn55)*state.rho - screened_rates(k_p_mn55_to_n_fe55_derived)*Y(Mn55)*state.rho - screened_rates(k_p_na23_to_he4_ne20)*Y(Na23)*state.rho - screened_rates(k_p_na23_to_mg24)*Y(Na23)*state.rho - screened_rates(k_p_ne21_to_he4_f18_derived)*Y(Ne21)*state.rho - screened_rates(k_p_ne21_to_na22)*Y(Ne21)*state.rho - screened_rates(k_p_ni58_to_cu59)*Y(Ni58)*state.rho - screened_rates(k_p_ni58_to_he4_co55_derived)*Y(Ni58)*state.rho - screened_rates(k_p_o16_to_he4_n13_derived)*Y(O16)*state.rho - screened_rates(k_p_p31_to_he4_si28)*Y(P31)*state.rho - screened_rates(k_p_p31_to_s32)*Y(P31)*state.rho;
    jac.set(H1, H1, scratch);

    scratch = screened_rates(k_he4_co55_to_p_ni58)*Y(Co55)*state.rho + screened_rates(k_he4_cr48_to_p_mn51)*Y(Cr48)*state.rho + screened_rates(k_he4_f18_to_p_ne21)*Y(F18)*state.rho + screened_rates(k_he4_fe52_to_p_co55)*Y(Fe52)*state.rho + screened_rates(k_he4_fe53_to_p_co56)*Y(Fe53)*state.rho + screened_rates(k_he4_fe54_to_p_co57_derived)*Y(Fe54)*state.rho + screened_rates(k_he4_mg24_to_p_al27_derived)*Y(Mg24)*state.rho + screened_rates(k_he4_mn51_to_p_fe54)*Y(Mn51)*state.rho + screened_rates(k_he4_n13_to_p_o16)*Y(N13)*state.rho + screened_rates(k_he4_ne20_to_p_na23_derived)*Y(Ne20)*state.rho + screened_rates(k_he4_ni56_to_p_cu59_derived)*Y(Ni56)*state.rho + screened_rates(k_he4_si28_to_p_p31_derived)*Y(Si28)*state.rho;
    jac.set(H1, He4, scratch);

    scratch = 1.0*screened_rates(k_c12_c12_to_p_na23)*Y(C12)*state.rho + screened_rates(k_c12_o16_to_p_al27)*Y(O16)*state.rho - screened_rates(k_p_c12_to_n13)*Y(H1)*state.rho;
    jac.set(H1, C12, scratch);

    scratch = screened_rates(k_he4_n13_to_p_o16)*Y(He4)*state.rho + screened_rates(k_n13_to_p_c12_derived);
    jac.set(H1, N13, scratch);

    scratch = screened_rates(k_c12_o16_to_p_al27)*Y(C12)*state.rho + 1.0*screened_rates(k_o16_o16_to_p_p31)*Y(O16)*state.rho - screened_rates(k_p_o16_to_he4_n13_derived)*Y(H1)*state.rho;
    jac.set(H1, O16, scratch);

    scratch = screened_rates(k_he4_f18_to_p_ne21)*Y(He4)*state.rho;
    jac.set(H1, F18, scratch);

    scratch = screened_rates(k_he4_ne20_to_p_na23_derived)*Y(He4)*state.rho;
    jac.set(H1, Ne20, scratch);

    scratch = -screened_rates(k_p_ne21_to_he4_f18_derived)*Y(H1)*state.rho - screened_rates(k_p_ne21_to_na22)*Y(H1)*state.rho;
    jac.set(H1, Ne21, scratch);

    scratch = screened_rates(k_na22_to_p_ne21_derived);
    jac.set(H1, Na22, scratch);

    scratch = -screened_rates(k_p_na23_to_he4_ne20)*Y(H1)*state.rho - screened_rates(k_p_na23_to_mg24)*Y(H1)*state.rho;
    jac.set(H1, Na23, scratch);

    scratch = screened_rates(k_he4_mg24_to_p_al27_derived)*Y(He4)*state.rho + screened_rates(k_mg24_to_p_na23_derived);
    jac.set(H1, Mg24, scratch);

    scratch = -screened_rates(k_p_al27_to_he4_mg24)*Y(H1)*state.rho - screened_rates(k_p_al27_to_si28)*Y(H1)*state.rho;
    jac.set(H1, Al27, scratch);

    scratch = screened_rates(k_he4_si28_to_p_p31_derived)*Y(He4)*state.rho + screened_rates(k_si28_to_p_al27_derived);
    jac.set(H1, Si28, scratch);

    scratch = -screened_rates(k_p_p31_to_he4_si28)*Y(H1)*state.rho - screened_rates(k_p_p31_to_s32)*Y(H1)*state.rho;
    jac.set(H1, P31, scratch);

    scratch = screened_rates(k_s32_to_p_p31_derived);
    jac.set(H1, S32, scratch);

    scratch = screened_rates(k_he4_cr48_to_p_mn51)*Y(He4)*state.rho;
    jac.set(H1, Cr48, scratch);

    scratch = screened_rates(k_he4_mn51_to_p_fe54)*Y(He4)*state.rho - screened_rates(k_p_mn51_to_fe52)*Y(H1)*state.rho - screened_rates(k_p_mn51_to_he4_cr48_derived)*Y(H1)*state.rho;
    jac.set(H1, Mn51, scratch);

    scratch = -screened_rates(k_p_mn55_to_fe56)*Y(H1)*state.rho - screened_rates(k_p_mn55_to_n_fe55_derived)*Y(H1)*state.rho;
    jac.set(H1, Mn55, scratch);

    scratch = screened_rates(k_fe52_to_p_mn51_derived) + screened_rates(k_he4_fe52_to_p_co55)*Y(He4)*state.rho;
    jac.set(H1, Fe52, scratch);

    scratch = screened_rates(k_he4_fe53_to_p_co56)*Y(He4)*state.rho;
    jac.set(H1, Fe53, scratch);

    scratch = screened_rates(k_he4_fe54_to_p_co57_derived)*Y(He4)*state.rho - screened_rates(k_p_fe54_to_co55)*Y(H1)*state.rho - screened_rates(k_p_fe54_to_he4_mn51_derived)*Y(H1)*state.rho;
    jac.set(H1, Fe54, scratch);

    scratch = screened_rates(k_n_fe55_to_p_mn55)*Y(N)*state.rho - screened_rates(k_p_fe55_to_co56)*Y(H1)*state.rho - screened_rates(k_p_fe55_to_n_co55_derived)*Y(H1)*state.rho;
    jac.set(H1, Fe55, scratch);

    scratch = screened_rates(k_fe56_to_p_mn55_derived) - screened_rates(k_p_fe56_to_co57)*Y(H1)*state.rho - screened_rates(k_p_fe56_to_n_co56_derived)*Y(H1)*state.rho;
    jac.set(H1, Fe56, scratch);

    scratch = screened_rates(k_co55_to_p_fe54_derived) + screened_rates(k_he4_co55_to_p_ni58)*Y(He4)*state.rho + screened_rates(k_n_co55_to_p_fe55)*Y(N)*state.rho - screened_rates(k_p_co55_to_he4_fe52_derived)*Y(H1)*state.rho - screened_rates(k_p_co55_to_ni56)*Y(H1)*state.rho;
    jac.set(H1, Co55, scratch);

    scratch = screened_rates(k_co56_to_p_fe55_derived) + screened_rates(k_n_co56_to_p_fe56)*Y(N)*state.rho - screened_rates(k_p_co56_to_he4_fe53_derived)*Y(H1)*state.rho - screened_rates(k_p_co56_to_n_ni56_derived)*Y(H1)*state.rho - screened_rates(k_p_co56_to_ni57)*Y(H1)*state.rho;
    jac.set(H1, Co56, scratch);

    scratch = screened_rates(k_co57_to_p_fe56_derived) - screened_rates(k_p_co57_to_he4_fe54)*Y(H1)*state.rho - screened_rates(k_p_co57_to_n_ni57_derived)*Y(H1)*state.rho - screened_rates(k_p_co57_to_ni58)*Y(H1)*state.rho;
    jac.set(H1, Co57, scratch);

    scratch = screened_rates(k_he4_ni56_to_p_cu59_derived)*Y(He4)*state.rho + screened_rates(k_n_ni56_to_p_co56)*Y(N)*state.rho + screened_rates(k_ni56_to_p_co55_derived);
    jac.set(H1, Ni56, scratch);

    scratch = screened_rates(k_n_ni57_to_p_co57)*Y(N)*state.rho + screened_rates(k_ni57_to_p_co56_derived);
    jac.set(H1, Ni57, scratch);

    scratch = screened_rates(k_ni58_to_p_co57_derived) - screened_rates(k_p_ni58_to_cu59)*Y(H1)*state.rho - screened_rates(k_p_ni58_to_he4_co55_derived)*Y(H1)*state.rho;
    jac.set(H1, Ni58, scratch);

    scratch = screened_rates(k_cu59_to_p_ni58_derived) - screened_rates(k_p_cu59_to_he4_ni56)*Y(H1)*state.rho - screened_rates(k_p_cu59_to_zn60)*Y(H1)*state.rho;
    jac.set(H1, Cu59, scratch);

    scratch = screened_rates(k_zn60_to_p_cu59_derived);
    jac.set(H1, Zn60, scratch);

    scratch = screened_rates(k_n_cu59_to_he4_co56)*Y(Cu59)*state.rho + screened_rates(k_n_ni56_to_he4_fe53)*Y(Ni56)*state.rho + screened_rates(k_n_ni57_to_he4_fe54)*Y(Ni57)*state.rho + screened_rates(k_n_ni58_to_he4_fe55)*Y(Ni58)*state.rho + screened_rates(k_n_zn60_to_he4_ni57)*Y(Zn60)*state.rho;
    jac.set(He4, N, scratch);

    scratch = screened_rates(k_p_al27_to_he4_mg24)*Y(Al27)*state.rho + screened_rates(k_p_co55_to_he4_fe52_derived)*Y(Co55)*state.rho + screened_rates(k_p_co56_to_he4_fe53_derived)*Y(Co56)*state.rho + screened_rates(k_p_co57_to_he4_fe54)*Y(Co57)*state.rho + screened_rates(k_p_cu59_to_he4_ni56)*Y(Cu59)*state.rho + screened_rates(k_p_fe54_to_he4_mn51_derived)*Y(Fe54)*state.rho + screened_rates(k_p_mn51_to_he4_cr48_derived)*Y(Mn51)*state.rho + screened_rates(k_p_na23_to_he4_ne20)*Y(Na23)*state.rho + screened_rates(k_p_ne21_to_he4_f18_derived)*Y(Ne21)*state.rho + screened_rates(k_p_ni58_to_he4_co55_derived)*Y(Ni58)*state.rho + screened_rates(k_p_o16_to_he4_n13_derived)*Y(O16)*state.rho + screened_rates(k_p_p31_to_he4_si28)*Y(P31)*state.rho;
    jac.set(He4, H1, scratch);

    scratch = -screened_rates(k_ar36_he4_to_ca40_approx)*Y(Ar36)*state.rho - screened_rates(k_ca40_he4_to_ti44_approx)*Y(Ca40)*state.rho - screened_rates(k_he4_c12_to_o16)*Y(C12)*state.rho - screened_rates(k_he4_co55_to_cu59)*Y(Co55)*state.rho - screened_rates(k_he4_co55_to_p_ni58)*Y(Co55)*state.rho - screened_rates(k_he4_co56_to_n_cu59_derived)*Y(Co56)*state.rho - screened_rates(k_he4_cr48_to_fe52)*Y(Cr48)*state.rho - screened_rates(k_he4_cr48_to_p_mn51)*Y(Cr48)*state.rho - screened_rates(k_he4_f18_to_na22)*Y(F18)*state.rho - screened_rates(k_he4_f18_to_p_ne21)*Y(F18)*state.rho - screened_rates(k_he4_fe52_to_ni56)*Y(Fe52)*state.rho - screened_rates(k_he4_fe52_to_p_co55)*Y(Fe52)*state.rho - screened_rates(k_he4_fe53_to_n_ni56_derived)*Y(Fe53)*state.rho - screened_rates(k_he4_fe53_to_ni57)*Y(Fe53)*state.rho - screened_rates(k_he4_fe53_to_p_co56)*Y(Fe53)*state.rho - screened_rates(k_he4_fe54_to_n_ni57_derived)*Y(Fe54)*state.rho - screened_rates(k_he4_fe54_to_ni58)*Y(Fe54)*state.rho - screened_rates(k_he4_fe54_to_p_co57_derived)*Y(Fe54)*state.rho - screened_rates(k_he4_fe55_to_n_ni58_derived)*Y(Fe55)*state.rho - 1.5*screened_rates(k_he4_he4_he4_to_c12)*std::pow(Y(He4), 2)*std::pow(state.rho, 2) - screened_rates(k_he4_mg24_to_p_al27_derived)*Y(Mg24)*state.rho - screened_rates(k_he4_mg24_to_si28)*Y(Mg24)*state.rho - screened_rates(k_he4_mn51_to_co55)*Y(Mn51)*state.rho - screened_rates(k_he4_mn51_to_p_fe54)*Y(Mn51)*state.rho - screened_rates(k_he4_n13_to_p_o16)*Y(N13)*state.rho - screened_rates(k_he4_n14_to_f18)*Y(N14)*state.rho - screened_rates(k_he4_ne20_to_mg24)*Y(Ne20)*state.rho - screened_rates(k_he4_ne20_to_p_na23_derived)*Y(Ne20)*state.rho - screened_rates(k_he4_ni56_to_p_cu59_derived)*Y(Ni56)*state.rho - screened_rates(k_he4_ni56_to_zn60)*Y(Ni56)*state.rho - screened_rates(k_he4_ni57_to_n_zn60_derived)*Y(Ni57)*state.rho - screened_rates(k_he4_o16_to_ne20)*Y(O16)*state.rho - screened_rates(k_he4_si28_to_p_p31_derived)*Y(Si28)*state.rho - screened_rates(k_he4_si28_to_s32)*Y(Si28)*state.rho - screened_rates(k_s32_he4_to_ar36_approx)*Y(S32)*state.rho - screened_rates(k_ti44_he4_to_cr48_approx)*Y(Ti44)*state.rho;
    jac.set(He4, He4, scratch);

    scratch = 1.0*screened_rates(k_c12_c12_to_he4_ne20)*Y(C12)*state.rho + screened_rates(k_c12_o16_to_he4_mg24)*Y(O16)*state.rho + 3.0*screened_rates(k_c12_to_he4_he4_he4_derived) - screened_rates(k_he4_c12_to_o16)*Y(He4)*state.rho;
    jac.set(He4, C12, scratch);

    scratch = -screened_rates(k_he4_n13_to_p_o16)*Y(He4)*state.rho;
    jac.set(He4, N13, scratch);

    scratch = -screened_rates(k_he4_n14_to_f18)*Y(He4)*state.rho;
    jac.set(He4, N14, scratch);

    scratch = screened_rates(k_c12_o16_to_he4_mg24)*Y(C12)*state.rho - screened_rates(k_he4_o16_to_ne20)*Y(He4)*state.rho + 1.0*screened_rates(k_o16_o16_to_he4_si28)*Y(O16)*state.rho + screened_rates(k_o16_to_he4_c12_derived) + screened_rates(k_p_o16_to_he4_n13_derived)*Y(H1)*state.rho;
    jac.set(He4, O16, scratch);

    scratch = screened_rates(k_f18_to_he4_n14_derived) - screened_rates(k_he4_f18_to_na22)*Y(He4)*state.rho - screened_rates(k_he4_f18_to_p_ne21)*Y(He4)*state.rho;
    jac.set(He4, F18, scratch);

    scratch = -screened_rates(k_he4_ne20_to_mg24)*Y(He4)*state.rho - screened_rates(k_he4_ne20_to_p_na23_derived)*Y(He4)*state.rho + screened_rates(k_ne20_to_he4_o16_derived);
    jac.set(He4, Ne20, scratch);

    scratch = screened_rates(k_p_ne21_to_he4_f18_derived)*Y(H1)*state.rho;
    jac.set(He4, Ne21, scratch);

    scratch = screened_rates(k_na22_to_he4_f18_derived);
    jac.set(He4, Na22, scratch);

    scratch = screened_rates(k_p_na23_to_he4_ne20)*Y(H1)*state.rho;
    jac.set(He4, Na23, scratch);

    scratch = -screened_rates(k_he4_mg24_to_p_al27_derived)*Y(He4)*state.rho - screened_rates(k_he4_mg24_to_si28)*Y(He4)*state.rho + screened_rates(k_mg24_to_he4_ne20_derived);
    jac.set(He4, Mg24, scratch);

    scratch = screened_rates(k_p_al27_to_he4_mg24)*Y(H1)*state.rho;
    jac.set(He4, Al27, scratch);

    scratch = -screened_rates(k_he4_si28_to_p_p31_derived)*Y(He4)*state.rho - screened_rates(k_he4_si28_to_s32)*Y(He4)*state.rho + screened_rates(k_si28_to_he4_mg24_derived);
    jac.set(He4, Si28, scratch);

    scratch = screened_rates(k_p_p31_to_he4_si28)*Y(H1)*state.rho;
    jac.set(He4, P31, scratch);

    scratch = -screened_rates(k_s32_he4_to_ar36_approx)*Y(He4)*state.rho + screened_rates(k_s32_to_he4_si28_derived);
    jac.set(He4, S32, scratch);

    scratch = -screened_rates(k_ar36_he4_to_ca40_approx)*Y(He4)*state.rho + screened_rates(k_ar36_to_s32_he4_approx);
    jac.set(He4, Ar36, scratch);

    scratch = -screened_rates(k_ca40_he4_to_ti44_approx)*Y(He4)*state.rho + screened_rates(k_ca40_to_ar36_he4_approx);
    jac.set(He4, Ca40, scratch);

    scratch = -screened_rates(k_ti44_he4_to_cr48_approx)*Y(He4)*state.rho + screened_rates(k_ti44_to_ca40_he4_approx);
    jac.set(He4, Ti44, scratch);

    scratch = screened_rates(k_cr48_to_ti44_he4_approx) - screened_rates(k_he4_cr48_to_fe52)*Y(He4)*state.rho - screened_rates(k_he4_cr48_to_p_mn51)*Y(He4)*state.rho;
    jac.set(He4, Cr48, scratch);

    scratch = -screened_rates(k_he4_mn51_to_co55)*Y(He4)*state.rho - screened_rates(k_he4_mn51_to_p_fe54)*Y(He4)*state.rho + screened_rates(k_p_mn51_to_he4_cr48_derived)*Y(H1)*state.rho;
    jac.set(He4, Mn51, scratch);

    scratch = screened_rates(k_fe52_to_he4_cr48_derived) - screened_rates(k_he4_fe52_to_ni56)*Y(He4)*state.rho - screened_rates(k_he4_fe52_to_p_co55)*Y(He4)*state.rho;
    jac.set(He4, Fe52, scratch);

    scratch = -screened_rates(k_he4_fe53_to_n_ni56_derived)*Y(He4)*state.rho - screened_rates(k_he4_fe53_to_ni57)*Y(He4)*state.rho - screened_rates(k_he4_fe53_to_p_co56)*Y(He4)*state.rho;
    jac.set(He4, Fe53, scratch);

    scratch = -screened_rates(k_he4_fe54_to_n_ni57_derived)*Y(He4)*state.rho - screened_rates(k_he4_fe54_to_ni58)*Y(He4)*state.rho - screened_rates(k_he4_fe54_to_p_co57_derived)*Y(He4)*state.rho + screened_rates(k_p_fe54_to_he4_mn51_derived)*Y(H1)*state.rho;
    jac.set(He4, Fe54, scratch);

    scratch = -screened_rates(k_he4_fe55_to_n_ni58_derived)*Y(He4)*state.rho;
    jac.set(He4, Fe55, scratch);

    scratch = screened_rates(k_co55_to_he4_mn51_derived) - screened_rates(k_he4_co55_to_cu59)*Y(He4)*state.rho - screened_rates(k_he4_co55_to_p_ni58)*Y(He4)*state.rho + screened_rates(k_p_co55_to_he4_fe52_derived)*Y(H1)*state.rho;
    jac.set(He4, Co55, scratch);

    scratch = -screened_rates(k_he4_co56_to_n_cu59_derived)*Y(He4)*state.rho + screened_rates(k_p_co56_to_he4_fe53_derived)*Y(H1)*state.rho;
    jac.set(He4, Co56, scratch);

    scratch = screened_rates(k_p_co57_to_he4_fe54)*Y(H1)*state.rho;
    jac.set(He4, Co57, scratch);

    scratch = -screened_rates(k_he4_ni56_to_p_cu59_derived)*Y(He4)*state.rho - screened_rates(k_he4_ni56_to_zn60)*Y(He4)*state.rho + screened_rates(k_n_ni56_to_he4_fe53)*Y(N)*state.rho + screened_rates(k_ni56_to_he4_fe52_derived);
    jac.set(He4, Ni56, scratch);

    scratch = -screened_rates(k_he4_ni57_to_n_zn60_derived)*Y(He4)*state.rho + screened_rates(k_n_ni57_to_he4_fe54)*Y(N)*state.rho + screened_rates(k_ni57_to_he4_fe53_derived);
    jac.set(He4, Ni57, scratch);

    scratch = screened_rates(k_n_ni58_to_he4_fe55)*Y(N)*state.rho + screened_rates(k_ni58_to_he4_fe54_derived) + screened_rates(k_p_ni58_to_he4_co55_derived)*Y(H1)*state.rho;
    jac.set(He4, Ni58, scratch);

    scratch = screened_rates(k_cu59_to_he4_co55_derived) + screened_rates(k_n_cu59_to_he4_co56)*Y(N)*state.rho + screened_rates(k_p_cu59_to_he4_ni56)*Y(H1)*state.rho;
    jac.set(He4, Cu59, scratch);

    scratch = screened_rates(k_n_zn60_to_he4_ni57)*Y(N)*state.rho + screened_rates(k_zn60_to_he4_ni56_derived);
    jac.set(He4, Zn60, scratch);

    scratch = -screened_rates(k_p_c12_to_n13)*Y(C12)*state.rho;
    jac.set(C12, H1, scratch);

    scratch = -screened_rates(k_he4_c12_to_o16)*Y(C12)*state.rho + 0.5*screened_rates(k_he4_he4_he4_to_c12)*std::pow(Y(He4), 2)*std::pow(state.rho, 2);
    jac.set(C12, He4, scratch);

    scratch = -2.0*screened_rates(k_c12_c12_to_he4_ne20)*Y(C12)*state.rho - 2.0*screened_rates(k_c12_c12_to_mg24_modified)*Y(C12)*state.rho - 2.0*screened_rates(k_c12_c12_to_p_na23)*Y(C12)*state.rho - screened_rates(k_c12_o16_to_he4_mg24)*Y(O16)*state.rho - screened_rates(k_c12_o16_to_p_al27)*Y(O16)*state.rho - screened_rates(k_c12_o16_to_si28_modified)*Y(O16)*state.rho - screened_rates(k_c12_to_he4_he4_he4_derived) - screened_rates(k_he4_c12_to_o16)*Y(He4)*state.rho - screened_rates(k_p_c12_to_n13)*Y(H1)*state.rho;
    jac.set(C12, C12, scratch);

    scratch = screened_rates(k_n13_to_p_c12_derived);
    jac.set(C12, N13, scratch);

    scratch = -screened_rates(k_c12_o16_to_he4_mg24)*Y(C12)*state.rho - screened_rates(k_c12_o16_to_p_al27)*Y(C12)*state.rho - screened_rates(k_c12_o16_to_si28_modified)*Y(C12)*state.rho + screened_rates(k_o16_to_he4_c12_derived);
    jac.set(C12, O16, scratch);

    scratch = screened_rates(k_p_c12_to_n13)*Y(C12)*state.rho + screened_rates(k_p_o16_to_he4_n13_derived)*Y(O16)*state.rho;
    jac.set(N13, H1, scratch);

    scratch = -screened_rates(k_he4_n13_to_p_o16)*Y(N13)*state.rho;
    jac.set(N13, He4, scratch);

    scratch = screened_rates(k_p_c12_to_n13)*Y(H1)*state.rho;
    jac.set(N13, C12, scratch);

    scratch = -screened_rates(k_he4_n13_to_p_o16)*Y(He4)*state.rho - screened_rates(k_n13_to_p_c12_derived);
    jac.set(N13, N13, scratch);

    scratch = screened_rates(k_p_o16_to_he4_n13_derived)*Y(H1)*state.rho;
    jac.set(N13, O16, scratch);

    scratch = -screened_rates(k_he4_n14_to_f18)*Y(N14)*state.rho;
    jac.set(N14, He4, scratch);

    scratch = -screened_rates(k_he4_n14_to_f18)*Y(He4)*state.rho;
    jac.set(N14, N14, scratch);

    scratch = screened_rates(k_f18_to_he4_n14_derived);
    jac.set(N14, F18, scratch);

    scratch = -screened_rates(k_p_o16_to_he4_n13_derived)*Y(O16)*state.rho;
    jac.set(O16, H1, scratch);

    scratch = screened_rates(k_he4_c12_to_o16)*Y(C12)*state.rho + screened_rates(k_he4_n13_to_p_o16)*Y(N13)*state.rho - screened_rates(k_he4_o16_to_ne20)*Y(O16)*state.rho;
    jac.set(O16, He4, scratch);

    scratch = -screened_rates(k_c12_o16_to_he4_mg24)*Y(O16)*state.rho - screened_rates(k_c12_o16_to_p_al27)*Y(O16)*state.rho - screened_rates(k_c12_o16_to_si28_modified)*Y(O16)*state.rho + screened_rates(k_he4_c12_to_o16)*Y(He4)*state.rho;
    jac.set(O16, C12, scratch);

    scratch = screened_rates(k_he4_n13_to_p_o16)*Y(He4)*state.rho;
    jac.set(O16, N13, scratch);

    scratch = -screened_rates(k_c12_o16_to_he4_mg24)*Y(C12)*state.rho - screened_rates(k_c12_o16_to_p_al27)*Y(C12)*state.rho - screened_rates(k_c12_o16_to_si28_modified)*Y(C12)*state.rho - screened_rates(k_he4_o16_to_ne20)*Y(He4)*state.rho - 2.0*screened_rates(k_o16_o16_to_he4_si28)*Y(O16)*state.rho - 2.0*screened_rates(k_o16_o16_to_p_p31)*Y(O16)*state.rho - 2.0*screened_rates(k_o16_o16_to_s32_modified)*Y(O16)*state.rho - screened_rates(k_o16_to_he4_c12_derived) - screened_rates(k_p_o16_to_he4_n13_derived)*Y(H1)*state.rho;
    jac.set(O16, O16, scratch);

    scratch = screened_rates(k_ne20_to_he4_o16_derived);
    jac.set(O16, Ne20, scratch);

    scratch = screened_rates(k_p_ne21_to_he4_f18_derived)*Y(Ne21)*state.rho;
    jac.set(F18, H1, scratch);

    scratch = -screened_rates(k_he4_f18_to_na22)*Y(F18)*state.rho - screened_rates(k_he4_f18_to_p_ne21)*Y(F18)*state.rho + screened_rates(k_he4_n14_to_f18)*Y(N14)*state.rho;
    jac.set(F18, He4, scratch);

    scratch = screened_rates(k_he4_n14_to_f18)*Y(He4)*state.rho;
    jac.set(F18, N14, scratch);

    scratch = -screened_rates(k_f18_to_he4_n14_derived) - screened_rates(k_he4_f18_to_na22)*Y(He4)*state.rho - screened_rates(k_he4_f18_to_p_ne21)*Y(He4)*state.rho;
    jac.set(F18, F18, scratch);

    scratch = screened_rates(k_p_ne21_to_he4_f18_derived)*Y(H1)*state.rho;
    jac.set(F18, Ne21, scratch);

    scratch = screened_rates(k_na22_to_he4_f18_derived);
    jac.set(F18, Na22, scratch);

    scratch = screened_rates(k_p_na23_to_he4_ne20)*Y(Na23)*state.rho;
    jac.set(Ne20, H1, scratch);

    scratch = -screened_rates(k_he4_ne20_to_mg24)*Y(Ne20)*state.rho - screened_rates(k_he4_ne20_to_p_na23_derived)*Y(Ne20)*state.rho + screened_rates(k_he4_o16_to_ne20)*Y(O16)*state.rho;
    jac.set(Ne20, He4, scratch);

    scratch = 1.0*screened_rates(k_c12_c12_to_he4_ne20)*Y(C12)*state.rho;
    jac.set(Ne20, C12, scratch);

    scratch = screened_rates(k_he4_o16_to_ne20)*Y(He4)*state.rho;
    jac.set(Ne20, O16, scratch);

    scratch = -screened_rates(k_he4_ne20_to_mg24)*Y(He4)*state.rho - screened_rates(k_he4_ne20_to_p_na23_derived)*Y(He4)*state.rho - screened_rates(k_ne20_to_he4_o16_derived);
    jac.set(Ne20, Ne20, scratch);

    scratch = screened_rates(k_p_na23_to_he4_ne20)*Y(H1)*state.rho;
    jac.set(Ne20, Na23, scratch);

    scratch = screened_rates(k_mg24_to_he4_ne20_derived);
    jac.set(Ne20, Mg24, scratch);

    scratch = -screened_rates(k_p_ne21_to_he4_f18_derived)*Y(Ne21)*state.rho - screened_rates(k_p_ne21_to_na22)*Y(Ne21)*state.rho;
    jac.set(Ne21, H1, scratch);

    scratch = screened_rates(k_he4_f18_to_p_ne21)*Y(F18)*state.rho;
    jac.set(Ne21, He4, scratch);

    scratch = screened_rates(k_he4_f18_to_p_ne21)*Y(He4)*state.rho;
    jac.set(Ne21, F18, scratch);

    scratch = -screened_rates(k_p_ne21_to_he4_f18_derived)*Y(H1)*state.rho - screened_rates(k_p_ne21_to_na22)*Y(H1)*state.rho;
    jac.set(Ne21, Ne21, scratch);

    scratch = screened_rates(k_na22_to_p_ne21_derived);
    jac.set(Ne21, Na22, scratch);

    scratch = screened_rates(k_p_ne21_to_na22)*Y(Ne21)*state.rho;
    jac.set(Na22, H1, scratch);

    scratch = screened_rates(k_he4_f18_to_na22)*Y(F18)*state.rho;
    jac.set(Na22, He4, scratch);

    scratch = screened_rates(k_he4_f18_to_na22)*Y(He4)*state.rho;
    jac.set(Na22, F18, scratch);

    scratch = screened_rates(k_p_ne21_to_na22)*Y(H1)*state.rho;
    jac.set(Na22, Ne21, scratch);

    scratch = -screened_rates(k_na22_to_he4_f18_derived) - screened_rates(k_na22_to_p_ne21_derived);
    jac.set(Na22, Na22, scratch);

    scratch = -screened_rates(k_p_na23_to_he4_ne20)*Y(Na23)*state.rho - screened_rates(k_p_na23_to_mg24)*Y(Na23)*state.rho;
    jac.set(Na23, H1, scratch);

    scratch = screened_rates(k_he4_ne20_to_p_na23_derived)*Y(Ne20)*state.rho;
    jac.set(Na23, He4, scratch);

    scratch = 1.0*screened_rates(k_c12_c12_to_p_na23)*Y(C12)*state.rho;
    jac.set(Na23, C12, scratch);

    scratch = screened_rates(k_he4_ne20_to_p_na23_derived)*Y(He4)*state.rho;
    jac.set(Na23, Ne20, scratch);

    scratch = -screened_rates(k_p_na23_to_he4_ne20)*Y(H1)*state.rho - screened_rates(k_p_na23_to_mg24)*Y(H1)*state.rho;
    jac.set(Na23, Na23, scratch);

    scratch = screened_rates(k_mg24_to_p_na23_derived);
    jac.set(Na23, Mg24, scratch);

    scratch = screened_rates(k_p_al27_to_he4_mg24)*Y(Al27)*state.rho + screened_rates(k_p_na23_to_mg24)*Y(Na23)*state.rho;
    jac.set(Mg24, H1, scratch);

    scratch = -screened_rates(k_he4_mg24_to_p_al27_derived)*Y(Mg24)*state.rho - screened_rates(k_he4_mg24_to_si28)*Y(Mg24)*state.rho + screened_rates(k_he4_ne20_to_mg24)*Y(Ne20)*state.rho;
    jac.set(Mg24, He4, scratch);

    scratch = 1.0*screened_rates(k_c12_c12_to_mg24_modified)*Y(C12)*state.rho + screened_rates(k_c12_o16_to_he4_mg24)*Y(O16)*state.rho;
    jac.set(Mg24, C12, scratch);

    scratch = screened_rates(k_c12_o16_to_he4_mg24)*Y(C12)*state.rho;
    jac.set(Mg24, O16, scratch);

    scratch = screened_rates(k_he4_ne20_to_mg24)*Y(He4)*state.rho;
    jac.set(Mg24, Ne20, scratch);

    scratch = screened_rates(k_p_na23_to_mg24)*Y(H1)*state.rho;
    jac.set(Mg24, Na23, scratch);

    scratch = -screened_rates(k_he4_mg24_to_p_al27_derived)*Y(He4)*state.rho - screened_rates(k_he4_mg24_to_si28)*Y(He4)*state.rho - screened_rates(k_mg24_to_he4_ne20_derived) - screened_rates(k_mg24_to_p_na23_derived);
    jac.set(Mg24, Mg24, scratch);

    scratch = screened_rates(k_p_al27_to_he4_mg24)*Y(H1)*state.rho;
    jac.set(Mg24, Al27, scratch);

    scratch = screened_rates(k_si28_to_he4_mg24_derived);
    jac.set(Mg24, Si28, scratch);

    scratch = -screened_rates(k_p_al27_to_he4_mg24)*Y(Al27)*state.rho - screened_rates(k_p_al27_to_si28)*Y(Al27)*state.rho;
    jac.set(Al27, H1, scratch);

    scratch = screened_rates(k_he4_mg24_to_p_al27_derived)*Y(Mg24)*state.rho;
    jac.set(Al27, He4, scratch);

    scratch = screened_rates(k_c12_o16_to_p_al27)*Y(O16)*state.rho;
    jac.set(Al27, C12, scratch);

    scratch = screened_rates(k_c12_o16_to_p_al27)*Y(C12)*state.rho;
    jac.set(Al27, O16, scratch);

    scratch = screened_rates(k_he4_mg24_to_p_al27_derived)*Y(He4)*state.rho;
    jac.set(Al27, Mg24, scratch);

    scratch = -screened_rates(k_p_al27_to_he4_mg24)*Y(H1)*state.rho - screened_rates(k_p_al27_to_si28)*Y(H1)*state.rho;
    jac.set(Al27, Al27, scratch);

    scratch = screened_rates(k_si28_to_p_al27_derived);
    jac.set(Al27, Si28, scratch);

    scratch = screened_rates(k_p_al27_to_si28)*Y(Al27)*state.rho + screened_rates(k_p_p31_to_he4_si28)*Y(P31)*state.rho;
    jac.set(Si28, H1, scratch);

    scratch = screened_rates(k_he4_mg24_to_si28)*Y(Mg24)*state.rho - screened_rates(k_he4_si28_to_p_p31_derived)*Y(Si28)*state.rho - screened_rates(k_he4_si28_to_s32)*Y(Si28)*state.rho;
    jac.set(Si28, He4, scratch);

    scratch = screened_rates(k_c12_o16_to_si28_modified)*Y(O16)*state.rho;
    jac.set(Si28, C12, scratch);

    scratch = screened_rates(k_c12_o16_to_si28_modified)*Y(C12)*state.rho + 1.0*screened_rates(k_o16_o16_to_he4_si28)*Y(O16)*state.rho;
    jac.set(Si28, O16, scratch);

    scratch = screened_rates(k_he4_mg24_to_si28)*Y(He4)*state.rho;
    jac.set(Si28, Mg24, scratch);

    scratch = screened_rates(k_p_al27_to_si28)*Y(H1)*state.rho;
    jac.set(Si28, Al27, scratch);

    scratch = -screened_rates(k_he4_si28_to_p_p31_derived)*Y(He4)*state.rho - screened_rates(k_he4_si28_to_s32)*Y(He4)*state.rho - screened_rates(k_si28_to_he4_mg24_derived) - screened_rates(k_si28_to_p_al27_derived);
    jac.set(Si28, Si28, scratch);

    scratch = screened_rates(k_p_p31_to_he4_si28)*Y(H1)*state.rho;
    jac.set(Si28, P31, scratch);

    scratch = screened_rates(k_s32_to_he4_si28_derived);
    jac.set(Si28, S32, scratch);

    scratch = -screened_rates(k_p_p31_to_he4_si28)*Y(P31)*state.rho - screened_rates(k_p_p31_to_s32)*Y(P31)*state.rho;
    jac.set(P31, H1, scratch);

    scratch = screened_rates(k_he4_si28_to_p_p31_derived)*Y(Si28)*state.rho;
    jac.set(P31, He4, scratch);

    scratch = 1.0*screened_rates(k_o16_o16_to_p_p31)*Y(O16)*state.rho;
    jac.set(P31, O16, scratch);

    scratch = screened_rates(k_he4_si28_to_p_p31_derived)*Y(He4)*state.rho;
    jac.set(P31, Si28, scratch);

    scratch = -screened_rates(k_p_p31_to_he4_si28)*Y(H1)*state.rho - screened_rates(k_p_p31_to_s32)*Y(H1)*state.rho;
    jac.set(P31, P31, scratch);

    scratch = screened_rates(k_s32_to_p_p31_derived);
    jac.set(P31, S32, scratch);

    scratch = screened_rates(k_p_p31_to_s32)*Y(P31)*state.rho;
    jac.set(S32, H1, scratch);

    scratch = screened_rates(k_he4_si28_to_s32)*Y(Si28)*state.rho - screened_rates(k_s32_he4_to_ar36_approx)*Y(S32)*state.rho;
    jac.set(S32, He4, scratch);

    scratch = 1.0*screened_rates(k_o16_o16_to_s32_modified)*Y(O16)*state.rho;
    jac.set(S32, O16, scratch);

    scratch = screened_rates(k_he4_si28_to_s32)*Y(He4)*state.rho;
    jac.set(S32, Si28, scratch);

    scratch = screened_rates(k_p_p31_to_s32)*Y(H1)*state.rho;
    jac.set(S32, P31, scratch);

    scratch = -screened_rates(k_s32_he4_to_ar36_approx)*Y(He4)*state.rho - screened_rates(k_s32_to_he4_si28_derived) - screened_rates(k_s32_to_p_p31_derived);
    jac.set(S32, S32, scratch);

    scratch = screened_rates(k_ar36_to_s32_he4_approx);
    jac.set(S32, Ar36, scratch);

    scratch = -screened_rates(k_ar36_he4_to_ca40_approx)*Y(Ar36)*state.rho + screened_rates(k_s32_he4_to_ar36_approx)*Y(S32)*state.rho;
    jac.set(Ar36, He4, scratch);

    scratch = screened_rates(k_s32_he4_to_ar36_approx)*Y(He4)*state.rho;
    jac.set(Ar36, S32, scratch);

    scratch = -screened_rates(k_ar36_he4_to_ca40_approx)*Y(He4)*state.rho - screened_rates(k_ar36_to_s32_he4_approx);
    jac.set(Ar36, Ar36, scratch);

    scratch = screened_rates(k_ca40_to_ar36_he4_approx);
    jac.set(Ar36, Ca40, scratch);

    scratch = screened_rates(k_ar36_he4_to_ca40_approx)*Y(Ar36)*state.rho - screened_rates(k_ca40_he4_to_ti44_approx)*Y(Ca40)*state.rho;
    jac.set(Ca40, He4, scratch);

    scratch = screened_rates(k_ar36_he4_to_ca40_approx)*Y(He4)*state.rho;
    jac.set(Ca40, Ar36, scratch);

    scratch = -screened_rates(k_ca40_he4_to_ti44_approx)*Y(He4)*state.rho - screened_rates(k_ca40_to_ar36_he4_approx);
    jac.set(Ca40, Ca40, scratch);

    scratch = screened_rates(k_ti44_to_ca40_he4_approx);
    jac.set(Ca40, Ti44, scratch);

    scratch = screened_rates(k_ca40_he4_to_ti44_approx)*Y(Ca40)*state.rho - screened_rates(k_ti44_he4_to_cr48_approx)*Y(Ti44)*state.rho;
    jac.set(Ti44, He4, scratch);

    scratch = screened_rates(k_ca40_he4_to_ti44_approx)*Y(He4)*state.rho;
    jac.set(Ti44, Ca40, scratch);

    scratch = -screened_rates(k_ti44_he4_to_cr48_approx)*Y(He4)*state.rho - screened_rates(k_ti44_to_ca40_he4_approx);
    jac.set(Ti44, Ti44, scratch);

    scratch = screened_rates(k_cr48_to_ti44_he4_approx);
    jac.set(Ti44, Cr48, scratch);

    scratch = screened_rates(k_p_mn51_to_he4_cr48_derived)*Y(Mn51)*state.rho;
    jac.set(Cr48, H1, scratch);

    scratch = -screened_rates(k_he4_cr48_to_fe52)*Y(Cr48)*state.rho - screened_rates(k_he4_cr48_to_p_mn51)*Y(Cr48)*state.rho + screened_rates(k_ti44_he4_to_cr48_approx)*Y(Ti44)*state.rho;
    jac.set(Cr48, He4, scratch);

    scratch = screened_rates(k_ti44_he4_to_cr48_approx)*Y(He4)*state.rho;
    jac.set(Cr48, Ti44, scratch);

    scratch = -screened_rates(k_cr48_to_ti44_he4_approx) - screened_rates(k_he4_cr48_to_fe52)*Y(He4)*state.rho - screened_rates(k_he4_cr48_to_p_mn51)*Y(He4)*state.rho;
    jac.set(Cr48, Cr48, scratch);

    scratch = screened_rates(k_p_mn51_to_he4_cr48_derived)*Y(H1)*state.rho;
    jac.set(Cr48, Mn51, scratch);

    scratch = screened_rates(k_fe52_to_he4_cr48_derived);
    jac.set(Cr48, Fe52, scratch);

    scratch = screened_rates(k_p_fe54_to_he4_mn51_derived)*Y(Fe54)*state.rho - screened_rates(k_p_mn51_to_fe52)*Y(Mn51)*state.rho - screened_rates(k_p_mn51_to_he4_cr48_derived)*Y(Mn51)*state.rho;
    jac.set(Mn51, H1, scratch);

    scratch = screened_rates(k_he4_cr48_to_p_mn51)*Y(Cr48)*state.rho - screened_rates(k_he4_mn51_to_co55)*Y(Mn51)*state.rho - screened_rates(k_he4_mn51_to_p_fe54)*Y(Mn51)*state.rho;
    jac.set(Mn51, He4, scratch);

    scratch = screened_rates(k_he4_cr48_to_p_mn51)*Y(He4)*state.rho;
    jac.set(Mn51, Cr48, scratch);

    scratch = -screened_rates(k_he4_mn51_to_co55)*Y(He4)*state.rho - screened_rates(k_he4_mn51_to_p_fe54)*Y(He4)*state.rho - screened_rates(k_p_mn51_to_fe52)*Y(H1)*state.rho - screened_rates(k_p_mn51_to_he4_cr48_derived)*Y(H1)*state.rho;
    jac.set(Mn51, Mn51, scratch);

    scratch = screened_rates(k_fe52_to_p_mn51_derived);
    jac.set(Mn51, Fe52, scratch);

    scratch = screened_rates(k_p_fe54_to_he4_mn51_derived)*Y(H1)*state.rho;
    jac.set(Mn51, Fe54, scratch);

    scratch = screened_rates(k_co55_to_he4_mn51_derived);
    jac.set(Mn51, Co55, scratch);

    scratch = screened_rates(k_n_fe55_to_p_mn55)*Y(Fe55)*state.rho;
    jac.set(Mn55, N, scratch);

    scratch = -screened_rates(k_p_mn55_to_fe56)*Y(Mn55)*state.rho - screened_rates(k_p_mn55_to_n_fe55_derived)*Y(Mn55)*state.rho;
    jac.set(Mn55, H1, scratch);

    scratch = -screened_rates(k_mn55_to_fe55) - screened_rates(k_p_mn55_to_fe56)*Y(H1)*state.rho - screened_rates(k_p_mn55_to_n_fe55_derived)*Y(H1)*state.rho;
    jac.set(Mn55, Mn55, scratch);

    scratch = screened_rates(k_fe55_to_mn55) + screened_rates(k_n_fe55_to_p_mn55)*Y(N)*state.rho;
    jac.set(Mn55, Fe55, scratch);

    scratch = screened_rates(k_fe56_to_p_mn55_derived);
    jac.set(Mn55, Fe56, scratch);

    scratch = -screened_rates(k_n_fe52_to_fe53)*Y(Fe52)*state.rho;
    jac.set(Fe52, N, scratch);

    scratch = screened_rates(k_p_co55_to_he4_fe52_derived)*Y(Co55)*state.rho + screened_rates(k_p_mn51_to_fe52)*Y(Mn51)*state.rho;
    jac.set(Fe52, H1, scratch);

    scratch = screened_rates(k_he4_cr48_to_fe52)*Y(Cr48)*state.rho - screened_rates(k_he4_fe52_to_ni56)*Y(Fe52)*state.rho - screened_rates(k_he4_fe52_to_p_co55)*Y(Fe52)*state.rho;
    jac.set(Fe52, He4, scratch);

    scratch = screened_rates(k_he4_cr48_to_fe52)*Y(He4)*state.rho;
    jac.set(Fe52, Cr48, scratch);

    scratch = screened_rates(k_p_mn51_to_fe52)*Y(H1)*state.rho;
    jac.set(Fe52, Mn51, scratch);

    scratch = -screened_rates(k_fe52_to_he4_cr48_derived) - screened_rates(k_fe52_to_p_mn51_derived) - screened_rates(k_he4_fe52_to_ni56)*Y(He4)*state.rho - screened_rates(k_he4_fe52_to_p_co55)*Y(He4)*state.rho - screened_rates(k_n_fe52_to_fe53)*Y(N)*state.rho;
    jac.set(Fe52, Fe52, scratch);

    scratch = screened_rates(k_fe53_to_n_fe52_derived);
    jac.set(Fe52, Fe53, scratch);

    scratch = screened_rates(k_p_co55_to_he4_fe52_derived)*Y(H1)*state.rho;
    jac.set(Fe52, Co55, scratch);

    scratch = screened_rates(k_ni56_to_he4_fe52_derived);
    jac.set(Fe52, Ni56, scratch);

    scratch = screened_rates(k_n_fe52_to_fe53)*Y(Fe52)*state.rho - screened_rates(k_n_fe53_to_fe54)*Y(Fe53)*state.rho + screened_rates(k_n_ni56_to_he4_fe53)*Y(Ni56)*state.rho;
    jac.set(Fe53, N, scratch);

    scratch = screened_rates(k_p_co56_to_he4_fe53_derived)*Y(Co56)*state.rho;
    jac.set(Fe53, H1, scratch);

    scratch = -screened_rates(k_he4_fe53_to_n_ni56_derived)*Y(Fe53)*state.rho - screened_rates(k_he4_fe53_to_ni57)*Y(Fe53)*state.rho - screened_rates(k_he4_fe53_to_p_co56)*Y(Fe53)*state.rho;
    jac.set(Fe53, He4, scratch);

    scratch = screened_rates(k_n_fe52_to_fe53)*Y(N)*state.rho;
    jac.set(Fe53, Fe52, scratch);

    scratch = -screened_rates(k_fe53_to_n_fe52_derived) - screened_rates(k_he4_fe53_to_n_ni56_derived)*Y(He4)*state.rho - screened_rates(k_he4_fe53_to_ni57)*Y(He4)*state.rho - screened_rates(k_he4_fe53_to_p_co56)*Y(He4)*state.rho - screened_rates(k_n_fe53_to_fe54)*Y(N)*state.rho;
    jac.set(Fe53, Fe53, scratch);

    scratch = screened_rates(k_fe54_to_n_fe53_derived);
    jac.set(Fe53, Fe54, scratch);

    scratch = screened_rates(k_p_co56_to_he4_fe53_derived)*Y(H1)*state.rho;
    jac.set(Fe53, Co56, scratch);

    scratch = screened_rates(k_n_ni56_to_he4_fe53)*Y(N)*state.rho;
    jac.set(Fe53, Ni56, scratch);

    scratch = screened_rates(k_ni57_to_he4_fe53_derived);
    jac.set(Fe53, Ni57, scratch);

    scratch = screened_rates(k_n_fe53_to_fe54)*Y(Fe53)*state.rho - screened_rates(k_n_fe54_to_fe55)*Y(Fe54)*state.rho + screened_rates(k_n_ni57_to_he4_fe54)*Y(Ni57)*state.rho;
    jac.set(Fe54, N, scratch);

    scratch = screened_rates(k_p_co57_to_he4_fe54)*Y(Co57)*state.rho - screened_rates(k_p_fe54_to_co55)*Y(Fe54)*state.rho - screened_rates(k_p_fe54_to_he4_mn51_derived)*Y(Fe54)*state.rho;
    jac.set(Fe54, H1, scratch);

    scratch = -screened_rates(k_he4_fe54_to_n_ni57_derived)*Y(Fe54)*state.rho - screened_rates(k_he4_fe54_to_ni58)*Y(Fe54)*state.rho - screened_rates(k_he4_fe54_to_p_co57_derived)*Y(Fe54)*state.rho + screened_rates(k_he4_mn51_to_p_fe54)*Y(Mn51)*state.rho;
    jac.set(Fe54, He4, scratch);

    scratch = screened_rates(k_he4_mn51_to_p_fe54)*Y(He4)*state.rho;
    jac.set(Fe54, Mn51, scratch);

    scratch = screened_rates(k_n_fe53_to_fe54)*Y(N)*state.rho;
    jac.set(Fe54, Fe53, scratch);

    scratch = -screened_rates(k_fe54_to_n_fe53_derived) - screened_rates(k_he4_fe54_to_n_ni57_derived)*Y(He4)*state.rho - screened_rates(k_he4_fe54_to_ni58)*Y(He4)*state.rho - screened_rates(k_he4_fe54_to_p_co57_derived)*Y(He4)*state.rho - screened_rates(k_n_fe54_to_fe55)*Y(N)*state.rho - screened_rates(k_p_fe54_to_co55)*Y(H1)*state.rho - screened_rates(k_p_fe54_to_he4_mn51_derived)*Y(H1)*state.rho;
    jac.set(Fe54, Fe54, scratch);

    scratch = screened_rates(k_fe55_to_n_fe54_derived);
    jac.set(Fe54, Fe55, scratch);

    scratch = screened_rates(k_co55_to_p_fe54_derived);
    jac.set(Fe54, Co55, scratch);

    scratch = screened_rates(k_p_co57_to_he4_fe54)*Y(H1)*state.rho;
    jac.set(Fe54, Co57, scratch);

    scratch = screened_rates(k_n_ni57_to_he4_fe54)*Y(N)*state.rho;
    jac.set(Fe54, Ni57, scratch);

    scratch = screened_rates(k_ni58_to_he4_fe54_derived);
    jac.set(Fe54, Ni58, scratch);

    scratch = screened_rates(k_n_co55_to_p_fe55)*Y(Co55)*state.rho + screened_rates(k_n_fe54_to_fe55)*Y(Fe54)*state.rho - screened_rates(k_n_fe55_to_fe56)*Y(Fe55)*state.rho - screened_rates(k_n_fe55_to_p_mn55)*Y(Fe55)*state.rho + screened_rates(k_n_ni58_to_he4_fe55)*Y(Ni58)*state.rho;
    jac.set(Fe55, N, scratch);

    scratch = -screened_rates(k_p_fe55_to_co56)*Y(Fe55)*state.rho - screened_rates(k_p_fe55_to_n_co55_derived)*Y(Fe55)*state.rho + screened_rates(k_p_mn55_to_n_fe55_derived)*Y(Mn55)*state.rho;
    jac.set(Fe55, H1, scratch);

    scratch = -screened_rates(k_he4_fe55_to_n_ni58_derived)*Y(Fe55)*state.rho;
    jac.set(Fe55, He4, scratch);

    scratch = screened_rates(k_mn55_to_fe55) + screened_rates(k_p_mn55_to_n_fe55_derived)*Y(H1)*state.rho;
    jac.set(Fe55, Mn55, scratch);

    scratch = screened_rates(k_n_fe54_to_fe55)*Y(N)*state.rho;
    jac.set(Fe55, Fe54, scratch);

    scratch = -screened_rates(k_fe55_to_co55) - screened_rates(k_fe55_to_mn55) - screened_rates(k_fe55_to_n_fe54_derived) - screened_rates(k_he4_fe55_to_n_ni58_derived)*Y(He4)*state.rho - screened_rates(k_n_fe55_to_fe56)*Y(N)*state.rho - screened_rates(k_n_fe55_to_p_mn55)*Y(N)*state.rho - screened_rates(k_p_fe55_to_co56)*Y(H1)*state.rho - screened_rates(k_p_fe55_to_n_co55_derived)*Y(H1)*state.rho;
    jac.set(Fe55, Fe55, scratch);

    scratch = screened_rates(k_fe56_to_n_fe55_derived);
    jac.set(Fe55, Fe56, scratch);

    scratch = screened_rates(k_co55_to_fe55) + screened_rates(k_n_co55_to_p_fe55)*Y(N)*state.rho;
    jac.set(Fe55, Co55, scratch);

    scratch = screened_rates(k_co56_to_p_fe55_derived);
    jac.set(Fe55, Co56, scratch);

    scratch = screened_rates(k_n_ni58_to_he4_fe55)*Y(N)*state.rho;
    jac.set(Fe55, Ni58, scratch);

    scratch = screened_rates(k_n_co56_to_p_fe56)*Y(Co56)*state.rho + screened_rates(k_n_fe55_to_fe56)*Y(Fe55)*state.rho;
    jac.set(Fe56, N, scratch);

    scratch = -screened_rates(k_p_fe56_to_co57)*Y(Fe56)*state.rho - screened_rates(k_p_fe56_to_n_co56_derived)*Y(Fe56)*state.rho + screened_rates(k_p_mn55_to_fe56)*Y(Mn55)*state.rho;
    jac.set(Fe56, H1, scratch);

    scratch = screened_rates(k_p_mn55_to_fe56)*Y(H1)*state.rho;
    jac.set(Fe56, Mn55, scratch);

    scratch = screened_rates(k_n_fe55_to_fe56)*Y(N)*state.rho;
    jac.set(Fe56, Fe55, scratch);

    scratch = -screened_rates(k_fe56_to_co56) - screened_rates(k_fe56_to_n_fe55_derived) - screened_rates(k_fe56_to_p_mn55_derived) - screened_rates(k_p_fe56_to_co57)*Y(H1)*state.rho - screened_rates(k_p_fe56_to_n_co56_derived)*Y(H1)*state.rho;
    jac.set(Fe56, Fe56, scratch);

    scratch = screened_rates(k_co56_to_fe56) + screened_rates(k_n_co56_to_p_fe56)*Y(N)*state.rho;
    jac.set(Fe56, Co56, scratch);

    scratch = screened_rates(k_co57_to_p_fe56_derived);
    jac.set(Fe56, Co57, scratch);

    scratch = -screened_rates(k_n_co55_to_co56)*Y(Co55)*state.rho - screened_rates(k_n_co55_to_p_fe55)*Y(Co55)*state.rho;
    jac.set(Co55, N, scratch);

    scratch = -screened_rates(k_p_co55_to_he4_fe52_derived)*Y(Co55)*state.rho - screened_rates(k_p_co55_to_ni56)*Y(Co55)*state.rho + screened_rates(k_p_fe54_to_co55)*Y(Fe54)*state.rho + screened_rates(k_p_fe55_to_n_co55_derived)*Y(Fe55)*state.rho + screened_rates(k_p_ni58_to_he4_co55_derived)*Y(Ni58)*state.rho;
    jac.set(Co55, H1, scratch);

    scratch = -screened_rates(k_he4_co55_to_cu59)*Y(Co55)*state.rho - screened_rates(k_he4_co55_to_p_ni58)*Y(Co55)*state.rho + screened_rates(k_he4_fe52_to_p_co55)*Y(Fe52)*state.rho + screened_rates(k_he4_mn51_to_co55)*Y(Mn51)*state.rho;
    jac.set(Co55, He4, scratch);

    scratch = screened_rates(k_he4_mn51_to_co55)*Y(He4)*state.rho;
    jac.set(Co55, Mn51, scratch);

    scratch = screened_rates(k_he4_fe52_to_p_co55)*Y(He4)*state.rho;
    jac.set(Co55, Fe52, scratch);

    scratch = screened_rates(k_p_fe54_to_co55)*Y(H1)*state.rho;
    jac.set(Co55, Fe54, scratch);

    scratch = screened_rates(k_fe55_to_co55) + screened_rates(k_p_fe55_to_n_co55_derived)*Y(H1)*state.rho;
    jac.set(Co55, Fe55, scratch);

    scratch = -screened_rates(k_co55_to_fe55) - screened_rates(k_co55_to_he4_mn51_derived) - screened_rates(k_co55_to_p_fe54_derived) - screened_rates(k_he4_co55_to_cu59)*Y(He4)*state.rho - screened_rates(k_he4_co55_to_p_ni58)*Y(He4)*state.rho - screened_rates(k_n_co55_to_co56)*Y(N)*state.rho - screened_rates(k_n_co55_to_p_fe55)*Y(N)*state.rho - screened_rates(k_p_co55_to_he4_fe52_derived)*Y(H1)*state.rho - screened_rates(k_p_co55_to_ni56)*Y(H1)*state.rho;
    jac.set(Co55, Co55, scratch);

    scratch = screened_rates(k_co56_to_n_co55_derived);
    jac.set(Co55, Co56, scratch);

    scratch = screened_rates(k_ni56_to_p_co55_derived);
    jac.set(Co55, Ni56, scratch);

    scratch = screened_rates(k_p_ni58_to_he4_co55_derived)*Y(H1)*state.rho;
    jac.set(Co55, Ni58, scratch);

    scratch = screened_rates(k_cu59_to_he4_co55_derived);
    jac.set(Co55, Cu59, scratch);

    scratch = screened_rates(k_n_co55_to_co56)*Y(Co55)*state.rho - screened_rates(k_n_co56_to_co57)*Y(Co56)*state.rho - screened_rates(k_n_co56_to_p_fe56)*Y(Co56)*state.rho + screened_rates(k_n_cu59_to_he4_co56)*Y(Cu59)*state.rho + screened_rates(k_n_ni56_to_p_co56)*Y(Ni56)*state.rho;
    jac.set(Co56, N, scratch);

    scratch = -screened_rates(k_p_co56_to_he4_fe53_derived)*Y(Co56)*state.rho - screened_rates(k_p_co56_to_n_ni56_derived)*Y(Co56)*state.rho - screened_rates(k_p_co56_to_ni57)*Y(Co56)*state.rho + screened_rates(k_p_fe55_to_co56)*Y(Fe55)*state.rho + screened_rates(k_p_fe56_to_n_co56_derived)*Y(Fe56)*state.rho;
    jac.set(Co56, H1, scratch);

    scratch = -screened_rates(k_he4_co56_to_n_cu59_derived)*Y(Co56)*state.rho + screened_rates(k_he4_fe53_to_p_co56)*Y(Fe53)*state.rho;
    jac.set(Co56, He4, scratch);

    scratch = screened_rates(k_he4_fe53_to_p_co56)*Y(He4)*state.rho;
    jac.set(Co56, Fe53, scratch);

    scratch = screened_rates(k_p_fe55_to_co56)*Y(H1)*state.rho;
    jac.set(Co56, Fe55, scratch);

    scratch = screened_rates(k_fe56_to_co56) + screened_rates(k_p_fe56_to_n_co56_derived)*Y(H1)*state.rho;
    jac.set(Co56, Fe56, scratch);

    scratch = screened_rates(k_n_co55_to_co56)*Y(N)*state.rho;
    jac.set(Co56, Co55, scratch);

    scratch = -screened_rates(k_co56_to_fe56) - screened_rates(k_co56_to_n_co55_derived) - screened_rates(k_co56_to_ni56) - screened_rates(k_co56_to_p_fe55_derived) - screened_rates(k_he4_co56_to_n_cu59_derived)*Y(He4)*state.rho - screened_rates(k_n_co56_to_co57)*Y(N)*state.rho - screened_rates(k_n_co56_to_p_fe56)*Y(N)*state.rho - screened_rates(k_p_co56_to_he4_fe53_derived)*Y(H1)*state.rho - screened_rates(k_p_co56_to_n_ni56_derived)*Y(H1)*state.rho - screened_rates(k_p_co56_to_ni57)*Y(H1)*state.rho;
    jac.set(Co56, Co56, scratch);

    scratch = screened_rates(k_co57_to_n_co56_derived);
    jac.set(Co56, Co57, scratch);

    scratch = screened_rates(k_n_ni56_to_p_co56)*Y(N)*state.rho + screened_rates(k_ni56_to_co56);
    jac.set(Co56, Ni56, scratch);

    scratch = screened_rates(k_ni57_to_p_co56_derived);
    jac.set(Co56, Ni57, scratch);

    scratch = screened_rates(k_n_cu59_to_he4_co56)*Y(N)*state.rho;
    jac.set(Co56, Cu59, scratch);

    scratch = screened_rates(k_n_co56_to_co57)*Y(Co56)*state.rho + screened_rates(k_n_ni57_to_p_co57)*Y(Ni57)*state.rho;
    jac.set(Co57, N, scratch);

    scratch = -screened_rates(k_p_co57_to_he4_fe54)*Y(Co57)*state.rho - screened_rates(k_p_co57_to_n_ni57_derived)*Y(Co57)*state.rho - screened_rates(k_p_co57_to_ni58)*Y(Co57)*state.rho + screened_rates(k_p_fe56_to_co57)*Y(Fe56)*state.rho;
    jac.set(Co57, H1, scratch);

    scratch = screened_rates(k_he4_fe54_to_p_co57_derived)*Y(Fe54)*state.rho;
    jac.set(Co57, He4, scratch);

    scratch = screened_rates(k_he4_fe54_to_p_co57_derived)*Y(He4)*state.rho;
    jac.set(Co57, Fe54, scratch);

    scratch = screened_rates(k_p_fe56_to_co57)*Y(H1)*state.rho;
    jac.set(Co57, Fe56, scratch);

    scratch = screened_rates(k_n_co56_to_co57)*Y(N)*state.rho;
    jac.set(Co57, Co56, scratch);

    scratch = -screened_rates(k_co57_to_n_co56_derived) - screened_rates(k_co57_to_ni57) - screened_rates(k_co57_to_p_fe56_derived) - screened_rates(k_p_co57_to_he4_fe54)*Y(H1)*state.rho - screened_rates(k_p_co57_to_n_ni57_derived)*Y(H1)*state.rho - screened_rates(k_p_co57_to_ni58)*Y(H1)*state.rho;
    jac.set(Co57, Co57, scratch);

    scratch = screened_rates(k_n_ni57_to_p_co57)*Y(N)*state.rho + screened_rates(k_ni57_to_co57);
    jac.set(Co57, Ni57, scratch);

    scratch = screened_rates(k_ni58_to_p_co57_derived);
    jac.set(Co57, Ni58, scratch);

    scratch = -screened_rates(k_n_ni56_to_he4_fe53)*Y(Ni56)*state.rho - screened_rates(k_n_ni56_to_ni57)*Y(Ni56)*state.rho - screened_rates(k_n_ni56_to_p_co56)*Y(Ni56)*state.rho;
    jac.set(Ni56, N, scratch);

    scratch = screened_rates(k_p_co55_to_ni56)*Y(Co55)*state.rho + screened_rates(k_p_co56_to_n_ni56_derived)*Y(Co56)*state.rho + screened_rates(k_p_cu59_to_he4_ni56)*Y(Cu59)*state.rho;
    jac.set(Ni56, H1, scratch);

    scratch = screened_rates(k_he4_fe52_to_ni56)*Y(Fe52)*state.rho + screened_rates(k_he4_fe53_to_n_ni56_derived)*Y(Fe53)*state.rho - screened_rates(k_he4_ni56_to_p_cu59_derived)*Y(Ni56)*state.rho - screened_rates(k_he4_ni56_to_zn60)*Y(Ni56)*state.rho;
    jac.set(Ni56, He4, scratch);

    scratch = screened_rates(k_he4_fe52_to_ni56)*Y(He4)*state.rho;
    jac.set(Ni56, Fe52, scratch);

    scratch = screened_rates(k_he4_fe53_to_n_ni56_derived)*Y(He4)*state.rho;
    jac.set(Ni56, Fe53, scratch);

    scratch = screened_rates(k_p_co55_to_ni56)*Y(H1)*state.rho;
    jac.set(Ni56, Co55, scratch);

    scratch = screened_rates(k_co56_to_ni56) + screened_rates(k_p_co56_to_n_ni56_derived)*Y(H1)*state.rho;
    jac.set(Ni56, Co56, scratch);

    scratch = -screened_rates(k_he4_ni56_to_p_cu59_derived)*Y(He4)*state.rho - screened_rates(k_he4_ni56_to_zn60)*Y(He4)*state.rho - screened_rates(k_n_ni56_to_he4_fe53)*Y(N)*state.rho - screened_rates(k_n_ni56_to_ni57)*Y(N)*state.rho - screened_rates(k_n_ni56_to_p_co56)*Y(N)*state.rho - screened_rates(k_ni56_to_co56) - screened_rates(k_ni56_to_he4_fe52_derived) - screened_rates(k_ni56_to_p_co55_derived);
    jac.set(Ni56, Ni56, scratch);

    scratch = screened_rates(k_ni57_to_n_ni56_derived);
    jac.set(Ni56, Ni57, scratch);

    scratch = screened_rates(k_p_cu59_to_he4_ni56)*Y(H1)*state.rho;
    jac.set(Ni56, Cu59, scratch);

    scratch = screened_rates(k_zn60_to_he4_ni56_derived);
    jac.set(Ni56, Zn60, scratch);

    scratch = screened_rates(k_n_ni56_to_ni57)*Y(Ni56)*state.rho - screened_rates(k_n_ni57_to_he4_fe54)*Y(Ni57)*state.rho - screened_rates(k_n_ni57_to_ni58)*Y(Ni57)*state.rho - screened_rates(k_n_ni57_to_p_co57)*Y(Ni57)*state.rho + screened_rates(k_n_zn60_to_he4_ni57)*Y(Zn60)*state.rho;
    jac.set(Ni57, N, scratch);

    scratch = screened_rates(k_p_co56_to_ni57)*Y(Co56)*state.rho + screened_rates(k_p_co57_to_n_ni57_derived)*Y(Co57)*state.rho;
    jac.set(Ni57, H1, scratch);

    scratch = screened_rates(k_he4_fe53_to_ni57)*Y(Fe53)*state.rho + screened_rates(k_he4_fe54_to_n_ni57_derived)*Y(Fe54)*state.rho - screened_rates(k_he4_ni57_to_n_zn60_derived)*Y(Ni57)*state.rho;
    jac.set(Ni57, He4, scratch);

    scratch = screened_rates(k_he4_fe53_to_ni57)*Y(He4)*state.rho;
    jac.set(Ni57, Fe53, scratch);

    scratch = screened_rates(k_he4_fe54_to_n_ni57_derived)*Y(He4)*state.rho;
    jac.set(Ni57, Fe54, scratch);

    scratch = screened_rates(k_p_co56_to_ni57)*Y(H1)*state.rho;
    jac.set(Ni57, Co56, scratch);

    scratch = screened_rates(k_co57_to_ni57) + screened_rates(k_p_co57_to_n_ni57_derived)*Y(H1)*state.rho;
    jac.set(Ni57, Co57, scratch);

    scratch = screened_rates(k_n_ni56_to_ni57)*Y(N)*state.rho;
    jac.set(Ni57, Ni56, scratch);

    scratch = -screened_rates(k_he4_ni57_to_n_zn60_derived)*Y(He4)*state.rho - screened_rates(k_n_ni57_to_he4_fe54)*Y(N)*state.rho - screened_rates(k_n_ni57_to_ni58)*Y(N)*state.rho - screened_rates(k_n_ni57_to_p_co57)*Y(N)*state.rho - screened_rates(k_ni57_to_co57) - screened_rates(k_ni57_to_he4_fe53_derived) - screened_rates(k_ni57_to_n_ni56_derived) - screened_rates(k_ni57_to_p_co56_derived);
    jac.set(Ni57, Ni57, scratch);

    scratch = screened_rates(k_ni58_to_n_ni57_derived);
    jac.set(Ni57, Ni58, scratch);

    scratch = screened_rates(k_n_zn60_to_he4_ni57)*Y(N)*state.rho;
    jac.set(Ni57, Zn60, scratch);

    scratch = screened_rates(k_n_ni57_to_ni58)*Y(Ni57)*state.rho - screened_rates(k_n_ni58_to_he4_fe55)*Y(Ni58)*state.rho;
    jac.set(Ni58, N, scratch);

    scratch = screened_rates(k_p_co57_to_ni58)*Y(Co57)*state.rho - screened_rates(k_p_ni58_to_cu59)*Y(Ni58)*state.rho - screened_rates(k_p_ni58_to_he4_co55_derived)*Y(Ni58)*state.rho;
    jac.set(Ni58, H1, scratch);

    scratch = screened_rates(k_he4_co55_to_p_ni58)*Y(Co55)*state.rho + screened_rates(k_he4_fe54_to_ni58)*Y(Fe54)*state.rho + screened_rates(k_he4_fe55_to_n_ni58_derived)*Y(Fe55)*state.rho;
    jac.set(Ni58, He4, scratch);

    scratch = screened_rates(k_he4_fe54_to_ni58)*Y(He4)*state.rho;
    jac.set(Ni58, Fe54, scratch);

    scratch = screened_rates(k_he4_fe55_to_n_ni58_derived)*Y(He4)*state.rho;
    jac.set(Ni58, Fe55, scratch);

    scratch = screened_rates(k_he4_co55_to_p_ni58)*Y(He4)*state.rho;
    jac.set(Ni58, Co55, scratch);

    scratch = screened_rates(k_p_co57_to_ni58)*Y(H1)*state.rho;
    jac.set(Ni58, Co57, scratch);

    scratch = screened_rates(k_n_ni57_to_ni58)*Y(N)*state.rho;
    jac.set(Ni58, Ni57, scratch);

    scratch = -screened_rates(k_n_ni58_to_he4_fe55)*Y(N)*state.rho - screened_rates(k_ni58_to_he4_fe54_derived) - screened_rates(k_ni58_to_n_ni57_derived) - screened_rates(k_ni58_to_p_co57_derived) - screened_rates(k_p_ni58_to_cu59)*Y(H1)*state.rho - screened_rates(k_p_ni58_to_he4_co55_derived)*Y(H1)*state.rho;
    jac.set(Ni58, Ni58, scratch);

    scratch = screened_rates(k_cu59_to_p_ni58_derived);
    jac.set(Ni58, Cu59, scratch);

    scratch = -screened_rates(k_n_cu59_to_he4_co56)*Y(Cu59)*state.rho;
    jac.set(Cu59, N, scratch);

    scratch = -screened_rates(k_p_cu59_to_he4_ni56)*Y(Cu59)*state.rho - screened_rates(k_p_cu59_to_zn60)*Y(Cu59)*state.rho + screened_rates(k_p_ni58_to_cu59)*Y(Ni58)*state.rho;
    jac.set(Cu59, H1, scratch);

    scratch = screened_rates(k_he4_co55_to_cu59)*Y(Co55)*state.rho + screened_rates(k_he4_co56_to_n_cu59_derived)*Y(Co56)*state.rho + screened_rates(k_he4_ni56_to_p_cu59_derived)*Y(Ni56)*state.rho;
    jac.set(Cu59, He4, scratch);

    scratch = screened_rates(k_he4_co55_to_cu59)*Y(He4)*state.rho;
    jac.set(Cu59, Co55, scratch);

    scratch = screened_rates(k_he4_co56_to_n_cu59_derived)*Y(He4)*state.rho;
    jac.set(Cu59, Co56, scratch);

    scratch = screened_rates(k_he4_ni56_to_p_cu59_derived)*Y(He4)*state.rho;
    jac.set(Cu59, Ni56, scratch);

    scratch = screened_rates(k_p_ni58_to_cu59)*Y(H1)*state.rho;
    jac.set(Cu59, Ni58, scratch);

    scratch = -screened_rates(k_cu59_to_he4_co55_derived) - screened_rates(k_cu59_to_p_ni58_derived) - screened_rates(k_n_cu59_to_he4_co56)*Y(N)*state.rho - screened_rates(k_p_cu59_to_he4_ni56)*Y(H1)*state.rho - screened_rates(k_p_cu59_to_zn60)*Y(H1)*state.rho;
    jac.set(Cu59, Cu59, scratch);

    scratch = screened_rates(k_zn60_to_p_cu59_derived);
    jac.set(Cu59, Zn60, scratch);

    scratch = -screened_rates(k_n_zn60_to_he4_ni57)*Y(Zn60)*state.rho;
    jac.set(Zn60, N, scratch);

    scratch = screened_rates(k_p_cu59_to_zn60)*Y(Cu59)*state.rho;
    jac.set(Zn60, H1, scratch);

    scratch = screened_rates(k_he4_ni56_to_zn60)*Y(Ni56)*state.rho + screened_rates(k_he4_ni57_to_n_zn60_derived)*Y(Ni57)*state.rho;
    jac.set(Zn60, He4, scratch);

    scratch = screened_rates(k_he4_ni56_to_zn60)*Y(He4)*state.rho;
    jac.set(Zn60, Ni56, scratch);

    scratch = screened_rates(k_he4_ni57_to_n_zn60_derived)*Y(He4)*state.rho;
    jac.set(Zn60, Ni57, scratch);

    scratch = screened_rates(k_p_cu59_to_zn60)*Y(H1)*state.rho;
    jac.set(Zn60, Cu59, scratch);

    scratch = -screened_rates(k_n_zn60_to_he4_ni57)*Y(N)*state.rho - screened_rates(k_zn60_to_he4_ni56_derived) - screened_rates(k_zn60_to_p_cu59_derived);
    jac.set(Zn60, Zn60, scratch);


}



template<class MatrixType>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void actual_jac(const burn_t& state, MatrixType& jac)
{

    // Set molar abundances
    Array1D<Real, 1, NumSpec> Y;
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
        auto jac_slice_2 = [&](int i) -> Real { return jac.get(i, j); };
        ener_gener_rate(jac_slice_2, jac(net_ienuc,j));
    }

    // Account for the thermal neutrino losses

    Real sneut, dsneutdt, dsneutdd, dsnuda, dsnudz;
    constexpr int do_derivatives{1};
    sneut5<do_derivatives>(state.T, state.rho, state.abar, state.zbar, sneut, dsneutdt, dsneutdd, dsnuda, dsnudz);

    for (int j = 1; j <= NumSpec; ++j) {
       Real b1 = (-state.abar * state.abar * dsnuda + (zion[j-1] - state.zbar) * state.abar * dsnudz);
       jac.add(net_ienuc, j, -b1);
    }


    // Evaluate the Jacobian elements with respect to energy by
    // calling the RHS using d(rate) / dT and then transform them
    // to our energy integration variable.

    Array1D<Real, 1, neqs>  yderivs;

    rhs_nuc(state, yderivs, Y, rate_eval.dscreened_rates_dT);

    for (int k = 1; k <= NumSpec; k++) {
        jac.set(k, net_ienuc, temperature_to_energy_jacobian(state, yderivs(k)));
    }


    // finally, d(de/dt)/de

    Real jac_e_T;
    ener_gener_rate(yderivs, jac_e_T);
    jac_e_T -= dsneutdt;
    jac.set(net_ienuc, net_ienuc, temperature_to_energy_jacobian(state, jac_e_T));

}


AMREX_INLINE
void actual_rhs_init () {

    init_tabular();

}


#endif
