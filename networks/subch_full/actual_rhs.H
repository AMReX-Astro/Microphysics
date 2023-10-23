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

    if (disable_p_C12_to_N13) {
        rate_eval.screened_rates(k_p_C12_to_N13) = 0.0;
        if constexpr (std::is_same<T, rate_derivs_t>::value) {
            rate_eval.dscreened_rates_dT(k_p_C12_to_N13) = 0.0;
        }
        rate_eval.screened_rates(k_N13_to_p_C12) = 0.0;
        if constexpr (std::is_same<T, rate_derivs_t>::value) {
            rate_eval.dscreened_rates_dT(k_N13_to_p_C12) = 0.0;
        }
    }

    if (disable_He4_N13_to_p_O16) {
        rate_eval.screened_rates(k_He4_N13_to_p_O16) = 0.0;
        if constexpr (std::is_same<T, rate_derivs_t>::value) {
            rate_eval.dscreened_rates_dT(k_He4_N13_to_p_O16) = 0.0;
        }
        rate_eval.screened_rates(k_p_O16_to_He4_N13) = 0.0;
        if constexpr (std::is_same<T, rate_derivs_t>::value) {
            rate_eval.dscreened_rates_dT(k_p_O16_to_He4_N13) = 0.0;
        }
    }



    // Evaluate screening factors

    Real ratraw, dratraw_dT;
    Real scor, dscor_dt;
    [[maybe_unused]] Real scor2, dscor2_dt;


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 6.0_rt, 12.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_N13_to_p_C12);
    rate_eval.screened_rates(k_N13_to_p_C12) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_N13_to_p_C12);
        rate_eval.dscreened_rates_dT(k_N13_to_p_C12) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_C12_to_N13);
    rate_eval.screened_rates(k_p_C12_to_N13) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_C12_to_N13);
        rate_eval.dscreened_rates_dT(k_p_C12_to_N13) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 6.0_rt, 12.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_O16_to_He4_C12);
    rate_eval.screened_rates(k_O16_to_He4_C12) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_O16_to_He4_C12);
        rate_eval.dscreened_rates_dT(k_O16_to_He4_C12) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_C12_to_O16);
    rate_eval.screened_rates(k_He4_C12_to_O16) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_C12_to_O16);
        rate_eval.dscreened_rates_dT(k_He4_C12_to_O16) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 7.0_rt, 14.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_F18_to_He4_N14);
    rate_eval.screened_rates(k_F18_to_He4_N14) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_F18_to_He4_N14);
        rate_eval.dscreened_rates_dT(k_F18_to_He4_N14) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_N14_to_F18);
    rate_eval.screened_rates(k_He4_N14_to_F18) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_N14_to_F18);
        rate_eval.dscreened_rates_dT(k_He4_N14_to_F18) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 8.0_rt, 16.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_Ne20_to_He4_O16);
    rate_eval.screened_rates(k_Ne20_to_He4_O16) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_Ne20_to_He4_O16);
        rate_eval.dscreened_rates_dT(k_Ne20_to_He4_O16) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_O16_to_Ne20);
    rate_eval.screened_rates(k_He4_O16_to_Ne20) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_O16_to_Ne20);
        rate_eval.dscreened_rates_dT(k_He4_O16_to_Ne20) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 10.0_rt, 21.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_Na22_to_p_Ne21);
    rate_eval.screened_rates(k_Na22_to_p_Ne21) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_Na22_to_p_Ne21);
        rate_eval.dscreened_rates_dT(k_Na22_to_p_Ne21) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ne21_to_Na22);
    rate_eval.screened_rates(k_p_Ne21_to_Na22) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ne21_to_Na22);
        rate_eval.dscreened_rates_dT(k_p_Ne21_to_Na22) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 9.0_rt, 18.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_Na22_to_He4_F18);
    rate_eval.screened_rates(k_Na22_to_He4_F18) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_Na22_to_He4_F18);
        rate_eval.dscreened_rates_dT(k_Na22_to_He4_F18) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_F18_to_Na22);
    rate_eval.screened_rates(k_He4_F18_to_Na22) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_F18_to_Na22);
        rate_eval.dscreened_rates_dT(k_He4_F18_to_Na22) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_F18_to_p_Ne21);
    rate_eval.screened_rates(k_He4_F18_to_p_Ne21) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_F18_to_p_Ne21);
        rate_eval.dscreened_rates_dT(k_He4_F18_to_p_Ne21) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ne21_to_He4_F18);
    rate_eval.screened_rates(k_p_Ne21_to_He4_F18) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ne21_to_He4_F18);
        rate_eval.dscreened_rates_dT(k_p_Ne21_to_He4_F18) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 11.0_rt, 23.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_Mg24_to_p_Na23);
    rate_eval.screened_rates(k_Mg24_to_p_Na23) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_Mg24_to_p_Na23);
        rate_eval.dscreened_rates_dT(k_Mg24_to_p_Na23) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Na23_to_Mg24);
    rate_eval.screened_rates(k_p_Na23_to_Mg24) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Na23_to_Mg24);
        rate_eval.dscreened_rates_dT(k_p_Na23_to_Mg24) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ne20_to_p_Na23);
    rate_eval.screened_rates(k_He4_Ne20_to_p_Na23) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ne20_to_p_Na23);
        rate_eval.dscreened_rates_dT(k_He4_Ne20_to_p_Na23) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Na23_to_He4_Ne20);
    rate_eval.screened_rates(k_p_Na23_to_He4_Ne20) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Na23_to_He4_Ne20);
        rate_eval.dscreened_rates_dT(k_p_Na23_to_He4_Ne20) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 10.0_rt, 20.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_Mg24_to_He4_Ne20);
    rate_eval.screened_rates(k_Mg24_to_He4_Ne20) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_Mg24_to_He4_Ne20);
        rate_eval.dscreened_rates_dT(k_Mg24_to_He4_Ne20) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ne20_to_Mg24);
    rate_eval.screened_rates(k_He4_Ne20_to_Mg24) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ne20_to_Mg24);
        rate_eval.dscreened_rates_dT(k_He4_Ne20_to_Mg24) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 11.0_rt, 23.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_Al27_to_He4_Na23);
    rate_eval.screened_rates(k_Al27_to_He4_Na23) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_Al27_to_He4_Na23);
        rate_eval.dscreened_rates_dT(k_Al27_to_He4_Na23) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Na23_to_Al27);
    rate_eval.screened_rates(k_He4_Na23_to_Al27) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Na23_to_Al27);
        rate_eval.dscreened_rates_dT(k_He4_Na23_to_Al27) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 13.0_rt, 27.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_Si28_to_p_Al27);
    rate_eval.screened_rates(k_Si28_to_p_Al27) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_Si28_to_p_Al27);
        rate_eval.dscreened_rates_dT(k_Si28_to_p_Al27) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Al27_to_Si28);
    rate_eval.screened_rates(k_p_Al27_to_Si28) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Al27_to_Si28);
        rate_eval.dscreened_rates_dT(k_p_Al27_to_Si28) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Mg24_to_p_Al27);
    rate_eval.screened_rates(k_He4_Mg24_to_p_Al27) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mg24_to_p_Al27);
        rate_eval.dscreened_rates_dT(k_He4_Mg24_to_p_Al27) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Al27_to_He4_Mg24);
    rate_eval.screened_rates(k_p_Al27_to_He4_Mg24) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Al27_to_He4_Mg24);
        rate_eval.dscreened_rates_dT(k_p_Al27_to_He4_Mg24) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 12.0_rt, 24.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_Si28_to_He4_Mg24);
    rate_eval.screened_rates(k_Si28_to_He4_Mg24) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_Si28_to_He4_Mg24);
        rate_eval.dscreened_rates_dT(k_Si28_to_He4_Mg24) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Mg24_to_Si28);
    rate_eval.screened_rates(k_He4_Mg24_to_Si28) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mg24_to_Si28);
        rate_eval.dscreened_rates_dT(k_He4_Mg24_to_Si28) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 13.0_rt, 27.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_P31_to_He4_Al27);
    rate_eval.screened_rates(k_P31_to_He4_Al27) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_P31_to_He4_Al27);
        rate_eval.dscreened_rates_dT(k_P31_to_He4_Al27) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Al27_to_P31);
    rate_eval.screened_rates(k_He4_Al27_to_P31) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Al27_to_P31);
        rate_eval.dscreened_rates_dT(k_He4_Al27_to_P31) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 15.0_rt, 31.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_S32_to_p_P31);
    rate_eval.screened_rates(k_S32_to_p_P31) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_S32_to_p_P31);
        rate_eval.dscreened_rates_dT(k_S32_to_p_P31) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_P31_to_S32);
    rate_eval.screened_rates(k_p_P31_to_S32) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_P31_to_S32);
        rate_eval.dscreened_rates_dT(k_p_P31_to_S32) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Si28_to_p_P31);
    rate_eval.screened_rates(k_He4_Si28_to_p_P31) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Si28_to_p_P31);
        rate_eval.dscreened_rates_dT(k_He4_Si28_to_p_P31) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_P31_to_He4_Si28);
    rate_eval.screened_rates(k_p_P31_to_He4_Si28) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_P31_to_He4_Si28);
        rate_eval.dscreened_rates_dT(k_p_P31_to_He4_Si28) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 14.0_rt, 28.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_S32_to_He4_Si28);
    rate_eval.screened_rates(k_S32_to_He4_Si28) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_S32_to_He4_Si28);
        rate_eval.dscreened_rates_dT(k_S32_to_He4_Si28) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Si28_to_S32);
    rate_eval.screened_rates(k_He4_Si28_to_S32) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Si28_to_S32);
        rate_eval.dscreened_rates_dT(k_He4_Si28_to_S32) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 15.0_rt, 31.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_Cl35_to_He4_P31);
    rate_eval.screened_rates(k_Cl35_to_He4_P31) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_Cl35_to_He4_P31);
        rate_eval.dscreened_rates_dT(k_Cl35_to_He4_P31) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_P31_to_Cl35);
    rate_eval.screened_rates(k_He4_P31_to_Cl35) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_P31_to_Cl35);
        rate_eval.dscreened_rates_dT(k_He4_P31_to_Cl35) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 17.0_rt, 35.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_Ar36_to_p_Cl35);
    rate_eval.screened_rates(k_Ar36_to_p_Cl35) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_Ar36_to_p_Cl35);
        rate_eval.dscreened_rates_dT(k_Ar36_to_p_Cl35) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cl35_to_Ar36);
    rate_eval.screened_rates(k_p_Cl35_to_Ar36) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cl35_to_Ar36);
        rate_eval.dscreened_rates_dT(k_p_Cl35_to_Ar36) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_S32_to_p_Cl35);
    rate_eval.screened_rates(k_He4_S32_to_p_Cl35) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_S32_to_p_Cl35);
        rate_eval.dscreened_rates_dT(k_He4_S32_to_p_Cl35) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cl35_to_He4_S32);
    rate_eval.screened_rates(k_p_Cl35_to_He4_S32) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cl35_to_He4_S32);
        rate_eval.dscreened_rates_dT(k_p_Cl35_to_He4_S32) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 16.0_rt, 32.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_Ar36_to_He4_S32);
    rate_eval.screened_rates(k_Ar36_to_He4_S32) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_Ar36_to_He4_S32);
        rate_eval.dscreened_rates_dT(k_Ar36_to_He4_S32) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_S32_to_Ar36);
    rate_eval.screened_rates(k_He4_S32_to_Ar36) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_S32_to_Ar36);
        rate_eval.dscreened_rates_dT(k_He4_S32_to_Ar36) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 17.0_rt, 35.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_K39_to_He4_Cl35);
    rate_eval.screened_rates(k_K39_to_He4_Cl35) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_K39_to_He4_Cl35);
        rate_eval.dscreened_rates_dT(k_K39_to_He4_Cl35) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cl35_to_K39);
    rate_eval.screened_rates(k_He4_Cl35_to_K39) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cl35_to_K39);
        rate_eval.dscreened_rates_dT(k_He4_Cl35_to_K39) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 19.0_rt, 39.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_Ca40_to_p_K39);
    rate_eval.screened_rates(k_Ca40_to_p_K39) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_Ca40_to_p_K39);
        rate_eval.dscreened_rates_dT(k_Ca40_to_p_K39) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_K39_to_Ca40);
    rate_eval.screened_rates(k_p_K39_to_Ca40) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_K39_to_Ca40);
        rate_eval.dscreened_rates_dT(k_p_K39_to_Ca40) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ar36_to_p_K39);
    rate_eval.screened_rates(k_He4_Ar36_to_p_K39) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ar36_to_p_K39);
        rate_eval.dscreened_rates_dT(k_He4_Ar36_to_p_K39) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_K39_to_He4_Ar36);
    rate_eval.screened_rates(k_p_K39_to_He4_Ar36) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_K39_to_He4_Ar36);
        rate_eval.dscreened_rates_dT(k_p_K39_to_He4_Ar36) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 18.0_rt, 36.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_Ca40_to_He4_Ar36);
    rate_eval.screened_rates(k_Ca40_to_He4_Ar36) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_Ca40_to_He4_Ar36);
        rate_eval.dscreened_rates_dT(k_Ca40_to_He4_Ar36) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ar36_to_Ca40);
    rate_eval.screened_rates(k_He4_Ar36_to_Ca40) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ar36_to_Ca40);
        rate_eval.dscreened_rates_dT(k_He4_Ar36_to_Ca40) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 19.0_rt, 39.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_Sc43_to_He4_K39);
    rate_eval.screened_rates(k_Sc43_to_He4_K39) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_Sc43_to_He4_K39);
        rate_eval.dscreened_rates_dT(k_Sc43_to_He4_K39) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_K39_to_Sc43);
    rate_eval.screened_rates(k_He4_K39_to_Sc43) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_K39_to_Sc43);
        rate_eval.dscreened_rates_dT(k_He4_K39_to_Sc43) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 21.0_rt, 43.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_Ti44_to_p_Sc43);
    rate_eval.screened_rates(k_Ti44_to_p_Sc43) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_Ti44_to_p_Sc43);
        rate_eval.dscreened_rates_dT(k_Ti44_to_p_Sc43) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Sc43_to_Ti44);
    rate_eval.screened_rates(k_p_Sc43_to_Ti44) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Sc43_to_Ti44);
        rate_eval.dscreened_rates_dT(k_p_Sc43_to_Ti44) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ca40_to_p_Sc43);
    rate_eval.screened_rates(k_He4_Ca40_to_p_Sc43) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca40_to_p_Sc43);
        rate_eval.dscreened_rates_dT(k_He4_Ca40_to_p_Sc43) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Sc43_to_He4_Ca40);
    rate_eval.screened_rates(k_p_Sc43_to_He4_Ca40) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Sc43_to_He4_Ca40);
        rate_eval.dscreened_rates_dT(k_p_Sc43_to_He4_Ca40) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 20.0_rt, 40.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_Ti44_to_He4_Ca40);
    rate_eval.screened_rates(k_Ti44_to_He4_Ca40) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_Ti44_to_He4_Ca40);
        rate_eval.dscreened_rates_dT(k_Ti44_to_He4_Ca40) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ca40_to_Ti44);
    rate_eval.screened_rates(k_He4_Ca40_to_Ti44) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca40_to_Ti44);
        rate_eval.dscreened_rates_dT(k_He4_Ca40_to_Ti44) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 21.0_rt, 43.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_V47_to_He4_Sc43);
    rate_eval.screened_rates(k_V47_to_He4_Sc43) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_V47_to_He4_Sc43);
        rate_eval.dscreened_rates_dT(k_V47_to_He4_Sc43) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Sc43_to_V47);
    rate_eval.screened_rates(k_He4_Sc43_to_V47) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Sc43_to_V47);
        rate_eval.dscreened_rates_dT(k_He4_Sc43_to_V47) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 23.0_rt, 47.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_Cr48_to_p_V47);
    rate_eval.screened_rates(k_Cr48_to_p_V47) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_Cr48_to_p_V47);
        rate_eval.dscreened_rates_dT(k_Cr48_to_p_V47) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_V47_to_Cr48);
    rate_eval.screened_rates(k_p_V47_to_Cr48) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_V47_to_Cr48);
        rate_eval.dscreened_rates_dT(k_p_V47_to_Cr48) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ti44_to_p_V47);
    rate_eval.screened_rates(k_He4_Ti44_to_p_V47) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ti44_to_p_V47);
        rate_eval.dscreened_rates_dT(k_He4_Ti44_to_p_V47) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_V47_to_He4_Ti44);
    rate_eval.screened_rates(k_p_V47_to_He4_Ti44) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_V47_to_He4_Ti44);
        rate_eval.dscreened_rates_dT(k_p_V47_to_He4_Ti44) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 22.0_rt, 44.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_Cr48_to_He4_Ti44);
    rate_eval.screened_rates(k_Cr48_to_He4_Ti44) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_Cr48_to_He4_Ti44);
        rate_eval.dscreened_rates_dT(k_Cr48_to_He4_Ti44) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ti44_to_Cr48);
    rate_eval.screened_rates(k_He4_Ti44_to_Cr48) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ti44_to_Cr48);
        rate_eval.dscreened_rates_dT(k_He4_Ti44_to_Cr48) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 23.0_rt, 47.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_Mn51_to_He4_V47);
    rate_eval.screened_rates(k_Mn51_to_He4_V47) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_Mn51_to_He4_V47);
        rate_eval.dscreened_rates_dT(k_Mn51_to_He4_V47) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_V47_to_Mn51);
    rate_eval.screened_rates(k_He4_V47_to_Mn51) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_V47_to_Mn51);
        rate_eval.dscreened_rates_dT(k_He4_V47_to_Mn51) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 25.0_rt, 51.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_Fe52_to_p_Mn51);
    rate_eval.screened_rates(k_Fe52_to_p_Mn51) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_Fe52_to_p_Mn51);
        rate_eval.dscreened_rates_dT(k_Fe52_to_p_Mn51) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Mn51_to_Fe52);
    rate_eval.screened_rates(k_p_Mn51_to_Fe52) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mn51_to_Fe52);
        rate_eval.dscreened_rates_dT(k_p_Mn51_to_Fe52) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 24.0_rt, 48.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_Fe52_to_He4_Cr48);
    rate_eval.screened_rates(k_Fe52_to_He4_Cr48) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_Fe52_to_He4_Cr48);
        rate_eval.dscreened_rates_dT(k_Fe52_to_He4_Cr48) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cr48_to_Fe52);
    rate_eval.screened_rates(k_He4_Cr48_to_Fe52) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cr48_to_Fe52);
        rate_eval.dscreened_rates_dT(k_He4_Cr48_to_Fe52) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cr48_to_p_Mn51);
    rate_eval.screened_rates(k_He4_Cr48_to_p_Mn51) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cr48_to_p_Mn51);
        rate_eval.dscreened_rates_dT(k_He4_Cr48_to_p_Mn51) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Mn51_to_He4_Cr48);
    rate_eval.screened_rates(k_p_Mn51_to_He4_Cr48) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mn51_to_He4_Cr48);
        rate_eval.dscreened_rates_dT(k_p_Mn51_to_He4_Cr48) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 25.0_rt, 51.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_Co55_to_He4_Mn51);
    rate_eval.screened_rates(k_Co55_to_He4_Mn51) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_Co55_to_He4_Mn51);
        rate_eval.dscreened_rates_dT(k_Co55_to_He4_Mn51) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Mn51_to_Co55);
    rate_eval.screened_rates(k_He4_Mn51_to_Co55) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mn51_to_Co55);
        rate_eval.dscreened_rates_dT(k_He4_Mn51_to_Co55) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 27.0_rt, 55.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_Ni56_to_p_Co55);
    rate_eval.screened_rates(k_Ni56_to_p_Co55) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_Ni56_to_p_Co55);
        rate_eval.dscreened_rates_dT(k_Ni56_to_p_Co55) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Co55_to_Ni56);
    rate_eval.screened_rates(k_p_Co55_to_Ni56) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Co55_to_Ni56);
        rate_eval.dscreened_rates_dT(k_p_Co55_to_Ni56) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 26.0_rt, 52.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_Ni56_to_He4_Fe52);
    rate_eval.screened_rates(k_Ni56_to_He4_Fe52) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_Ni56_to_He4_Fe52);
        rate_eval.dscreened_rates_dT(k_Ni56_to_He4_Fe52) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Fe52_to_Ni56);
    rate_eval.screened_rates(k_He4_Fe52_to_Ni56) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe52_to_Ni56);
        rate_eval.dscreened_rates_dT(k_He4_Fe52_to_Ni56) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Fe52_to_p_Co55);
    rate_eval.screened_rates(k_He4_Fe52_to_p_Co55) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe52_to_p_Co55);
        rate_eval.dscreened_rates_dT(k_He4_Fe52_to_p_Co55) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Co55_to_He4_Fe52);
    rate_eval.screened_rates(k_p_Co55_to_He4_Fe52) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Co55_to_He4_Fe52);
        rate_eval.dscreened_rates_dT(k_p_Co55_to_He4_Fe52) = ratraw * dscor_dt + dratraw_dT * scor;
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


    ratraw = rate_eval.screened_rates(k_C12_to_He4_He4_He4);
    rate_eval.screened_rates(k_C12_to_He4_He4_He4) *= scor * scor2;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_to_He4_He4_He4);
        rate_eval.dscreened_rates_dT(k_C12_to_He4_He4_He4) = ratraw * (scor * dscor2_dt + dscor_dt * scor2) + dratraw_dT * scor * scor2;
    }

    ratraw = rate_eval.screened_rates(k_He4_He4_He4_to_C12);
    rate_eval.screened_rates(k_He4_He4_He4_to_C12) *= scor * scor2;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_He4_He4_to_C12);
        rate_eval.dscreened_rates_dT(k_He4_He4_He4_to_C12) = ratraw * (scor * dscor2_dt + dscor_dt * scor2) + dratraw_dT * scor * scor2;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(6.0_rt, 12.0_rt, 6.0_rt, 12.0_rt);


        static_assert(scn_fac.z1 == 6.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_C12_C12_to_p_Na23);
    rate_eval.screened_rates(k_C12_C12_to_p_Na23) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_C12_to_p_Na23);
        rate_eval.dscreened_rates_dT(k_C12_C12_to_p_Na23) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_C12_C12_to_He4_Ne20);
    rate_eval.screened_rates(k_C12_C12_to_He4_Ne20) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_C12_to_He4_Ne20);
        rate_eval.dscreened_rates_dT(k_C12_C12_to_He4_Ne20) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ne20_to_C12_C12);
    rate_eval.screened_rates(k_He4_Ne20_to_C12_C12) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ne20_to_C12_C12);
        rate_eval.dscreened_rates_dT(k_He4_Ne20_to_C12_C12) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Na23_to_C12_C12);
    rate_eval.screened_rates(k_p_Na23_to_C12_C12) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Na23_to_C12_C12);
        rate_eval.dscreened_rates_dT(k_p_Na23_to_C12_C12) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_C12_C12_to_Mg24_modified);
    rate_eval.screened_rates(k_C12_C12_to_Mg24_modified) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_C12_to_Mg24_modified);
        rate_eval.dscreened_rates_dT(k_C12_C12_to_Mg24_modified) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 7.0_rt, 13.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_N13_to_p_O16);
    rate_eval.screened_rates(k_He4_N13_to_p_O16) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_N13_to_p_O16);
        rate_eval.dscreened_rates_dT(k_He4_N13_to_p_O16) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_O16_to_He4_N13);
    rate_eval.screened_rates(k_p_O16_to_He4_N13) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_O16_to_He4_N13);
        rate_eval.dscreened_rates_dT(k_p_O16_to_He4_N13) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(6.0_rt, 12.0_rt, 8.0_rt, 16.0_rt);


        static_assert(scn_fac.z1 == 6.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_C12_O16_to_p_Al27);
    rate_eval.screened_rates(k_C12_O16_to_p_Al27) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_O16_to_p_Al27);
        rate_eval.dscreened_rates_dT(k_C12_O16_to_p_Al27) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_C12_O16_to_He4_Mg24);
    rate_eval.screened_rates(k_C12_O16_to_He4_Mg24) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_O16_to_He4_Mg24);
        rate_eval.dscreened_rates_dT(k_C12_O16_to_He4_Mg24) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Mg24_to_C12_O16);
    rate_eval.screened_rates(k_He4_Mg24_to_C12_O16) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mg24_to_C12_O16);
        rate_eval.dscreened_rates_dT(k_He4_Mg24_to_C12_O16) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Al27_to_C12_O16);
    rate_eval.screened_rates(k_p_Al27_to_C12_O16) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Al27_to_C12_O16);
        rate_eval.dscreened_rates_dT(k_p_Al27_to_C12_O16) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_C12_O16_to_Si28_modified);
    rate_eval.screened_rates(k_C12_O16_to_Si28_modified) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_O16_to_Si28_modified);
        rate_eval.dscreened_rates_dT(k_C12_O16_to_Si28_modified) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(8.0_rt, 16.0_rt, 8.0_rt, 16.0_rt);


        static_assert(scn_fac.z1 == 8.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_O16_O16_to_p_P31);
    rate_eval.screened_rates(k_O16_O16_to_p_P31) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_O16_O16_to_p_P31);
        rate_eval.dscreened_rates_dT(k_O16_O16_to_p_P31) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_O16_O16_to_He4_Si28);
    rate_eval.screened_rates(k_O16_O16_to_He4_Si28) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_O16_O16_to_He4_Si28);
        rate_eval.dscreened_rates_dT(k_O16_O16_to_He4_Si28) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Si28_to_O16_O16);
    rate_eval.screened_rates(k_He4_Si28_to_O16_O16) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Si28_to_O16_O16);
        rate_eval.dscreened_rates_dT(k_He4_Si28_to_O16_O16) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_P31_to_O16_O16);
    rate_eval.screened_rates(k_p_P31_to_O16_O16) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_P31_to_O16_O16);
        rate_eval.dscreened_rates_dT(k_p_P31_to_O16_O16) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_O16_O16_to_S32_modified);
    rate_eval.screened_rates(k_O16_O16_to_S32_modified) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_O16_O16_to_S32_modified);
        rate_eval.dscreened_rates_dT(k_O16_O16_to_S32_modified) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(6.0_rt, 12.0_rt, 10.0_rt, 20.0_rt);


        static_assert(scn_fac.z1 == 6.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_C12_Ne20_to_p_P31);
    rate_eval.screened_rates(k_C12_Ne20_to_p_P31) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_Ne20_to_p_P31);
        rate_eval.dscreened_rates_dT(k_C12_Ne20_to_p_P31) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_C12_Ne20_to_He4_Si28);
    rate_eval.screened_rates(k_C12_Ne20_to_He4_Si28) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_Ne20_to_He4_Si28);
        rate_eval.dscreened_rates_dT(k_C12_Ne20_to_He4_Si28) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Si28_to_C12_Ne20);
    rate_eval.screened_rates(k_He4_Si28_to_C12_Ne20) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Si28_to_C12_Ne20);
        rate_eval.dscreened_rates_dT(k_He4_Si28_to_C12_Ne20) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_P31_to_C12_Ne20);
    rate_eval.screened_rates(k_p_P31_to_C12_Ne20) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_P31_to_C12_Ne20);
        rate_eval.dscreened_rates_dT(k_p_P31_to_C12_Ne20) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    // Fill approximate rates

    fill_approx_rates<do_T_derivatives, T>(tfactors, rate_eval);

    // Calculate tabular rates

    [[maybe_unused]] Real rate, drate_dt, edot_nu, edot_gamma;

    rate_eval.enuc_weak = 0.0;


}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rhs_nuc(const burn_t& state,
             Array1D<Real, 1, neqs>& ydot_nuc,
             const Array1D<Real, 1, NumSpec>& Y,
             const Array1D<Real, 1, NumRates>& screened_rates) {

    using namespace Rates;

    ydot_nuc(H1) =
        (-screened_rates(k_p_C12_to_N13)*Y(C12)*Y(H1)*state.rho + screened_rates(k_N13_to_p_C12)*Y(N13)) +
        (-screened_rates(k_p_Ne21_to_Na22)*Y(Ne21)*Y(H1)*state.rho + screened_rates(k_Na22_to_p_Ne21)*Y(Na22)) +
        (-screened_rates(k_p_Na23_to_Mg24)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_Mg24_to_p_Na23)*Y(Mg24)) +
        (-screened_rates(k_p_Al27_to_Si28)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_Si28_to_p_Al27)*Y(Si28)) +
        (-screened_rates(k_p_P31_to_S32)*Y(P31)*Y(H1)*state.rho + screened_rates(k_S32_to_p_P31)*Y(S32)) +
        (-screened_rates(k_p_Cl35_to_Ar36)*Y(Cl35)*Y(H1)*state.rho + screened_rates(k_Ar36_to_p_Cl35)*Y(Ar36)) +
        (-screened_rates(k_p_K39_to_Ca40)*Y(K39)*Y(H1)*state.rho + screened_rates(k_Ca40_to_p_K39)*Y(Ca40)) +
        (-screened_rates(k_p_Sc43_to_Ti44)*Y(Sc43)*Y(H1)*state.rho + screened_rates(k_Ti44_to_p_Sc43)*Y(Ti44)) +
        (-screened_rates(k_p_V47_to_Cr48)*Y(V47)*Y(H1)*state.rho + screened_rates(k_Cr48_to_p_V47)*Y(Cr48)) +
        (-screened_rates(k_p_Mn51_to_Fe52)*Y(Mn51)*Y(H1)*state.rho + screened_rates(k_Fe52_to_p_Mn51)*Y(Fe52)) +
        (-screened_rates(k_p_Co55_to_Ni56)*Y(Co55)*Y(H1)*state.rho + screened_rates(k_Ni56_to_p_Co55)*Y(Ni56)) +
        (0.5*screened_rates(k_C12_C12_to_p_Na23)*std::pow(Y(C12), 2)*state.rho + -screened_rates(k_p_Na23_to_C12_C12)*Y(Na23)*Y(H1)*state.rho) +
        (screened_rates(k_He4_N13_to_p_O16)*Y(He4)*Y(N13)*state.rho + -screened_rates(k_p_O16_to_He4_N13)*Y(O16)*Y(H1)*state.rho) +
        (screened_rates(k_C12_O16_to_p_Al27)*Y(C12)*Y(O16)*state.rho + -screened_rates(k_p_Al27_to_C12_O16)*Y(Al27)*Y(H1)*state.rho) +
        (0.5*screened_rates(k_O16_O16_to_p_P31)*std::pow(Y(O16), 2)*state.rho + -screened_rates(k_p_P31_to_O16_O16)*Y(P31)*Y(H1)*state.rho) +
        (screened_rates(k_He4_F18_to_p_Ne21)*Y(F18)*Y(He4)*state.rho + -screened_rates(k_p_Ne21_to_He4_F18)*Y(Ne21)*Y(H1)*state.rho) +
        (screened_rates(k_C12_Ne20_to_p_P31)*Y(C12)*Y(Ne20)*state.rho + -screened_rates(k_p_P31_to_C12_Ne20)*Y(P31)*Y(H1)*state.rho) +
        (-screened_rates(k_p_Na23_to_He4_Ne20)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_He4_Ne20_to_p_Na23)*Y(He4)*Y(Ne20)*state.rho) +
        (-screened_rates(k_p_Al27_to_He4_Mg24)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_He4_Mg24_to_p_Al27)*Y(He4)*Y(Mg24)*state.rho) +
        (-screened_rates(k_p_P31_to_He4_Si28)*Y(P31)*Y(H1)*state.rho + screened_rates(k_He4_Si28_to_p_P31)*Y(He4)*Y(Si28)*state.rho) +
        (-screened_rates(k_p_Cl35_to_He4_S32)*Y(Cl35)*Y(H1)*state.rho + screened_rates(k_He4_S32_to_p_Cl35)*Y(He4)*Y(S32)*state.rho) +
        (-screened_rates(k_p_K39_to_He4_Ar36)*Y(K39)*Y(H1)*state.rho + screened_rates(k_He4_Ar36_to_p_K39)*Y(Ar36)*Y(He4)*state.rho) +
        (-screened_rates(k_p_Sc43_to_He4_Ca40)*Y(Sc43)*Y(H1)*state.rho + screened_rates(k_He4_Ca40_to_p_Sc43)*Y(Ca40)*Y(He4)*state.rho) +
        (-screened_rates(k_p_V47_to_He4_Ti44)*Y(V47)*Y(H1)*state.rho + screened_rates(k_He4_Ti44_to_p_V47)*Y(He4)*Y(Ti44)*state.rho) +
        (screened_rates(k_He4_Cr48_to_p_Mn51)*Y(Cr48)*Y(He4)*state.rho + -screened_rates(k_p_Mn51_to_He4_Cr48)*Y(Mn51)*Y(H1)*state.rho) +
        (screened_rates(k_He4_Fe52_to_p_Co55)*Y(Fe52)*Y(He4)*state.rho + -screened_rates(k_p_Co55_to_He4_Fe52)*Y(Co55)*Y(H1)*state.rho);

    ydot_nuc(He4) =
        (-screened_rates(k_He4_C12_to_O16)*Y(C12)*Y(He4)*state.rho + screened_rates(k_O16_to_He4_C12)*Y(O16)) +
        (-screened_rates(k_He4_N14_to_F18)*Y(He4)*Y(N14)*state.rho + screened_rates(k_F18_to_He4_N14)*Y(F18)) +
        (-screened_rates(k_He4_O16_to_Ne20)*Y(He4)*Y(O16)*state.rho + screened_rates(k_Ne20_to_He4_O16)*Y(Ne20)) +
        (-screened_rates(k_He4_F18_to_Na22)*Y(F18)*Y(He4)*state.rho + screened_rates(k_Na22_to_He4_F18)*Y(Na22)) +
        (-screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*Y(Ne20)*state.rho + screened_rates(k_Mg24_to_He4_Ne20)*Y(Mg24)) +
        (-screened_rates(k_He4_Na23_to_Al27)*Y(He4)*Y(Na23)*state.rho + screened_rates(k_Al27_to_He4_Na23)*Y(Al27)) +
        (-screened_rates(k_He4_Mg24_to_Si28)*Y(He4)*Y(Mg24)*state.rho + screened_rates(k_Si28_to_He4_Mg24)*Y(Si28)) +
        (-screened_rates(k_He4_Al27_to_P31)*Y(Al27)*Y(He4)*state.rho + screened_rates(k_P31_to_He4_Al27)*Y(P31)) +
        (-screened_rates(k_He4_Si28_to_S32)*Y(He4)*Y(Si28)*state.rho + screened_rates(k_S32_to_He4_Si28)*Y(S32)) +
        (-screened_rates(k_He4_P31_to_Cl35)*Y(He4)*Y(P31)*state.rho + screened_rates(k_Cl35_to_He4_P31)*Y(Cl35)) +
        (-screened_rates(k_He4_S32_to_Ar36)*Y(He4)*Y(S32)*state.rho + screened_rates(k_Ar36_to_He4_S32)*Y(Ar36)) +
        (-screened_rates(k_He4_Cl35_to_K39)*Y(Cl35)*Y(He4)*state.rho + screened_rates(k_K39_to_He4_Cl35)*Y(K39)) +
        (-screened_rates(k_He4_Ar36_to_Ca40)*Y(Ar36)*Y(He4)*state.rho + screened_rates(k_Ca40_to_He4_Ar36)*Y(Ca40)) +
        (-screened_rates(k_He4_K39_to_Sc43)*Y(He4)*Y(K39)*state.rho + screened_rates(k_Sc43_to_He4_K39)*Y(Sc43)) +
        (-screened_rates(k_He4_Ca40_to_Ti44)*Y(Ca40)*Y(He4)*state.rho + screened_rates(k_Ti44_to_He4_Ca40)*Y(Ti44)) +
        (-screened_rates(k_He4_Sc43_to_V47)*Y(He4)*Y(Sc43)*state.rho + screened_rates(k_V47_to_He4_Sc43)*Y(V47)) +
        (-screened_rates(k_He4_Ti44_to_Cr48)*Y(He4)*Y(Ti44)*state.rho + screened_rates(k_Cr48_to_He4_Ti44)*Y(Cr48)) +
        (-screened_rates(k_He4_V47_to_Mn51)*Y(He4)*Y(V47)*state.rho + screened_rates(k_Mn51_to_He4_V47)*Y(Mn51)) +
        (-screened_rates(k_He4_Cr48_to_Fe52)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_Fe52_to_He4_Cr48)*Y(Fe52)) +
        (-screened_rates(k_He4_Mn51_to_Co55)*Y(He4)*Y(Mn51)*state.rho + screened_rates(k_Co55_to_He4_Mn51)*Y(Co55)) +
        (-screened_rates(k_He4_Fe52_to_Ni56)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_Ni56_to_He4_Fe52)*Y(Ni56)) +
        (0.5*screened_rates(k_C12_C12_to_He4_Ne20)*std::pow(Y(C12), 2)*state.rho + -screened_rates(k_He4_Ne20_to_C12_C12)*Y(He4)*Y(Ne20)*state.rho) +
        (-screened_rates(k_He4_N13_to_p_O16)*Y(He4)*Y(N13)*state.rho + screened_rates(k_p_O16_to_He4_N13)*Y(O16)*Y(H1)*state.rho) +
        (screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*Y(O16)*state.rho + -screened_rates(k_He4_Mg24_to_C12_O16)*Y(He4)*Y(Mg24)*state.rho) +
        (0.5*screened_rates(k_O16_O16_to_He4_Si28)*std::pow(Y(O16), 2)*state.rho + -screened_rates(k_He4_Si28_to_O16_O16)*Y(He4)*Y(Si28)*state.rho) +
        (-screened_rates(k_He4_F18_to_p_Ne21)*Y(F18)*Y(He4)*state.rho + screened_rates(k_p_Ne21_to_He4_F18)*Y(Ne21)*Y(H1)*state.rho) +
        (screened_rates(k_C12_Ne20_to_He4_Si28)*Y(C12)*Y(Ne20)*state.rho + -screened_rates(k_He4_Si28_to_C12_Ne20)*Y(He4)*Y(Si28)*state.rho) +
        (screened_rates(k_p_Na23_to_He4_Ne20)*Y(Na23)*Y(H1)*state.rho + -screened_rates(k_He4_Ne20_to_p_Na23)*Y(He4)*Y(Ne20)*state.rho) +
        (screened_rates(k_p_Al27_to_He4_Mg24)*Y(Al27)*Y(H1)*state.rho + -screened_rates(k_He4_Mg24_to_p_Al27)*Y(He4)*Y(Mg24)*state.rho) +
        (screened_rates(k_p_P31_to_He4_Si28)*Y(P31)*Y(H1)*state.rho + -screened_rates(k_He4_Si28_to_p_P31)*Y(He4)*Y(Si28)*state.rho) +
        (screened_rates(k_p_Cl35_to_He4_S32)*Y(Cl35)*Y(H1)*state.rho + -screened_rates(k_He4_S32_to_p_Cl35)*Y(He4)*Y(S32)*state.rho) +
        (screened_rates(k_p_K39_to_He4_Ar36)*Y(K39)*Y(H1)*state.rho + -screened_rates(k_He4_Ar36_to_p_K39)*Y(Ar36)*Y(He4)*state.rho) +
        (screened_rates(k_p_Sc43_to_He4_Ca40)*Y(Sc43)*Y(H1)*state.rho + -screened_rates(k_He4_Ca40_to_p_Sc43)*Y(Ca40)*Y(He4)*state.rho) +
        (screened_rates(k_p_V47_to_He4_Ti44)*Y(V47)*Y(H1)*state.rho + -screened_rates(k_He4_Ti44_to_p_V47)*Y(He4)*Y(Ti44)*state.rho) +
        (-screened_rates(k_He4_Cr48_to_p_Mn51)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_p_Mn51_to_He4_Cr48)*Y(Mn51)*Y(H1)*state.rho) +
        (-screened_rates(k_He4_Fe52_to_p_Co55)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_p_Co55_to_He4_Fe52)*Y(Co55)*Y(H1)*state.rho) +
        (-0.5*screened_rates(k_He4_He4_He4_to_C12)*std::pow(Y(He4), 3)*std::pow(state.rho, 2) + 3.0*screened_rates(k_C12_to_He4_He4_He4)*Y(C12));

    ydot_nuc(C12) =
        (-screened_rates(k_p_C12_to_N13)*Y(C12)*Y(H1)*state.rho + screened_rates(k_N13_to_p_C12)*Y(N13)) +
        (-screened_rates(k_He4_C12_to_O16)*Y(C12)*Y(He4)*state.rho + screened_rates(k_O16_to_He4_C12)*Y(O16)) +
        (-screened_rates(k_C12_C12_to_p_Na23)*std::pow(Y(C12), 2)*state.rho + 2.0*screened_rates(k_p_Na23_to_C12_C12)*Y(Na23)*Y(H1)*state.rho) +
        (-screened_rates(k_C12_C12_to_He4_Ne20)*std::pow(Y(C12), 2)*state.rho + 2.0*screened_rates(k_He4_Ne20_to_C12_C12)*Y(He4)*Y(Ne20)*state.rho) +
        (-screened_rates(k_C12_O16_to_p_Al27)*Y(C12)*Y(O16)*state.rho + screened_rates(k_p_Al27_to_C12_O16)*Y(Al27)*Y(H1)*state.rho) +
        (-screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*Y(O16)*state.rho + screened_rates(k_He4_Mg24_to_C12_O16)*Y(He4)*Y(Mg24)*state.rho) +
        (-screened_rates(k_C12_Ne20_to_p_P31)*Y(C12)*Y(Ne20)*state.rho + screened_rates(k_p_P31_to_C12_Ne20)*Y(P31)*Y(H1)*state.rho) +
        (-screened_rates(k_C12_Ne20_to_He4_Si28)*Y(C12)*Y(Ne20)*state.rho + screened_rates(k_He4_Si28_to_C12_Ne20)*Y(He4)*Y(Si28)*state.rho) +
        (0.16666666666666667*screened_rates(k_He4_He4_He4_to_C12)*std::pow(Y(He4), 3)*std::pow(state.rho, 2) + -screened_rates(k_C12_to_He4_He4_He4)*Y(C12)) +
        -screened_rates(k_C12_C12_to_Mg24_modified)*std::pow(Y(C12), 2)*state.rho +
        -screened_rates(k_C12_O16_to_Si28_modified)*Y(C12)*Y(O16)*state.rho;

    ydot_nuc(N13) =
        (screened_rates(k_p_C12_to_N13)*Y(C12)*Y(H1)*state.rho + -screened_rates(k_N13_to_p_C12)*Y(N13)) +
        (-screened_rates(k_He4_N13_to_p_O16)*Y(He4)*Y(N13)*state.rho + screened_rates(k_p_O16_to_He4_N13)*Y(O16)*Y(H1)*state.rho);

    ydot_nuc(N14) =
        (-screened_rates(k_He4_N14_to_F18)*Y(He4)*Y(N14)*state.rho + screened_rates(k_F18_to_He4_N14)*Y(F18));

    ydot_nuc(O16) =
        (screened_rates(k_He4_C12_to_O16)*Y(C12)*Y(He4)*state.rho + -screened_rates(k_O16_to_He4_C12)*Y(O16)) +
        (-screened_rates(k_He4_O16_to_Ne20)*Y(He4)*Y(O16)*state.rho + screened_rates(k_Ne20_to_He4_O16)*Y(Ne20)) +
        (screened_rates(k_He4_N13_to_p_O16)*Y(He4)*Y(N13)*state.rho + -screened_rates(k_p_O16_to_He4_N13)*Y(O16)*Y(H1)*state.rho) +
        (-screened_rates(k_C12_O16_to_p_Al27)*Y(C12)*Y(O16)*state.rho + screened_rates(k_p_Al27_to_C12_O16)*Y(Al27)*Y(H1)*state.rho) +
        (-screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*Y(O16)*state.rho + screened_rates(k_He4_Mg24_to_C12_O16)*Y(He4)*Y(Mg24)*state.rho) +
        (-screened_rates(k_O16_O16_to_p_P31)*std::pow(Y(O16), 2)*state.rho + 2.0*screened_rates(k_p_P31_to_O16_O16)*Y(P31)*Y(H1)*state.rho) +
        (-screened_rates(k_O16_O16_to_He4_Si28)*std::pow(Y(O16), 2)*state.rho + 2.0*screened_rates(k_He4_Si28_to_O16_O16)*Y(He4)*Y(Si28)*state.rho) +
        -screened_rates(k_O16_O16_to_S32_modified)*std::pow(Y(O16), 2)*state.rho +
        -screened_rates(k_C12_O16_to_Si28_modified)*Y(C12)*Y(O16)*state.rho;

    ydot_nuc(F18) =
        (screened_rates(k_He4_N14_to_F18)*Y(He4)*Y(N14)*state.rho + -screened_rates(k_F18_to_He4_N14)*Y(F18)) +
        (-screened_rates(k_He4_F18_to_Na22)*Y(F18)*Y(He4)*state.rho + screened_rates(k_Na22_to_He4_F18)*Y(Na22)) +
        (-screened_rates(k_He4_F18_to_p_Ne21)*Y(F18)*Y(He4)*state.rho + screened_rates(k_p_Ne21_to_He4_F18)*Y(Ne21)*Y(H1)*state.rho);

    ydot_nuc(Ne20) =
        (screened_rates(k_He4_O16_to_Ne20)*Y(He4)*Y(O16)*state.rho + -screened_rates(k_Ne20_to_He4_O16)*Y(Ne20)) +
        (-screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*Y(Ne20)*state.rho + screened_rates(k_Mg24_to_He4_Ne20)*Y(Mg24)) +
        (0.5*screened_rates(k_C12_C12_to_He4_Ne20)*std::pow(Y(C12), 2)*state.rho + -screened_rates(k_He4_Ne20_to_C12_C12)*Y(He4)*Y(Ne20)*state.rho) +
        (-screened_rates(k_C12_Ne20_to_p_P31)*Y(C12)*Y(Ne20)*state.rho + screened_rates(k_p_P31_to_C12_Ne20)*Y(P31)*Y(H1)*state.rho) +
        (-screened_rates(k_C12_Ne20_to_He4_Si28)*Y(C12)*Y(Ne20)*state.rho + screened_rates(k_He4_Si28_to_C12_Ne20)*Y(He4)*Y(Si28)*state.rho) +
        (screened_rates(k_p_Na23_to_He4_Ne20)*Y(Na23)*Y(H1)*state.rho + -screened_rates(k_He4_Ne20_to_p_Na23)*Y(He4)*Y(Ne20)*state.rho);

    ydot_nuc(Ne21) =
        (-screened_rates(k_p_Ne21_to_Na22)*Y(Ne21)*Y(H1)*state.rho + screened_rates(k_Na22_to_p_Ne21)*Y(Na22)) +
        (screened_rates(k_He4_F18_to_p_Ne21)*Y(F18)*Y(He4)*state.rho + -screened_rates(k_p_Ne21_to_He4_F18)*Y(Ne21)*Y(H1)*state.rho);

    ydot_nuc(Na22) =
        (screened_rates(k_He4_F18_to_Na22)*Y(F18)*Y(He4)*state.rho + -screened_rates(k_Na22_to_He4_F18)*Y(Na22)) +
        (screened_rates(k_p_Ne21_to_Na22)*Y(Ne21)*Y(H1)*state.rho + -screened_rates(k_Na22_to_p_Ne21)*Y(Na22));

    ydot_nuc(Na23) =
        (-screened_rates(k_p_Na23_to_Mg24)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_Mg24_to_p_Na23)*Y(Mg24)) +
        (-screened_rates(k_He4_Na23_to_Al27)*Y(He4)*Y(Na23)*state.rho + screened_rates(k_Al27_to_He4_Na23)*Y(Al27)) +
        (0.5*screened_rates(k_C12_C12_to_p_Na23)*std::pow(Y(C12), 2)*state.rho + -screened_rates(k_p_Na23_to_C12_C12)*Y(Na23)*Y(H1)*state.rho) +
        (-screened_rates(k_p_Na23_to_He4_Ne20)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_He4_Ne20_to_p_Na23)*Y(He4)*Y(Ne20)*state.rho);

    ydot_nuc(Mg24) =
        (screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*Y(Ne20)*state.rho + -screened_rates(k_Mg24_to_He4_Ne20)*Y(Mg24)) +
        (screened_rates(k_p_Na23_to_Mg24)*Y(Na23)*Y(H1)*state.rho + -screened_rates(k_Mg24_to_p_Na23)*Y(Mg24)) +
        (-screened_rates(k_He4_Mg24_to_Si28)*Y(He4)*Y(Mg24)*state.rho + screened_rates(k_Si28_to_He4_Mg24)*Y(Si28)) +
        (screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*Y(O16)*state.rho + -screened_rates(k_He4_Mg24_to_C12_O16)*Y(He4)*Y(Mg24)*state.rho) +
        (screened_rates(k_p_Al27_to_He4_Mg24)*Y(Al27)*Y(H1)*state.rho + -screened_rates(k_He4_Mg24_to_p_Al27)*Y(He4)*Y(Mg24)*state.rho) +
        0.5*screened_rates(k_C12_C12_to_Mg24_modified)*std::pow(Y(C12), 2)*state.rho;

    ydot_nuc(Al27) =
        (screened_rates(k_He4_Na23_to_Al27)*Y(He4)*Y(Na23)*state.rho + -screened_rates(k_Al27_to_He4_Na23)*Y(Al27)) +
        (-screened_rates(k_p_Al27_to_Si28)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_Si28_to_p_Al27)*Y(Si28)) +
        (-screened_rates(k_He4_Al27_to_P31)*Y(Al27)*Y(He4)*state.rho + screened_rates(k_P31_to_He4_Al27)*Y(P31)) +
        (screened_rates(k_C12_O16_to_p_Al27)*Y(C12)*Y(O16)*state.rho + -screened_rates(k_p_Al27_to_C12_O16)*Y(Al27)*Y(H1)*state.rho) +
        (-screened_rates(k_p_Al27_to_He4_Mg24)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_He4_Mg24_to_p_Al27)*Y(He4)*Y(Mg24)*state.rho);

    ydot_nuc(Si28) =
        (screened_rates(k_He4_Mg24_to_Si28)*Y(He4)*Y(Mg24)*state.rho + -screened_rates(k_Si28_to_He4_Mg24)*Y(Si28)) +
        (screened_rates(k_p_Al27_to_Si28)*Y(Al27)*Y(H1)*state.rho + -screened_rates(k_Si28_to_p_Al27)*Y(Si28)) +
        (-screened_rates(k_He4_Si28_to_S32)*Y(He4)*Y(Si28)*state.rho + screened_rates(k_S32_to_He4_Si28)*Y(S32)) +
        (0.5*screened_rates(k_O16_O16_to_He4_Si28)*std::pow(Y(O16), 2)*state.rho + -screened_rates(k_He4_Si28_to_O16_O16)*Y(He4)*Y(Si28)*state.rho) +
        (screened_rates(k_C12_Ne20_to_He4_Si28)*Y(C12)*Y(Ne20)*state.rho + -screened_rates(k_He4_Si28_to_C12_Ne20)*Y(He4)*Y(Si28)*state.rho) +
        (screened_rates(k_p_P31_to_He4_Si28)*Y(P31)*Y(H1)*state.rho + -screened_rates(k_He4_Si28_to_p_P31)*Y(He4)*Y(Si28)*state.rho) +
        screened_rates(k_C12_O16_to_Si28_modified)*Y(C12)*Y(O16)*state.rho;

    ydot_nuc(P31) =
        (screened_rates(k_He4_Al27_to_P31)*Y(Al27)*Y(He4)*state.rho + -screened_rates(k_P31_to_He4_Al27)*Y(P31)) +
        (-screened_rates(k_p_P31_to_S32)*Y(P31)*Y(H1)*state.rho + screened_rates(k_S32_to_p_P31)*Y(S32)) +
        (-screened_rates(k_He4_P31_to_Cl35)*Y(He4)*Y(P31)*state.rho + screened_rates(k_Cl35_to_He4_P31)*Y(Cl35)) +
        (0.5*screened_rates(k_O16_O16_to_p_P31)*std::pow(Y(O16), 2)*state.rho + -screened_rates(k_p_P31_to_O16_O16)*Y(P31)*Y(H1)*state.rho) +
        (screened_rates(k_C12_Ne20_to_p_P31)*Y(C12)*Y(Ne20)*state.rho + -screened_rates(k_p_P31_to_C12_Ne20)*Y(P31)*Y(H1)*state.rho) +
        (-screened_rates(k_p_P31_to_He4_Si28)*Y(P31)*Y(H1)*state.rho + screened_rates(k_He4_Si28_to_p_P31)*Y(He4)*Y(Si28)*state.rho);

    ydot_nuc(S32) =
        (screened_rates(k_He4_Si28_to_S32)*Y(He4)*Y(Si28)*state.rho + -screened_rates(k_S32_to_He4_Si28)*Y(S32)) +
        (screened_rates(k_p_P31_to_S32)*Y(P31)*Y(H1)*state.rho + -screened_rates(k_S32_to_p_P31)*Y(S32)) +
        (-screened_rates(k_He4_S32_to_Ar36)*Y(He4)*Y(S32)*state.rho + screened_rates(k_Ar36_to_He4_S32)*Y(Ar36)) +
        (screened_rates(k_p_Cl35_to_He4_S32)*Y(Cl35)*Y(H1)*state.rho + -screened_rates(k_He4_S32_to_p_Cl35)*Y(He4)*Y(S32)*state.rho) +
        0.5*screened_rates(k_O16_O16_to_S32_modified)*std::pow(Y(O16), 2)*state.rho;

    ydot_nuc(Cl35) =
        (screened_rates(k_He4_P31_to_Cl35)*Y(He4)*Y(P31)*state.rho + -screened_rates(k_Cl35_to_He4_P31)*Y(Cl35)) +
        (-screened_rates(k_p_Cl35_to_Ar36)*Y(Cl35)*Y(H1)*state.rho + screened_rates(k_Ar36_to_p_Cl35)*Y(Ar36)) +
        (-screened_rates(k_He4_Cl35_to_K39)*Y(Cl35)*Y(He4)*state.rho + screened_rates(k_K39_to_He4_Cl35)*Y(K39)) +
        (-screened_rates(k_p_Cl35_to_He4_S32)*Y(Cl35)*Y(H1)*state.rho + screened_rates(k_He4_S32_to_p_Cl35)*Y(He4)*Y(S32)*state.rho);

    ydot_nuc(Ar36) =
        (screened_rates(k_He4_S32_to_Ar36)*Y(He4)*Y(S32)*state.rho + -screened_rates(k_Ar36_to_He4_S32)*Y(Ar36)) +
        (screened_rates(k_p_Cl35_to_Ar36)*Y(Cl35)*Y(H1)*state.rho + -screened_rates(k_Ar36_to_p_Cl35)*Y(Ar36)) +
        (-screened_rates(k_He4_Ar36_to_Ca40)*Y(Ar36)*Y(He4)*state.rho + screened_rates(k_Ca40_to_He4_Ar36)*Y(Ca40)) +
        (screened_rates(k_p_K39_to_He4_Ar36)*Y(K39)*Y(H1)*state.rho + -screened_rates(k_He4_Ar36_to_p_K39)*Y(Ar36)*Y(He4)*state.rho);

    ydot_nuc(K39) =
        (screened_rates(k_He4_Cl35_to_K39)*Y(Cl35)*Y(He4)*state.rho + -screened_rates(k_K39_to_He4_Cl35)*Y(K39)) +
        (-screened_rates(k_p_K39_to_Ca40)*Y(K39)*Y(H1)*state.rho + screened_rates(k_Ca40_to_p_K39)*Y(Ca40)) +
        (-screened_rates(k_He4_K39_to_Sc43)*Y(He4)*Y(K39)*state.rho + screened_rates(k_Sc43_to_He4_K39)*Y(Sc43)) +
        (-screened_rates(k_p_K39_to_He4_Ar36)*Y(K39)*Y(H1)*state.rho + screened_rates(k_He4_Ar36_to_p_K39)*Y(Ar36)*Y(He4)*state.rho);

    ydot_nuc(Ca40) =
        (screened_rates(k_He4_Ar36_to_Ca40)*Y(Ar36)*Y(He4)*state.rho + -screened_rates(k_Ca40_to_He4_Ar36)*Y(Ca40)) +
        (screened_rates(k_p_K39_to_Ca40)*Y(K39)*Y(H1)*state.rho + -screened_rates(k_Ca40_to_p_K39)*Y(Ca40)) +
        (-screened_rates(k_He4_Ca40_to_Ti44)*Y(Ca40)*Y(He4)*state.rho + screened_rates(k_Ti44_to_He4_Ca40)*Y(Ti44)) +
        (screened_rates(k_p_Sc43_to_He4_Ca40)*Y(Sc43)*Y(H1)*state.rho + -screened_rates(k_He4_Ca40_to_p_Sc43)*Y(Ca40)*Y(He4)*state.rho);

    ydot_nuc(Sc43) =
        (screened_rates(k_He4_K39_to_Sc43)*Y(He4)*Y(K39)*state.rho + -screened_rates(k_Sc43_to_He4_K39)*Y(Sc43)) +
        (-screened_rates(k_p_Sc43_to_Ti44)*Y(Sc43)*Y(H1)*state.rho + screened_rates(k_Ti44_to_p_Sc43)*Y(Ti44)) +
        (-screened_rates(k_He4_Sc43_to_V47)*Y(He4)*Y(Sc43)*state.rho + screened_rates(k_V47_to_He4_Sc43)*Y(V47)) +
        (-screened_rates(k_p_Sc43_to_He4_Ca40)*Y(Sc43)*Y(H1)*state.rho + screened_rates(k_He4_Ca40_to_p_Sc43)*Y(Ca40)*Y(He4)*state.rho);

    ydot_nuc(Ti44) =
        (screened_rates(k_He4_Ca40_to_Ti44)*Y(Ca40)*Y(He4)*state.rho + -screened_rates(k_Ti44_to_He4_Ca40)*Y(Ti44)) +
        (screened_rates(k_p_Sc43_to_Ti44)*Y(Sc43)*Y(H1)*state.rho + -screened_rates(k_Ti44_to_p_Sc43)*Y(Ti44)) +
        (-screened_rates(k_He4_Ti44_to_Cr48)*Y(He4)*Y(Ti44)*state.rho + screened_rates(k_Cr48_to_He4_Ti44)*Y(Cr48)) +
        (screened_rates(k_p_V47_to_He4_Ti44)*Y(V47)*Y(H1)*state.rho + -screened_rates(k_He4_Ti44_to_p_V47)*Y(He4)*Y(Ti44)*state.rho);

    ydot_nuc(V47) =
        (screened_rates(k_He4_Sc43_to_V47)*Y(He4)*Y(Sc43)*state.rho + -screened_rates(k_V47_to_He4_Sc43)*Y(V47)) +
        (-screened_rates(k_p_V47_to_Cr48)*Y(V47)*Y(H1)*state.rho + screened_rates(k_Cr48_to_p_V47)*Y(Cr48)) +
        (-screened_rates(k_He4_V47_to_Mn51)*Y(He4)*Y(V47)*state.rho + screened_rates(k_Mn51_to_He4_V47)*Y(Mn51)) +
        (-screened_rates(k_p_V47_to_He4_Ti44)*Y(V47)*Y(H1)*state.rho + screened_rates(k_He4_Ti44_to_p_V47)*Y(He4)*Y(Ti44)*state.rho);

    ydot_nuc(Cr48) =
        (screened_rates(k_He4_Ti44_to_Cr48)*Y(He4)*Y(Ti44)*state.rho + -screened_rates(k_Cr48_to_He4_Ti44)*Y(Cr48)) +
        (screened_rates(k_p_V47_to_Cr48)*Y(V47)*Y(H1)*state.rho + -screened_rates(k_Cr48_to_p_V47)*Y(Cr48)) +
        (-screened_rates(k_He4_Cr48_to_Fe52)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_Fe52_to_He4_Cr48)*Y(Fe52)) +
        (-screened_rates(k_He4_Cr48_to_p_Mn51)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_p_Mn51_to_He4_Cr48)*Y(Mn51)*Y(H1)*state.rho);

    ydot_nuc(Mn51) =
        (screened_rates(k_He4_V47_to_Mn51)*Y(He4)*Y(V47)*state.rho + -screened_rates(k_Mn51_to_He4_V47)*Y(Mn51)) +
        (-screened_rates(k_p_Mn51_to_Fe52)*Y(Mn51)*Y(H1)*state.rho + screened_rates(k_Fe52_to_p_Mn51)*Y(Fe52)) +
        (-screened_rates(k_He4_Mn51_to_Co55)*Y(He4)*Y(Mn51)*state.rho + screened_rates(k_Co55_to_He4_Mn51)*Y(Co55)) +
        (screened_rates(k_He4_Cr48_to_p_Mn51)*Y(Cr48)*Y(He4)*state.rho + -screened_rates(k_p_Mn51_to_He4_Cr48)*Y(Mn51)*Y(H1)*state.rho);

    ydot_nuc(Fe52) =
        (screened_rates(k_He4_Cr48_to_Fe52)*Y(Cr48)*Y(He4)*state.rho + -screened_rates(k_Fe52_to_He4_Cr48)*Y(Fe52)) +
        (screened_rates(k_p_Mn51_to_Fe52)*Y(Mn51)*Y(H1)*state.rho + -screened_rates(k_Fe52_to_p_Mn51)*Y(Fe52)) +
        (-screened_rates(k_He4_Fe52_to_Ni56)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_Ni56_to_He4_Fe52)*Y(Ni56)) +
        (-screened_rates(k_He4_Fe52_to_p_Co55)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_p_Co55_to_He4_Fe52)*Y(Co55)*Y(H1)*state.rho);

    ydot_nuc(Co55) =
        (screened_rates(k_He4_Mn51_to_Co55)*Y(He4)*Y(Mn51)*state.rho + -screened_rates(k_Co55_to_He4_Mn51)*Y(Co55)) +
        (-screened_rates(k_p_Co55_to_Ni56)*Y(Co55)*Y(H1)*state.rho + screened_rates(k_Ni56_to_p_Co55)*Y(Ni56)) +
        (screened_rates(k_He4_Fe52_to_p_Co55)*Y(Fe52)*Y(He4)*state.rho + -screened_rates(k_p_Co55_to_He4_Fe52)*Y(Co55)*Y(H1)*state.rho);

    ydot_nuc(Ni56) =
        (screened_rates(k_He4_Fe52_to_Ni56)*Y(Fe52)*Y(He4)*state.rho + -screened_rates(k_Ni56_to_He4_Fe52)*Y(Ni56)) +
        (screened_rates(k_p_Co55_to_Ni56)*Y(Co55)*Y(H1)*state.rho + -screened_rates(k_Ni56_to_p_Co55)*Y(Ni56));

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

    // include any weak rate neutrino losses
    enuc += rate_eval.enuc_weak;

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

    scratch = -screened_rates(k_p_Al27_to_C12_O16)*Y(Al27)*state.rho - screened_rates(k_p_Al27_to_He4_Mg24)*Y(Al27)*state.rho - screened_rates(k_p_Al27_to_Si28)*Y(Al27)*state.rho - screened_rates(k_p_C12_to_N13)*Y(C12)*state.rho - screened_rates(k_p_Cl35_to_Ar36)*Y(Cl35)*state.rho - screened_rates(k_p_Cl35_to_He4_S32)*Y(Cl35)*state.rho - screened_rates(k_p_Co55_to_He4_Fe52)*Y(Co55)*state.rho - screened_rates(k_p_Co55_to_Ni56)*Y(Co55)*state.rho - screened_rates(k_p_K39_to_Ca40)*Y(K39)*state.rho - screened_rates(k_p_K39_to_He4_Ar36)*Y(K39)*state.rho - screened_rates(k_p_Mn51_to_Fe52)*Y(Mn51)*state.rho - screened_rates(k_p_Mn51_to_He4_Cr48)*Y(Mn51)*state.rho - screened_rates(k_p_Na23_to_C12_C12)*Y(Na23)*state.rho - screened_rates(k_p_Na23_to_He4_Ne20)*Y(Na23)*state.rho - screened_rates(k_p_Na23_to_Mg24)*Y(Na23)*state.rho - screened_rates(k_p_Ne21_to_He4_F18)*Y(Ne21)*state.rho - screened_rates(k_p_Ne21_to_Na22)*Y(Ne21)*state.rho - screened_rates(k_p_O16_to_He4_N13)*Y(O16)*state.rho - screened_rates(k_p_P31_to_C12_Ne20)*Y(P31)*state.rho - screened_rates(k_p_P31_to_He4_Si28)*Y(P31)*state.rho - screened_rates(k_p_P31_to_O16_O16)*Y(P31)*state.rho - screened_rates(k_p_P31_to_S32)*Y(P31)*state.rho - screened_rates(k_p_Sc43_to_He4_Ca40)*Y(Sc43)*state.rho - screened_rates(k_p_Sc43_to_Ti44)*Y(Sc43)*state.rho - screened_rates(k_p_V47_to_Cr48)*Y(V47)*state.rho - screened_rates(k_p_V47_to_He4_Ti44)*Y(V47)*state.rho;
    jac.set(H1, H1, scratch);

    scratch = screened_rates(k_He4_Ar36_to_p_K39)*Y(Ar36)*state.rho + screened_rates(k_He4_Ca40_to_p_Sc43)*Y(Ca40)*state.rho + screened_rates(k_He4_Cr48_to_p_Mn51)*Y(Cr48)*state.rho + screened_rates(k_He4_F18_to_p_Ne21)*Y(F18)*state.rho + screened_rates(k_He4_Fe52_to_p_Co55)*Y(Fe52)*state.rho + screened_rates(k_He4_Mg24_to_p_Al27)*Y(Mg24)*state.rho + screened_rates(k_He4_N13_to_p_O16)*Y(N13)*state.rho + screened_rates(k_He4_Ne20_to_p_Na23)*Y(Ne20)*state.rho + screened_rates(k_He4_S32_to_p_Cl35)*Y(S32)*state.rho + screened_rates(k_He4_Si28_to_p_P31)*Y(Si28)*state.rho + screened_rates(k_He4_Ti44_to_p_V47)*Y(Ti44)*state.rho;
    jac.set(H1, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_p_Na23)*Y(C12)*state.rho + screened_rates(k_C12_Ne20_to_p_P31)*Y(Ne20)*state.rho + screened_rates(k_C12_O16_to_p_Al27)*Y(O16)*state.rho - screened_rates(k_p_C12_to_N13)*Y(H1)*state.rho;
    jac.set(H1, C12, scratch);

    scratch = screened_rates(k_He4_N13_to_p_O16)*Y(He4)*state.rho + screened_rates(k_N13_to_p_C12);
    jac.set(H1, N13, scratch);

    scratch = screened_rates(k_C12_O16_to_p_Al27)*Y(C12)*state.rho + 1.0*screened_rates(k_O16_O16_to_p_P31)*Y(O16)*state.rho - screened_rates(k_p_O16_to_He4_N13)*Y(H1)*state.rho;
    jac.set(H1, O16, scratch);

    scratch = screened_rates(k_He4_F18_to_p_Ne21)*Y(He4)*state.rho;
    jac.set(H1, F18, scratch);

    scratch = screened_rates(k_C12_Ne20_to_p_P31)*Y(C12)*state.rho + screened_rates(k_He4_Ne20_to_p_Na23)*Y(He4)*state.rho;
    jac.set(H1, Ne20, scratch);

    scratch = -screened_rates(k_p_Ne21_to_He4_F18)*Y(H1)*state.rho - screened_rates(k_p_Ne21_to_Na22)*Y(H1)*state.rho;
    jac.set(H1, Ne21, scratch);

    scratch = screened_rates(k_Na22_to_p_Ne21);
    jac.set(H1, Na22, scratch);

    scratch = -screened_rates(k_p_Na23_to_C12_C12)*Y(H1)*state.rho - screened_rates(k_p_Na23_to_He4_Ne20)*Y(H1)*state.rho - screened_rates(k_p_Na23_to_Mg24)*Y(H1)*state.rho;
    jac.set(H1, Na23, scratch);

    scratch = screened_rates(k_He4_Mg24_to_p_Al27)*Y(He4)*state.rho + screened_rates(k_Mg24_to_p_Na23);
    jac.set(H1, Mg24, scratch);

    scratch = -screened_rates(k_p_Al27_to_C12_O16)*Y(H1)*state.rho - screened_rates(k_p_Al27_to_He4_Mg24)*Y(H1)*state.rho - screened_rates(k_p_Al27_to_Si28)*Y(H1)*state.rho;
    jac.set(H1, Al27, scratch);

    scratch = screened_rates(k_He4_Si28_to_p_P31)*Y(He4)*state.rho + screened_rates(k_Si28_to_p_Al27);
    jac.set(H1, Si28, scratch);

    scratch = -screened_rates(k_p_P31_to_C12_Ne20)*Y(H1)*state.rho - screened_rates(k_p_P31_to_He4_Si28)*Y(H1)*state.rho - screened_rates(k_p_P31_to_O16_O16)*Y(H1)*state.rho - screened_rates(k_p_P31_to_S32)*Y(H1)*state.rho;
    jac.set(H1, P31, scratch);

    scratch = screened_rates(k_He4_S32_to_p_Cl35)*Y(He4)*state.rho + screened_rates(k_S32_to_p_P31);
    jac.set(H1, S32, scratch);

    scratch = -screened_rates(k_p_Cl35_to_Ar36)*Y(H1)*state.rho - screened_rates(k_p_Cl35_to_He4_S32)*Y(H1)*state.rho;
    jac.set(H1, Cl35, scratch);

    scratch = screened_rates(k_Ar36_to_p_Cl35) + screened_rates(k_He4_Ar36_to_p_K39)*Y(He4)*state.rho;
    jac.set(H1, Ar36, scratch);

    scratch = -screened_rates(k_p_K39_to_Ca40)*Y(H1)*state.rho - screened_rates(k_p_K39_to_He4_Ar36)*Y(H1)*state.rho;
    jac.set(H1, K39, scratch);

    scratch = screened_rates(k_Ca40_to_p_K39) + screened_rates(k_He4_Ca40_to_p_Sc43)*Y(He4)*state.rho;
    jac.set(H1, Ca40, scratch);

    scratch = -screened_rates(k_p_Sc43_to_He4_Ca40)*Y(H1)*state.rho - screened_rates(k_p_Sc43_to_Ti44)*Y(H1)*state.rho;
    jac.set(H1, Sc43, scratch);

    scratch = screened_rates(k_He4_Ti44_to_p_V47)*Y(He4)*state.rho + screened_rates(k_Ti44_to_p_Sc43);
    jac.set(H1, Ti44, scratch);

    scratch = -screened_rates(k_p_V47_to_Cr48)*Y(H1)*state.rho - screened_rates(k_p_V47_to_He4_Ti44)*Y(H1)*state.rho;
    jac.set(H1, V47, scratch);

    scratch = screened_rates(k_Cr48_to_p_V47) + screened_rates(k_He4_Cr48_to_p_Mn51)*Y(He4)*state.rho;
    jac.set(H1, Cr48, scratch);

    scratch = -screened_rates(k_p_Mn51_to_Fe52)*Y(H1)*state.rho - screened_rates(k_p_Mn51_to_He4_Cr48)*Y(H1)*state.rho;
    jac.set(H1, Mn51, scratch);

    scratch = screened_rates(k_Fe52_to_p_Mn51) + screened_rates(k_He4_Fe52_to_p_Co55)*Y(He4)*state.rho;
    jac.set(H1, Fe52, scratch);

    scratch = -screened_rates(k_p_Co55_to_He4_Fe52)*Y(H1)*state.rho - screened_rates(k_p_Co55_to_Ni56)*Y(H1)*state.rho;
    jac.set(H1, Co55, scratch);

    scratch = screened_rates(k_Ni56_to_p_Co55);
    jac.set(H1, Ni56, scratch);

    scratch = screened_rates(k_p_Al27_to_He4_Mg24)*Y(Al27)*state.rho + screened_rates(k_p_Cl35_to_He4_S32)*Y(Cl35)*state.rho + screened_rates(k_p_Co55_to_He4_Fe52)*Y(Co55)*state.rho + screened_rates(k_p_K39_to_He4_Ar36)*Y(K39)*state.rho + screened_rates(k_p_Mn51_to_He4_Cr48)*Y(Mn51)*state.rho + screened_rates(k_p_Na23_to_He4_Ne20)*Y(Na23)*state.rho + screened_rates(k_p_Ne21_to_He4_F18)*Y(Ne21)*state.rho + screened_rates(k_p_O16_to_He4_N13)*Y(O16)*state.rho + screened_rates(k_p_P31_to_He4_Si28)*Y(P31)*state.rho + screened_rates(k_p_Sc43_to_He4_Ca40)*Y(Sc43)*state.rho + screened_rates(k_p_V47_to_He4_Ti44)*Y(V47)*state.rho;
    jac.set(He4, H1, scratch);

    scratch = -screened_rates(k_He4_Al27_to_P31)*Y(Al27)*state.rho - screened_rates(k_He4_Ar36_to_Ca40)*Y(Ar36)*state.rho - screened_rates(k_He4_Ar36_to_p_K39)*Y(Ar36)*state.rho - screened_rates(k_He4_C12_to_O16)*Y(C12)*state.rho - screened_rates(k_He4_Ca40_to_Ti44)*Y(Ca40)*state.rho - screened_rates(k_He4_Ca40_to_p_Sc43)*Y(Ca40)*state.rho - screened_rates(k_He4_Cl35_to_K39)*Y(Cl35)*state.rho - screened_rates(k_He4_Cr48_to_Fe52)*Y(Cr48)*state.rho - screened_rates(k_He4_Cr48_to_p_Mn51)*Y(Cr48)*state.rho - screened_rates(k_He4_F18_to_Na22)*Y(F18)*state.rho - screened_rates(k_He4_F18_to_p_Ne21)*Y(F18)*state.rho - screened_rates(k_He4_Fe52_to_Ni56)*Y(Fe52)*state.rho - screened_rates(k_He4_Fe52_to_p_Co55)*Y(Fe52)*state.rho - 1.5*screened_rates(k_He4_He4_He4_to_C12)*std::pow(Y(He4), 2)*std::pow(state.rho, 2) - screened_rates(k_He4_K39_to_Sc43)*Y(K39)*state.rho - screened_rates(k_He4_Mg24_to_C12_O16)*Y(Mg24)*state.rho - screened_rates(k_He4_Mg24_to_Si28)*Y(Mg24)*state.rho - screened_rates(k_He4_Mg24_to_p_Al27)*Y(Mg24)*state.rho - screened_rates(k_He4_Mn51_to_Co55)*Y(Mn51)*state.rho - screened_rates(k_He4_N13_to_p_O16)*Y(N13)*state.rho - screened_rates(k_He4_N14_to_F18)*Y(N14)*state.rho - screened_rates(k_He4_Na23_to_Al27)*Y(Na23)*state.rho - screened_rates(k_He4_Ne20_to_C12_C12)*Y(Ne20)*state.rho - screened_rates(k_He4_Ne20_to_Mg24)*Y(Ne20)*state.rho - screened_rates(k_He4_Ne20_to_p_Na23)*Y(Ne20)*state.rho - screened_rates(k_He4_O16_to_Ne20)*Y(O16)*state.rho - screened_rates(k_He4_P31_to_Cl35)*Y(P31)*state.rho - screened_rates(k_He4_S32_to_Ar36)*Y(S32)*state.rho - screened_rates(k_He4_S32_to_p_Cl35)*Y(S32)*state.rho - screened_rates(k_He4_Sc43_to_V47)*Y(Sc43)*state.rho - screened_rates(k_He4_Si28_to_C12_Ne20)*Y(Si28)*state.rho - screened_rates(k_He4_Si28_to_O16_O16)*Y(Si28)*state.rho - screened_rates(k_He4_Si28_to_S32)*Y(Si28)*state.rho - screened_rates(k_He4_Si28_to_p_P31)*Y(Si28)*state.rho - screened_rates(k_He4_Ti44_to_Cr48)*Y(Ti44)*state.rho - screened_rates(k_He4_Ti44_to_p_V47)*Y(Ti44)*state.rho - screened_rates(k_He4_V47_to_Mn51)*Y(V47)*state.rho;
    jac.set(He4, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_He4_Ne20)*Y(C12)*state.rho + screened_rates(k_C12_Ne20_to_He4_Si28)*Y(Ne20)*state.rho + screened_rates(k_C12_O16_to_He4_Mg24)*Y(O16)*state.rho + 3.0*screened_rates(k_C12_to_He4_He4_He4) - screened_rates(k_He4_C12_to_O16)*Y(He4)*state.rho;
    jac.set(He4, C12, scratch);

    scratch = -screened_rates(k_He4_N13_to_p_O16)*Y(He4)*state.rho;
    jac.set(He4, N13, scratch);

    scratch = -screened_rates(k_He4_N14_to_F18)*Y(He4)*state.rho;
    jac.set(He4, N14, scratch);

    scratch = screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*state.rho - screened_rates(k_He4_O16_to_Ne20)*Y(He4)*state.rho + 1.0*screened_rates(k_O16_O16_to_He4_Si28)*Y(O16)*state.rho + screened_rates(k_O16_to_He4_C12) + screened_rates(k_p_O16_to_He4_N13)*Y(H1)*state.rho;
    jac.set(He4, O16, scratch);

    scratch = screened_rates(k_F18_to_He4_N14) - screened_rates(k_He4_F18_to_Na22)*Y(He4)*state.rho - screened_rates(k_He4_F18_to_p_Ne21)*Y(He4)*state.rho;
    jac.set(He4, F18, scratch);

    scratch = screened_rates(k_C12_Ne20_to_He4_Si28)*Y(C12)*state.rho - screened_rates(k_He4_Ne20_to_C12_C12)*Y(He4)*state.rho - screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*state.rho - screened_rates(k_He4_Ne20_to_p_Na23)*Y(He4)*state.rho + screened_rates(k_Ne20_to_He4_O16);
    jac.set(He4, Ne20, scratch);

    scratch = screened_rates(k_p_Ne21_to_He4_F18)*Y(H1)*state.rho;
    jac.set(He4, Ne21, scratch);

    scratch = screened_rates(k_Na22_to_He4_F18);
    jac.set(He4, Na22, scratch);

    scratch = -screened_rates(k_He4_Na23_to_Al27)*Y(He4)*state.rho + screened_rates(k_p_Na23_to_He4_Ne20)*Y(H1)*state.rho;
    jac.set(He4, Na23, scratch);

    scratch = -screened_rates(k_He4_Mg24_to_C12_O16)*Y(He4)*state.rho - screened_rates(k_He4_Mg24_to_Si28)*Y(He4)*state.rho - screened_rates(k_He4_Mg24_to_p_Al27)*Y(He4)*state.rho + screened_rates(k_Mg24_to_He4_Ne20);
    jac.set(He4, Mg24, scratch);

    scratch = screened_rates(k_Al27_to_He4_Na23) - screened_rates(k_He4_Al27_to_P31)*Y(He4)*state.rho + screened_rates(k_p_Al27_to_He4_Mg24)*Y(H1)*state.rho;
    jac.set(He4, Al27, scratch);

    scratch = -screened_rates(k_He4_Si28_to_C12_Ne20)*Y(He4)*state.rho - screened_rates(k_He4_Si28_to_O16_O16)*Y(He4)*state.rho - screened_rates(k_He4_Si28_to_S32)*Y(He4)*state.rho - screened_rates(k_He4_Si28_to_p_P31)*Y(He4)*state.rho + screened_rates(k_Si28_to_He4_Mg24);
    jac.set(He4, Si28, scratch);

    scratch = -screened_rates(k_He4_P31_to_Cl35)*Y(He4)*state.rho + screened_rates(k_P31_to_He4_Al27) + screened_rates(k_p_P31_to_He4_Si28)*Y(H1)*state.rho;
    jac.set(He4, P31, scratch);

    scratch = -screened_rates(k_He4_S32_to_Ar36)*Y(He4)*state.rho - screened_rates(k_He4_S32_to_p_Cl35)*Y(He4)*state.rho + screened_rates(k_S32_to_He4_Si28);
    jac.set(He4, S32, scratch);

    scratch = screened_rates(k_Cl35_to_He4_P31) - screened_rates(k_He4_Cl35_to_K39)*Y(He4)*state.rho + screened_rates(k_p_Cl35_to_He4_S32)*Y(H1)*state.rho;
    jac.set(He4, Cl35, scratch);

    scratch = screened_rates(k_Ar36_to_He4_S32) - screened_rates(k_He4_Ar36_to_Ca40)*Y(He4)*state.rho - screened_rates(k_He4_Ar36_to_p_K39)*Y(He4)*state.rho;
    jac.set(He4, Ar36, scratch);

    scratch = -screened_rates(k_He4_K39_to_Sc43)*Y(He4)*state.rho + screened_rates(k_K39_to_He4_Cl35) + screened_rates(k_p_K39_to_He4_Ar36)*Y(H1)*state.rho;
    jac.set(He4, K39, scratch);

    scratch = screened_rates(k_Ca40_to_He4_Ar36) - screened_rates(k_He4_Ca40_to_Ti44)*Y(He4)*state.rho - screened_rates(k_He4_Ca40_to_p_Sc43)*Y(He4)*state.rho;
    jac.set(He4, Ca40, scratch);

    scratch = -screened_rates(k_He4_Sc43_to_V47)*Y(He4)*state.rho + screened_rates(k_Sc43_to_He4_K39) + screened_rates(k_p_Sc43_to_He4_Ca40)*Y(H1)*state.rho;
    jac.set(He4, Sc43, scratch);

    scratch = -screened_rates(k_He4_Ti44_to_Cr48)*Y(He4)*state.rho - screened_rates(k_He4_Ti44_to_p_V47)*Y(He4)*state.rho + screened_rates(k_Ti44_to_He4_Ca40);
    jac.set(He4, Ti44, scratch);

    scratch = -screened_rates(k_He4_V47_to_Mn51)*Y(He4)*state.rho + screened_rates(k_V47_to_He4_Sc43) + screened_rates(k_p_V47_to_He4_Ti44)*Y(H1)*state.rho;
    jac.set(He4, V47, scratch);

    scratch = screened_rates(k_Cr48_to_He4_Ti44) - screened_rates(k_He4_Cr48_to_Fe52)*Y(He4)*state.rho - screened_rates(k_He4_Cr48_to_p_Mn51)*Y(He4)*state.rho;
    jac.set(He4, Cr48, scratch);

    scratch = -screened_rates(k_He4_Mn51_to_Co55)*Y(He4)*state.rho + screened_rates(k_Mn51_to_He4_V47) + screened_rates(k_p_Mn51_to_He4_Cr48)*Y(H1)*state.rho;
    jac.set(He4, Mn51, scratch);

    scratch = screened_rates(k_Fe52_to_He4_Cr48) - screened_rates(k_He4_Fe52_to_Ni56)*Y(He4)*state.rho - screened_rates(k_He4_Fe52_to_p_Co55)*Y(He4)*state.rho;
    jac.set(He4, Fe52, scratch);

    scratch = screened_rates(k_Co55_to_He4_Mn51) + screened_rates(k_p_Co55_to_He4_Fe52)*Y(H1)*state.rho;
    jac.set(He4, Co55, scratch);

    scratch = screened_rates(k_Ni56_to_He4_Fe52);
    jac.set(He4, Ni56, scratch);

    scratch = screened_rates(k_p_Al27_to_C12_O16)*Y(Al27)*state.rho - screened_rates(k_p_C12_to_N13)*Y(C12)*state.rho + 2.0*screened_rates(k_p_Na23_to_C12_C12)*Y(Na23)*state.rho + screened_rates(k_p_P31_to_C12_Ne20)*Y(P31)*state.rho;
    jac.set(C12, H1, scratch);

    scratch = -screened_rates(k_He4_C12_to_O16)*Y(C12)*state.rho + 0.5*screened_rates(k_He4_He4_He4_to_C12)*std::pow(Y(He4), 2)*std::pow(state.rho, 2) + screened_rates(k_He4_Mg24_to_C12_O16)*Y(Mg24)*state.rho + 2.0*screened_rates(k_He4_Ne20_to_C12_C12)*Y(Ne20)*state.rho + screened_rates(k_He4_Si28_to_C12_Ne20)*Y(Si28)*state.rho;
    jac.set(C12, He4, scratch);

    scratch = -2.0*screened_rates(k_C12_C12_to_He4_Ne20)*Y(C12)*state.rho - 2.0*screened_rates(k_C12_C12_to_Mg24_modified)*Y(C12)*state.rho - 2.0*screened_rates(k_C12_C12_to_p_Na23)*Y(C12)*state.rho - screened_rates(k_C12_Ne20_to_He4_Si28)*Y(Ne20)*state.rho - screened_rates(k_C12_Ne20_to_p_P31)*Y(Ne20)*state.rho - screened_rates(k_C12_O16_to_He4_Mg24)*Y(O16)*state.rho - screened_rates(k_C12_O16_to_Si28_modified)*Y(O16)*state.rho - screened_rates(k_C12_O16_to_p_Al27)*Y(O16)*state.rho - screened_rates(k_C12_to_He4_He4_He4) - screened_rates(k_He4_C12_to_O16)*Y(He4)*state.rho - screened_rates(k_p_C12_to_N13)*Y(H1)*state.rho;
    jac.set(C12, C12, scratch);

    scratch = screened_rates(k_N13_to_p_C12);
    jac.set(C12, N13, scratch);

    scratch = -screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*state.rho - screened_rates(k_C12_O16_to_Si28_modified)*Y(C12)*state.rho - screened_rates(k_C12_O16_to_p_Al27)*Y(C12)*state.rho + screened_rates(k_O16_to_He4_C12);
    jac.set(C12, O16, scratch);

    scratch = -screened_rates(k_C12_Ne20_to_He4_Si28)*Y(C12)*state.rho - screened_rates(k_C12_Ne20_to_p_P31)*Y(C12)*state.rho + 2.0*screened_rates(k_He4_Ne20_to_C12_C12)*Y(He4)*state.rho;
    jac.set(C12, Ne20, scratch);

    scratch = 2.0*screened_rates(k_p_Na23_to_C12_C12)*Y(H1)*state.rho;
    jac.set(C12, Na23, scratch);

    scratch = screened_rates(k_He4_Mg24_to_C12_O16)*Y(He4)*state.rho;
    jac.set(C12, Mg24, scratch);

    scratch = screened_rates(k_p_Al27_to_C12_O16)*Y(H1)*state.rho;
    jac.set(C12, Al27, scratch);

    scratch = screened_rates(k_He4_Si28_to_C12_Ne20)*Y(He4)*state.rho;
    jac.set(C12, Si28, scratch);

    scratch = screened_rates(k_p_P31_to_C12_Ne20)*Y(H1)*state.rho;
    jac.set(C12, P31, scratch);

    scratch = screened_rates(k_p_C12_to_N13)*Y(C12)*state.rho + screened_rates(k_p_O16_to_He4_N13)*Y(O16)*state.rho;
    jac.set(N13, H1, scratch);

    scratch = -screened_rates(k_He4_N13_to_p_O16)*Y(N13)*state.rho;
    jac.set(N13, He4, scratch);

    scratch = screened_rates(k_p_C12_to_N13)*Y(H1)*state.rho;
    jac.set(N13, C12, scratch);

    scratch = -screened_rates(k_He4_N13_to_p_O16)*Y(He4)*state.rho - screened_rates(k_N13_to_p_C12);
    jac.set(N13, N13, scratch);

    scratch = screened_rates(k_p_O16_to_He4_N13)*Y(H1)*state.rho;
    jac.set(N13, O16, scratch);

    scratch = -screened_rates(k_He4_N14_to_F18)*Y(N14)*state.rho;
    jac.set(N14, He4, scratch);

    scratch = -screened_rates(k_He4_N14_to_F18)*Y(He4)*state.rho;
    jac.set(N14, N14, scratch);

    scratch = screened_rates(k_F18_to_He4_N14);
    jac.set(N14, F18, scratch);

    scratch = screened_rates(k_p_Al27_to_C12_O16)*Y(Al27)*state.rho - screened_rates(k_p_O16_to_He4_N13)*Y(O16)*state.rho + 2.0*screened_rates(k_p_P31_to_O16_O16)*Y(P31)*state.rho;
    jac.set(O16, H1, scratch);

    scratch = screened_rates(k_He4_C12_to_O16)*Y(C12)*state.rho + screened_rates(k_He4_Mg24_to_C12_O16)*Y(Mg24)*state.rho + screened_rates(k_He4_N13_to_p_O16)*Y(N13)*state.rho - screened_rates(k_He4_O16_to_Ne20)*Y(O16)*state.rho + 2.0*screened_rates(k_He4_Si28_to_O16_O16)*Y(Si28)*state.rho;
    jac.set(O16, He4, scratch);

    scratch = -screened_rates(k_C12_O16_to_He4_Mg24)*Y(O16)*state.rho - screened_rates(k_C12_O16_to_Si28_modified)*Y(O16)*state.rho - screened_rates(k_C12_O16_to_p_Al27)*Y(O16)*state.rho + screened_rates(k_He4_C12_to_O16)*Y(He4)*state.rho;
    jac.set(O16, C12, scratch);

    scratch = screened_rates(k_He4_N13_to_p_O16)*Y(He4)*state.rho;
    jac.set(O16, N13, scratch);

    scratch = -screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*state.rho - screened_rates(k_C12_O16_to_Si28_modified)*Y(C12)*state.rho - screened_rates(k_C12_O16_to_p_Al27)*Y(C12)*state.rho - screened_rates(k_He4_O16_to_Ne20)*Y(He4)*state.rho - 2.0*screened_rates(k_O16_O16_to_He4_Si28)*Y(O16)*state.rho - 2.0*screened_rates(k_O16_O16_to_S32_modified)*Y(O16)*state.rho - 2.0*screened_rates(k_O16_O16_to_p_P31)*Y(O16)*state.rho - screened_rates(k_O16_to_He4_C12) - screened_rates(k_p_O16_to_He4_N13)*Y(H1)*state.rho;
    jac.set(O16, O16, scratch);

    scratch = screened_rates(k_Ne20_to_He4_O16);
    jac.set(O16, Ne20, scratch);

    scratch = screened_rates(k_He4_Mg24_to_C12_O16)*Y(He4)*state.rho;
    jac.set(O16, Mg24, scratch);

    scratch = screened_rates(k_p_Al27_to_C12_O16)*Y(H1)*state.rho;
    jac.set(O16, Al27, scratch);

    scratch = 2.0*screened_rates(k_He4_Si28_to_O16_O16)*Y(He4)*state.rho;
    jac.set(O16, Si28, scratch);

    scratch = 2.0*screened_rates(k_p_P31_to_O16_O16)*Y(H1)*state.rho;
    jac.set(O16, P31, scratch);

    scratch = screened_rates(k_p_Ne21_to_He4_F18)*Y(Ne21)*state.rho;
    jac.set(F18, H1, scratch);

    scratch = -screened_rates(k_He4_F18_to_Na22)*Y(F18)*state.rho - screened_rates(k_He4_F18_to_p_Ne21)*Y(F18)*state.rho + screened_rates(k_He4_N14_to_F18)*Y(N14)*state.rho;
    jac.set(F18, He4, scratch);

    scratch = screened_rates(k_He4_N14_to_F18)*Y(He4)*state.rho;
    jac.set(F18, N14, scratch);

    scratch = -screened_rates(k_F18_to_He4_N14) - screened_rates(k_He4_F18_to_Na22)*Y(He4)*state.rho - screened_rates(k_He4_F18_to_p_Ne21)*Y(He4)*state.rho;
    jac.set(F18, F18, scratch);

    scratch = screened_rates(k_p_Ne21_to_He4_F18)*Y(H1)*state.rho;
    jac.set(F18, Ne21, scratch);

    scratch = screened_rates(k_Na22_to_He4_F18);
    jac.set(F18, Na22, scratch);

    scratch = screened_rates(k_p_Na23_to_He4_Ne20)*Y(Na23)*state.rho + screened_rates(k_p_P31_to_C12_Ne20)*Y(P31)*state.rho;
    jac.set(Ne20, H1, scratch);

    scratch = -screened_rates(k_He4_Ne20_to_C12_C12)*Y(Ne20)*state.rho - screened_rates(k_He4_Ne20_to_Mg24)*Y(Ne20)*state.rho - screened_rates(k_He4_Ne20_to_p_Na23)*Y(Ne20)*state.rho + screened_rates(k_He4_O16_to_Ne20)*Y(O16)*state.rho + screened_rates(k_He4_Si28_to_C12_Ne20)*Y(Si28)*state.rho;
    jac.set(Ne20, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_He4_Ne20)*Y(C12)*state.rho - screened_rates(k_C12_Ne20_to_He4_Si28)*Y(Ne20)*state.rho - screened_rates(k_C12_Ne20_to_p_P31)*Y(Ne20)*state.rho;
    jac.set(Ne20, C12, scratch);

    scratch = screened_rates(k_He4_O16_to_Ne20)*Y(He4)*state.rho;
    jac.set(Ne20, O16, scratch);

    scratch = -screened_rates(k_C12_Ne20_to_He4_Si28)*Y(C12)*state.rho - screened_rates(k_C12_Ne20_to_p_P31)*Y(C12)*state.rho - screened_rates(k_He4_Ne20_to_C12_C12)*Y(He4)*state.rho - screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*state.rho - screened_rates(k_He4_Ne20_to_p_Na23)*Y(He4)*state.rho - screened_rates(k_Ne20_to_He4_O16);
    jac.set(Ne20, Ne20, scratch);

    scratch = screened_rates(k_p_Na23_to_He4_Ne20)*Y(H1)*state.rho;
    jac.set(Ne20, Na23, scratch);

    scratch = screened_rates(k_Mg24_to_He4_Ne20);
    jac.set(Ne20, Mg24, scratch);

    scratch = screened_rates(k_He4_Si28_to_C12_Ne20)*Y(He4)*state.rho;
    jac.set(Ne20, Si28, scratch);

    scratch = screened_rates(k_p_P31_to_C12_Ne20)*Y(H1)*state.rho;
    jac.set(Ne20, P31, scratch);

    scratch = -screened_rates(k_p_Ne21_to_He4_F18)*Y(Ne21)*state.rho - screened_rates(k_p_Ne21_to_Na22)*Y(Ne21)*state.rho;
    jac.set(Ne21, H1, scratch);

    scratch = screened_rates(k_He4_F18_to_p_Ne21)*Y(F18)*state.rho;
    jac.set(Ne21, He4, scratch);

    scratch = screened_rates(k_He4_F18_to_p_Ne21)*Y(He4)*state.rho;
    jac.set(Ne21, F18, scratch);

    scratch = -screened_rates(k_p_Ne21_to_He4_F18)*Y(H1)*state.rho - screened_rates(k_p_Ne21_to_Na22)*Y(H1)*state.rho;
    jac.set(Ne21, Ne21, scratch);

    scratch = screened_rates(k_Na22_to_p_Ne21);
    jac.set(Ne21, Na22, scratch);

    scratch = screened_rates(k_p_Ne21_to_Na22)*Y(Ne21)*state.rho;
    jac.set(Na22, H1, scratch);

    scratch = screened_rates(k_He4_F18_to_Na22)*Y(F18)*state.rho;
    jac.set(Na22, He4, scratch);

    scratch = screened_rates(k_He4_F18_to_Na22)*Y(He4)*state.rho;
    jac.set(Na22, F18, scratch);

    scratch = screened_rates(k_p_Ne21_to_Na22)*Y(H1)*state.rho;
    jac.set(Na22, Ne21, scratch);

    scratch = -screened_rates(k_Na22_to_He4_F18) - screened_rates(k_Na22_to_p_Ne21);
    jac.set(Na22, Na22, scratch);

    scratch = -screened_rates(k_p_Na23_to_C12_C12)*Y(Na23)*state.rho - screened_rates(k_p_Na23_to_He4_Ne20)*Y(Na23)*state.rho - screened_rates(k_p_Na23_to_Mg24)*Y(Na23)*state.rho;
    jac.set(Na23, H1, scratch);

    scratch = -screened_rates(k_He4_Na23_to_Al27)*Y(Na23)*state.rho + screened_rates(k_He4_Ne20_to_p_Na23)*Y(Ne20)*state.rho;
    jac.set(Na23, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_p_Na23)*Y(C12)*state.rho;
    jac.set(Na23, C12, scratch);

    scratch = screened_rates(k_He4_Ne20_to_p_Na23)*Y(He4)*state.rho;
    jac.set(Na23, Ne20, scratch);

    scratch = -screened_rates(k_He4_Na23_to_Al27)*Y(He4)*state.rho - screened_rates(k_p_Na23_to_C12_C12)*Y(H1)*state.rho - screened_rates(k_p_Na23_to_He4_Ne20)*Y(H1)*state.rho - screened_rates(k_p_Na23_to_Mg24)*Y(H1)*state.rho;
    jac.set(Na23, Na23, scratch);

    scratch = screened_rates(k_Mg24_to_p_Na23);
    jac.set(Na23, Mg24, scratch);

    scratch = screened_rates(k_Al27_to_He4_Na23);
    jac.set(Na23, Al27, scratch);

    scratch = screened_rates(k_p_Al27_to_He4_Mg24)*Y(Al27)*state.rho + screened_rates(k_p_Na23_to_Mg24)*Y(Na23)*state.rho;
    jac.set(Mg24, H1, scratch);

    scratch = -screened_rates(k_He4_Mg24_to_C12_O16)*Y(Mg24)*state.rho - screened_rates(k_He4_Mg24_to_Si28)*Y(Mg24)*state.rho - screened_rates(k_He4_Mg24_to_p_Al27)*Y(Mg24)*state.rho + screened_rates(k_He4_Ne20_to_Mg24)*Y(Ne20)*state.rho;
    jac.set(Mg24, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_Mg24_modified)*Y(C12)*state.rho + screened_rates(k_C12_O16_to_He4_Mg24)*Y(O16)*state.rho;
    jac.set(Mg24, C12, scratch);

    scratch = screened_rates(k_C12_O16_to_He4_Mg24)*Y(C12)*state.rho;
    jac.set(Mg24, O16, scratch);

    scratch = screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*state.rho;
    jac.set(Mg24, Ne20, scratch);

    scratch = screened_rates(k_p_Na23_to_Mg24)*Y(H1)*state.rho;
    jac.set(Mg24, Na23, scratch);

    scratch = -screened_rates(k_He4_Mg24_to_C12_O16)*Y(He4)*state.rho - screened_rates(k_He4_Mg24_to_Si28)*Y(He4)*state.rho - screened_rates(k_He4_Mg24_to_p_Al27)*Y(He4)*state.rho - screened_rates(k_Mg24_to_He4_Ne20) - screened_rates(k_Mg24_to_p_Na23);
    jac.set(Mg24, Mg24, scratch);

    scratch = screened_rates(k_p_Al27_to_He4_Mg24)*Y(H1)*state.rho;
    jac.set(Mg24, Al27, scratch);

    scratch = screened_rates(k_Si28_to_He4_Mg24);
    jac.set(Mg24, Si28, scratch);

    scratch = -screened_rates(k_p_Al27_to_C12_O16)*Y(Al27)*state.rho - screened_rates(k_p_Al27_to_He4_Mg24)*Y(Al27)*state.rho - screened_rates(k_p_Al27_to_Si28)*Y(Al27)*state.rho;
    jac.set(Al27, H1, scratch);

    scratch = -screened_rates(k_He4_Al27_to_P31)*Y(Al27)*state.rho + screened_rates(k_He4_Mg24_to_p_Al27)*Y(Mg24)*state.rho + screened_rates(k_He4_Na23_to_Al27)*Y(Na23)*state.rho;
    jac.set(Al27, He4, scratch);

    scratch = screened_rates(k_C12_O16_to_p_Al27)*Y(O16)*state.rho;
    jac.set(Al27, C12, scratch);

    scratch = screened_rates(k_C12_O16_to_p_Al27)*Y(C12)*state.rho;
    jac.set(Al27, O16, scratch);

    scratch = screened_rates(k_He4_Na23_to_Al27)*Y(He4)*state.rho;
    jac.set(Al27, Na23, scratch);

    scratch = screened_rates(k_He4_Mg24_to_p_Al27)*Y(He4)*state.rho;
    jac.set(Al27, Mg24, scratch);

    scratch = -screened_rates(k_Al27_to_He4_Na23) - screened_rates(k_He4_Al27_to_P31)*Y(He4)*state.rho - screened_rates(k_p_Al27_to_C12_O16)*Y(H1)*state.rho - screened_rates(k_p_Al27_to_He4_Mg24)*Y(H1)*state.rho - screened_rates(k_p_Al27_to_Si28)*Y(H1)*state.rho;
    jac.set(Al27, Al27, scratch);

    scratch = screened_rates(k_Si28_to_p_Al27);
    jac.set(Al27, Si28, scratch);

    scratch = screened_rates(k_P31_to_He4_Al27);
    jac.set(Al27, P31, scratch);

    scratch = screened_rates(k_p_Al27_to_Si28)*Y(Al27)*state.rho + screened_rates(k_p_P31_to_He4_Si28)*Y(P31)*state.rho;
    jac.set(Si28, H1, scratch);

    scratch = screened_rates(k_He4_Mg24_to_Si28)*Y(Mg24)*state.rho - screened_rates(k_He4_Si28_to_C12_Ne20)*Y(Si28)*state.rho - screened_rates(k_He4_Si28_to_O16_O16)*Y(Si28)*state.rho - screened_rates(k_He4_Si28_to_S32)*Y(Si28)*state.rho - screened_rates(k_He4_Si28_to_p_P31)*Y(Si28)*state.rho;
    jac.set(Si28, He4, scratch);

    scratch = screened_rates(k_C12_Ne20_to_He4_Si28)*Y(Ne20)*state.rho + screened_rates(k_C12_O16_to_Si28_modified)*Y(O16)*state.rho;
    jac.set(Si28, C12, scratch);

    scratch = screened_rates(k_C12_O16_to_Si28_modified)*Y(C12)*state.rho + 1.0*screened_rates(k_O16_O16_to_He4_Si28)*Y(O16)*state.rho;
    jac.set(Si28, O16, scratch);

    scratch = screened_rates(k_C12_Ne20_to_He4_Si28)*Y(C12)*state.rho;
    jac.set(Si28, Ne20, scratch);

    scratch = screened_rates(k_He4_Mg24_to_Si28)*Y(He4)*state.rho;
    jac.set(Si28, Mg24, scratch);

    scratch = screened_rates(k_p_Al27_to_Si28)*Y(H1)*state.rho;
    jac.set(Si28, Al27, scratch);

    scratch = -screened_rates(k_He4_Si28_to_C12_Ne20)*Y(He4)*state.rho - screened_rates(k_He4_Si28_to_O16_O16)*Y(He4)*state.rho - screened_rates(k_He4_Si28_to_S32)*Y(He4)*state.rho - screened_rates(k_He4_Si28_to_p_P31)*Y(He4)*state.rho - screened_rates(k_Si28_to_He4_Mg24) - screened_rates(k_Si28_to_p_Al27);
    jac.set(Si28, Si28, scratch);

    scratch = screened_rates(k_p_P31_to_He4_Si28)*Y(H1)*state.rho;
    jac.set(Si28, P31, scratch);

    scratch = screened_rates(k_S32_to_He4_Si28);
    jac.set(Si28, S32, scratch);

    scratch = -screened_rates(k_p_P31_to_C12_Ne20)*Y(P31)*state.rho - screened_rates(k_p_P31_to_He4_Si28)*Y(P31)*state.rho - screened_rates(k_p_P31_to_O16_O16)*Y(P31)*state.rho - screened_rates(k_p_P31_to_S32)*Y(P31)*state.rho;
    jac.set(P31, H1, scratch);

    scratch = screened_rates(k_He4_Al27_to_P31)*Y(Al27)*state.rho - screened_rates(k_He4_P31_to_Cl35)*Y(P31)*state.rho + screened_rates(k_He4_Si28_to_p_P31)*Y(Si28)*state.rho;
    jac.set(P31, He4, scratch);

    scratch = screened_rates(k_C12_Ne20_to_p_P31)*Y(Ne20)*state.rho;
    jac.set(P31, C12, scratch);

    scratch = 1.0*screened_rates(k_O16_O16_to_p_P31)*Y(O16)*state.rho;
    jac.set(P31, O16, scratch);

    scratch = screened_rates(k_C12_Ne20_to_p_P31)*Y(C12)*state.rho;
    jac.set(P31, Ne20, scratch);

    scratch = screened_rates(k_He4_Al27_to_P31)*Y(He4)*state.rho;
    jac.set(P31, Al27, scratch);

    scratch = screened_rates(k_He4_Si28_to_p_P31)*Y(He4)*state.rho;
    jac.set(P31, Si28, scratch);

    scratch = -screened_rates(k_He4_P31_to_Cl35)*Y(He4)*state.rho - screened_rates(k_P31_to_He4_Al27) - screened_rates(k_p_P31_to_C12_Ne20)*Y(H1)*state.rho - screened_rates(k_p_P31_to_He4_Si28)*Y(H1)*state.rho - screened_rates(k_p_P31_to_O16_O16)*Y(H1)*state.rho - screened_rates(k_p_P31_to_S32)*Y(H1)*state.rho;
    jac.set(P31, P31, scratch);

    scratch = screened_rates(k_S32_to_p_P31);
    jac.set(P31, S32, scratch);

    scratch = screened_rates(k_Cl35_to_He4_P31);
    jac.set(P31, Cl35, scratch);

    scratch = screened_rates(k_p_Cl35_to_He4_S32)*Y(Cl35)*state.rho + screened_rates(k_p_P31_to_S32)*Y(P31)*state.rho;
    jac.set(S32, H1, scratch);

    scratch = -screened_rates(k_He4_S32_to_Ar36)*Y(S32)*state.rho - screened_rates(k_He4_S32_to_p_Cl35)*Y(S32)*state.rho + screened_rates(k_He4_Si28_to_S32)*Y(Si28)*state.rho;
    jac.set(S32, He4, scratch);

    scratch = 1.0*screened_rates(k_O16_O16_to_S32_modified)*Y(O16)*state.rho;
    jac.set(S32, O16, scratch);

    scratch = screened_rates(k_He4_Si28_to_S32)*Y(He4)*state.rho;
    jac.set(S32, Si28, scratch);

    scratch = screened_rates(k_p_P31_to_S32)*Y(H1)*state.rho;
    jac.set(S32, P31, scratch);

    scratch = -screened_rates(k_He4_S32_to_Ar36)*Y(He4)*state.rho - screened_rates(k_He4_S32_to_p_Cl35)*Y(He4)*state.rho - screened_rates(k_S32_to_He4_Si28) - screened_rates(k_S32_to_p_P31);
    jac.set(S32, S32, scratch);

    scratch = screened_rates(k_p_Cl35_to_He4_S32)*Y(H1)*state.rho;
    jac.set(S32, Cl35, scratch);

    scratch = screened_rates(k_Ar36_to_He4_S32);
    jac.set(S32, Ar36, scratch);

    scratch = -screened_rates(k_p_Cl35_to_Ar36)*Y(Cl35)*state.rho - screened_rates(k_p_Cl35_to_He4_S32)*Y(Cl35)*state.rho;
    jac.set(Cl35, H1, scratch);

    scratch = -screened_rates(k_He4_Cl35_to_K39)*Y(Cl35)*state.rho + screened_rates(k_He4_P31_to_Cl35)*Y(P31)*state.rho + screened_rates(k_He4_S32_to_p_Cl35)*Y(S32)*state.rho;
    jac.set(Cl35, He4, scratch);

    scratch = screened_rates(k_He4_P31_to_Cl35)*Y(He4)*state.rho;
    jac.set(Cl35, P31, scratch);

    scratch = screened_rates(k_He4_S32_to_p_Cl35)*Y(He4)*state.rho;
    jac.set(Cl35, S32, scratch);

    scratch = -screened_rates(k_Cl35_to_He4_P31) - screened_rates(k_He4_Cl35_to_K39)*Y(He4)*state.rho - screened_rates(k_p_Cl35_to_Ar36)*Y(H1)*state.rho - screened_rates(k_p_Cl35_to_He4_S32)*Y(H1)*state.rho;
    jac.set(Cl35, Cl35, scratch);

    scratch = screened_rates(k_Ar36_to_p_Cl35);
    jac.set(Cl35, Ar36, scratch);

    scratch = screened_rates(k_K39_to_He4_Cl35);
    jac.set(Cl35, K39, scratch);

    scratch = screened_rates(k_p_Cl35_to_Ar36)*Y(Cl35)*state.rho + screened_rates(k_p_K39_to_He4_Ar36)*Y(K39)*state.rho;
    jac.set(Ar36, H1, scratch);

    scratch = -screened_rates(k_He4_Ar36_to_Ca40)*Y(Ar36)*state.rho - screened_rates(k_He4_Ar36_to_p_K39)*Y(Ar36)*state.rho + screened_rates(k_He4_S32_to_Ar36)*Y(S32)*state.rho;
    jac.set(Ar36, He4, scratch);

    scratch = screened_rates(k_He4_S32_to_Ar36)*Y(He4)*state.rho;
    jac.set(Ar36, S32, scratch);

    scratch = screened_rates(k_p_Cl35_to_Ar36)*Y(H1)*state.rho;
    jac.set(Ar36, Cl35, scratch);

    scratch = -screened_rates(k_Ar36_to_He4_S32) - screened_rates(k_Ar36_to_p_Cl35) - screened_rates(k_He4_Ar36_to_Ca40)*Y(He4)*state.rho - screened_rates(k_He4_Ar36_to_p_K39)*Y(He4)*state.rho;
    jac.set(Ar36, Ar36, scratch);

    scratch = screened_rates(k_p_K39_to_He4_Ar36)*Y(H1)*state.rho;
    jac.set(Ar36, K39, scratch);

    scratch = screened_rates(k_Ca40_to_He4_Ar36);
    jac.set(Ar36, Ca40, scratch);

    scratch = -screened_rates(k_p_K39_to_Ca40)*Y(K39)*state.rho - screened_rates(k_p_K39_to_He4_Ar36)*Y(K39)*state.rho;
    jac.set(K39, H1, scratch);

    scratch = screened_rates(k_He4_Ar36_to_p_K39)*Y(Ar36)*state.rho + screened_rates(k_He4_Cl35_to_K39)*Y(Cl35)*state.rho - screened_rates(k_He4_K39_to_Sc43)*Y(K39)*state.rho;
    jac.set(K39, He4, scratch);

    scratch = screened_rates(k_He4_Cl35_to_K39)*Y(He4)*state.rho;
    jac.set(K39, Cl35, scratch);

    scratch = screened_rates(k_He4_Ar36_to_p_K39)*Y(He4)*state.rho;
    jac.set(K39, Ar36, scratch);

    scratch = -screened_rates(k_He4_K39_to_Sc43)*Y(He4)*state.rho - screened_rates(k_K39_to_He4_Cl35) - screened_rates(k_p_K39_to_Ca40)*Y(H1)*state.rho - screened_rates(k_p_K39_to_He4_Ar36)*Y(H1)*state.rho;
    jac.set(K39, K39, scratch);

    scratch = screened_rates(k_Ca40_to_p_K39);
    jac.set(K39, Ca40, scratch);

    scratch = screened_rates(k_Sc43_to_He4_K39);
    jac.set(K39, Sc43, scratch);

    scratch = screened_rates(k_p_K39_to_Ca40)*Y(K39)*state.rho + screened_rates(k_p_Sc43_to_He4_Ca40)*Y(Sc43)*state.rho;
    jac.set(Ca40, H1, scratch);

    scratch = screened_rates(k_He4_Ar36_to_Ca40)*Y(Ar36)*state.rho - screened_rates(k_He4_Ca40_to_Ti44)*Y(Ca40)*state.rho - screened_rates(k_He4_Ca40_to_p_Sc43)*Y(Ca40)*state.rho;
    jac.set(Ca40, He4, scratch);

    scratch = screened_rates(k_He4_Ar36_to_Ca40)*Y(He4)*state.rho;
    jac.set(Ca40, Ar36, scratch);

    scratch = screened_rates(k_p_K39_to_Ca40)*Y(H1)*state.rho;
    jac.set(Ca40, K39, scratch);

    scratch = -screened_rates(k_Ca40_to_He4_Ar36) - screened_rates(k_Ca40_to_p_K39) - screened_rates(k_He4_Ca40_to_Ti44)*Y(He4)*state.rho - screened_rates(k_He4_Ca40_to_p_Sc43)*Y(He4)*state.rho;
    jac.set(Ca40, Ca40, scratch);

    scratch = screened_rates(k_p_Sc43_to_He4_Ca40)*Y(H1)*state.rho;
    jac.set(Ca40, Sc43, scratch);

    scratch = screened_rates(k_Ti44_to_He4_Ca40);
    jac.set(Ca40, Ti44, scratch);

    scratch = -screened_rates(k_p_Sc43_to_He4_Ca40)*Y(Sc43)*state.rho - screened_rates(k_p_Sc43_to_Ti44)*Y(Sc43)*state.rho;
    jac.set(Sc43, H1, scratch);

    scratch = screened_rates(k_He4_Ca40_to_p_Sc43)*Y(Ca40)*state.rho + screened_rates(k_He4_K39_to_Sc43)*Y(K39)*state.rho - screened_rates(k_He4_Sc43_to_V47)*Y(Sc43)*state.rho;
    jac.set(Sc43, He4, scratch);

    scratch = screened_rates(k_He4_K39_to_Sc43)*Y(He4)*state.rho;
    jac.set(Sc43, K39, scratch);

    scratch = screened_rates(k_He4_Ca40_to_p_Sc43)*Y(He4)*state.rho;
    jac.set(Sc43, Ca40, scratch);

    scratch = -screened_rates(k_He4_Sc43_to_V47)*Y(He4)*state.rho - screened_rates(k_Sc43_to_He4_K39) - screened_rates(k_p_Sc43_to_He4_Ca40)*Y(H1)*state.rho - screened_rates(k_p_Sc43_to_Ti44)*Y(H1)*state.rho;
    jac.set(Sc43, Sc43, scratch);

    scratch = screened_rates(k_Ti44_to_p_Sc43);
    jac.set(Sc43, Ti44, scratch);

    scratch = screened_rates(k_V47_to_He4_Sc43);
    jac.set(Sc43, V47, scratch);

    scratch = screened_rates(k_p_Sc43_to_Ti44)*Y(Sc43)*state.rho + screened_rates(k_p_V47_to_He4_Ti44)*Y(V47)*state.rho;
    jac.set(Ti44, H1, scratch);

    scratch = screened_rates(k_He4_Ca40_to_Ti44)*Y(Ca40)*state.rho - screened_rates(k_He4_Ti44_to_Cr48)*Y(Ti44)*state.rho - screened_rates(k_He4_Ti44_to_p_V47)*Y(Ti44)*state.rho;
    jac.set(Ti44, He4, scratch);

    scratch = screened_rates(k_He4_Ca40_to_Ti44)*Y(He4)*state.rho;
    jac.set(Ti44, Ca40, scratch);

    scratch = screened_rates(k_p_Sc43_to_Ti44)*Y(H1)*state.rho;
    jac.set(Ti44, Sc43, scratch);

    scratch = -screened_rates(k_He4_Ti44_to_Cr48)*Y(He4)*state.rho - screened_rates(k_He4_Ti44_to_p_V47)*Y(He4)*state.rho - screened_rates(k_Ti44_to_He4_Ca40) - screened_rates(k_Ti44_to_p_Sc43);
    jac.set(Ti44, Ti44, scratch);

    scratch = screened_rates(k_p_V47_to_He4_Ti44)*Y(H1)*state.rho;
    jac.set(Ti44, V47, scratch);

    scratch = screened_rates(k_Cr48_to_He4_Ti44);
    jac.set(Ti44, Cr48, scratch);

    scratch = -screened_rates(k_p_V47_to_Cr48)*Y(V47)*state.rho - screened_rates(k_p_V47_to_He4_Ti44)*Y(V47)*state.rho;
    jac.set(V47, H1, scratch);

    scratch = screened_rates(k_He4_Sc43_to_V47)*Y(Sc43)*state.rho + screened_rates(k_He4_Ti44_to_p_V47)*Y(Ti44)*state.rho - screened_rates(k_He4_V47_to_Mn51)*Y(V47)*state.rho;
    jac.set(V47, He4, scratch);

    scratch = screened_rates(k_He4_Sc43_to_V47)*Y(He4)*state.rho;
    jac.set(V47, Sc43, scratch);

    scratch = screened_rates(k_He4_Ti44_to_p_V47)*Y(He4)*state.rho;
    jac.set(V47, Ti44, scratch);

    scratch = -screened_rates(k_He4_V47_to_Mn51)*Y(He4)*state.rho - screened_rates(k_V47_to_He4_Sc43) - screened_rates(k_p_V47_to_Cr48)*Y(H1)*state.rho - screened_rates(k_p_V47_to_He4_Ti44)*Y(H1)*state.rho;
    jac.set(V47, V47, scratch);

    scratch = screened_rates(k_Cr48_to_p_V47);
    jac.set(V47, Cr48, scratch);

    scratch = screened_rates(k_Mn51_to_He4_V47);
    jac.set(V47, Mn51, scratch);

    scratch = screened_rates(k_p_Mn51_to_He4_Cr48)*Y(Mn51)*state.rho + screened_rates(k_p_V47_to_Cr48)*Y(V47)*state.rho;
    jac.set(Cr48, H1, scratch);

    scratch = -screened_rates(k_He4_Cr48_to_Fe52)*Y(Cr48)*state.rho - screened_rates(k_He4_Cr48_to_p_Mn51)*Y(Cr48)*state.rho + screened_rates(k_He4_Ti44_to_Cr48)*Y(Ti44)*state.rho;
    jac.set(Cr48, He4, scratch);

    scratch = screened_rates(k_He4_Ti44_to_Cr48)*Y(He4)*state.rho;
    jac.set(Cr48, Ti44, scratch);

    scratch = screened_rates(k_p_V47_to_Cr48)*Y(H1)*state.rho;
    jac.set(Cr48, V47, scratch);

    scratch = -screened_rates(k_Cr48_to_He4_Ti44) - screened_rates(k_Cr48_to_p_V47) - screened_rates(k_He4_Cr48_to_Fe52)*Y(He4)*state.rho - screened_rates(k_He4_Cr48_to_p_Mn51)*Y(He4)*state.rho;
    jac.set(Cr48, Cr48, scratch);

    scratch = screened_rates(k_p_Mn51_to_He4_Cr48)*Y(H1)*state.rho;
    jac.set(Cr48, Mn51, scratch);

    scratch = screened_rates(k_Fe52_to_He4_Cr48);
    jac.set(Cr48, Fe52, scratch);

    scratch = -screened_rates(k_p_Mn51_to_Fe52)*Y(Mn51)*state.rho - screened_rates(k_p_Mn51_to_He4_Cr48)*Y(Mn51)*state.rho;
    jac.set(Mn51, H1, scratch);

    scratch = screened_rates(k_He4_Cr48_to_p_Mn51)*Y(Cr48)*state.rho - screened_rates(k_He4_Mn51_to_Co55)*Y(Mn51)*state.rho + screened_rates(k_He4_V47_to_Mn51)*Y(V47)*state.rho;
    jac.set(Mn51, He4, scratch);

    scratch = screened_rates(k_He4_V47_to_Mn51)*Y(He4)*state.rho;
    jac.set(Mn51, V47, scratch);

    scratch = screened_rates(k_He4_Cr48_to_p_Mn51)*Y(He4)*state.rho;
    jac.set(Mn51, Cr48, scratch);

    scratch = -screened_rates(k_He4_Mn51_to_Co55)*Y(He4)*state.rho - screened_rates(k_Mn51_to_He4_V47) - screened_rates(k_p_Mn51_to_Fe52)*Y(H1)*state.rho - screened_rates(k_p_Mn51_to_He4_Cr48)*Y(H1)*state.rho;
    jac.set(Mn51, Mn51, scratch);

    scratch = screened_rates(k_Fe52_to_p_Mn51);
    jac.set(Mn51, Fe52, scratch);

    scratch = screened_rates(k_Co55_to_He4_Mn51);
    jac.set(Mn51, Co55, scratch);

    scratch = screened_rates(k_p_Co55_to_He4_Fe52)*Y(Co55)*state.rho + screened_rates(k_p_Mn51_to_Fe52)*Y(Mn51)*state.rho;
    jac.set(Fe52, H1, scratch);

    scratch = screened_rates(k_He4_Cr48_to_Fe52)*Y(Cr48)*state.rho - screened_rates(k_He4_Fe52_to_Ni56)*Y(Fe52)*state.rho - screened_rates(k_He4_Fe52_to_p_Co55)*Y(Fe52)*state.rho;
    jac.set(Fe52, He4, scratch);

    scratch = screened_rates(k_He4_Cr48_to_Fe52)*Y(He4)*state.rho;
    jac.set(Fe52, Cr48, scratch);

    scratch = screened_rates(k_p_Mn51_to_Fe52)*Y(H1)*state.rho;
    jac.set(Fe52, Mn51, scratch);

    scratch = -screened_rates(k_Fe52_to_He4_Cr48) - screened_rates(k_Fe52_to_p_Mn51) - screened_rates(k_He4_Fe52_to_Ni56)*Y(He4)*state.rho - screened_rates(k_He4_Fe52_to_p_Co55)*Y(He4)*state.rho;
    jac.set(Fe52, Fe52, scratch);

    scratch = screened_rates(k_p_Co55_to_He4_Fe52)*Y(H1)*state.rho;
    jac.set(Fe52, Co55, scratch);

    scratch = screened_rates(k_Ni56_to_He4_Fe52);
    jac.set(Fe52, Ni56, scratch);

    scratch = -screened_rates(k_p_Co55_to_He4_Fe52)*Y(Co55)*state.rho - screened_rates(k_p_Co55_to_Ni56)*Y(Co55)*state.rho;
    jac.set(Co55, H1, scratch);

    scratch = screened_rates(k_He4_Fe52_to_p_Co55)*Y(Fe52)*state.rho + screened_rates(k_He4_Mn51_to_Co55)*Y(Mn51)*state.rho;
    jac.set(Co55, He4, scratch);

    scratch = screened_rates(k_He4_Mn51_to_Co55)*Y(He4)*state.rho;
    jac.set(Co55, Mn51, scratch);

    scratch = screened_rates(k_He4_Fe52_to_p_Co55)*Y(He4)*state.rho;
    jac.set(Co55, Fe52, scratch);

    scratch = -screened_rates(k_Co55_to_He4_Mn51) - screened_rates(k_p_Co55_to_He4_Fe52)*Y(H1)*state.rho - screened_rates(k_p_Co55_to_Ni56)*Y(H1)*state.rho;
    jac.set(Co55, Co55, scratch);

    scratch = screened_rates(k_Ni56_to_p_Co55);
    jac.set(Co55, Ni56, scratch);

    scratch = screened_rates(k_p_Co55_to_Ni56)*Y(Co55)*state.rho;
    jac.set(Ni56, H1, scratch);

    scratch = screened_rates(k_He4_Fe52_to_Ni56)*Y(Fe52)*state.rho;
    jac.set(Ni56, He4, scratch);

    scratch = screened_rates(k_He4_Fe52_to_Ni56)*Y(He4)*state.rho;
    jac.set(Ni56, Fe52, scratch);

    scratch = screened_rates(k_p_Co55_to_Ni56)*Y(H1)*state.rho;
    jac.set(Ni56, Co55, scratch);

    scratch = -screened_rates(k_Ni56_to_He4_Fe52) - screened_rates(k_Ni56_to_p_Co55);
    jac.set(Ni56, Ni56, scratch);


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
