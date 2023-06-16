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

    if (disable_p_c12_to_n13) {
        rate_eval.screened_rates(k_p_c12_to_n13) = 0.0;
        if constexpr (std::is_same<T, rate_derivs_t>::value) {
            rate_eval.dscreened_rates_dT(k_p_c12_to_n13) = 0.0;
        }
        rate_eval.screened_rates(k_n13_to_p_c12_derived) = 0.0;
        if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_n13_to_p_c12_derived) = 0.0;
        }
    }

    if (disable_he4_n13_to_p_o16) {
        rate_eval.screened_rates(k_he4_n13_to_p_o16) = 0.0;
        if constexpr (std::is_same<T, rate_derivs_t>::value) {
            rate_eval.dscreened_rates_dT(k_he4_n13_to_p_o16) = 0.0;
        }
        rate_eval.screened_rates(k_p_o16_to_he4_n13_derived) = 0.0;
        if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_o16_to_he4_n13_derived) = 0.0;
        }
    }



    // Evaluate screening factors

    Real ratraw, dratraw_dT;
    Real scor, dscor_dt;
    Real scor2, dscor2_dt;


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
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 10.0_rt, 21.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_ne21_to_n_mg24);
    rate_eval.screened_rates(k_he4_ne21_to_n_mg24) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_ne21_to_n_mg24);
        rate_eval.dscreened_rates_dT(k_he4_ne21_to_n_mg24) = ratraw * dscor_dt + dratraw_dT * scor;
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


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 24.0_rt, 48.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_cr48_to_fe52_removed);
    rate_eval.screened_rates(k_he4_cr48_to_fe52_removed) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_cr48_to_fe52_removed);
        rate_eval.dscreened_rates_dT(k_he4_cr48_to_fe52_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_he4_cr48_to_p_mn51_removed);
    rate_eval.screened_rates(k_he4_cr48_to_p_mn51_removed) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_cr48_to_p_mn51_removed);
        rate_eval.dscreened_rates_dT(k_he4_cr48_to_p_mn51_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 25.0_rt, 51.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_mn51_to_fe52_removed);
    rate_eval.screened_rates(k_p_mn51_to_fe52_removed) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_mn51_to_fe52_removed);
        rate_eval.dscreened_rates_dT(k_p_mn51_to_fe52_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_mn51_to_he4_cr48_derived_removed);
    rate_eval.screened_rates(k_p_mn51_to_he4_cr48_derived_removed) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_mn51_to_he4_cr48_derived_removed);
        rate_eval.dscreened_rates_dT(k_p_mn51_to_he4_cr48_derived_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 26.0_rt, 52.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_fe52_to_ni56_removed);
    rate_eval.screened_rates(k_he4_fe52_to_ni56_removed) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_fe52_to_ni56_removed);
        rate_eval.dscreened_rates_dT(k_he4_fe52_to_ni56_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_he4_fe52_to_p_co55_removed);
    rate_eval.screened_rates(k_he4_fe52_to_p_co55_removed) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_fe52_to_p_co55_removed);
        rate_eval.dscreened_rates_dT(k_he4_fe52_to_p_co55_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 27.0_rt, 55.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_co55_to_ni56_removed);
    rate_eval.screened_rates(k_p_co55_to_ni56_removed) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_co55_to_ni56_removed);
        rate_eval.dscreened_rates_dT(k_p_co55_to_ni56_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_co55_to_he4_fe52_derived_removed);
    rate_eval.screened_rates(k_p_co55_to_he4_fe52_derived_removed) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_co55_to_he4_fe52_derived_removed);
        rate_eval.dscreened_rates_dT(k_p_co55_to_he4_fe52_derived_removed) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    // Fill approximate rates

    fill_approx_rates<do_T_derivatives, T>(tfactors, rate_eval);

    // Calculate tabular rates

    [[maybe_unused]] Real rate, drate_dt, edot_nu, edot_gamma;


}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rhs_nuc(const burn_t& state,
             Array1D<Real, 1, neqs>& ydot_nuc,
             const Array1D<Real, 1, NumSpec>& Y,
             const Array1D<Real, 1, NumRates>& screened_rates) {

    using namespace Rates;

    ydot_nuc(N) =
        (-screened_rates(k_n_ne20_to_ne21)*Y(N)*Y(Ne20)*state.rho + screened_rates(k_ne21_to_n_ne20_derived)*Y(Ne21)) +
        (-screened_rates(k_n_na22_to_na23)*Y(N)*Y(Na22)*state.rho + screened_rates(k_na23_to_n_na22_derived)*Y(Na23)) +
        (screened_rates(k_he4_ne21_to_n_mg24)*Y(He4)*Y(Ne21)*state.rho + -screened_rates(k_n_mg24_to_he4_ne21_derived)*Y(Mg24)*Y(N)*state.rho);

    ydot_nuc(H1) =
        (-screened_rates(k_p_c12_to_n13)*Y(C12)*Y(H1)*state.rho + screened_rates(k_n13_to_p_c12_derived)*Y(N13)) +
        (-screened_rates(k_p_ne21_to_na22)*Y(Ne21)*Y(H1)*state.rho + screened_rates(k_na22_to_p_ne21_derived)*Y(Na22)) +
        (-screened_rates(k_p_na23_to_mg24)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_mg24_to_p_na23_derived)*Y(Mg24)) +
        (-screened_rates(k_p_al27_to_si28)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_si28_to_p_al27_derived)*Y(Si28)) +
        (-screened_rates(k_p_p31_to_s32)*Y(P31)*Y(H1)*state.rho + screened_rates(k_s32_to_p_p31_derived)*Y(S32)) +
        0.5*screened_rates(k_c12_c12_to_p_na23)*std::pow(Y(C12), 2)*state.rho +
        (screened_rates(k_he4_n13_to_p_o16)*Y(He4)*Y(N13)*state.rho + -screened_rates(k_p_o16_to_he4_n13_derived)*Y(O16)*Y(H1)*state.rho) +
        screened_rates(k_c12_o16_to_p_al27)*Y(C12)*Y(O16)*state.rho +
        0.5*screened_rates(k_o16_o16_to_p_p31)*std::pow(Y(O16), 2)*state.rho +
        (screened_rates(k_he4_f18_to_p_ne21)*Y(F18)*Y(He4)*state.rho + -screened_rates(k_p_ne21_to_he4_f18_derived)*Y(Ne21)*Y(H1)*state.rho) +
        (-screened_rates(k_p_na23_to_he4_ne20)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_he4_ne20_to_p_na23_derived)*Y(He4)*Y(Ne20)*state.rho) +
        (-screened_rates(k_p_al27_to_he4_mg24)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_he4_mg24_to_p_al27_derived)*Y(He4)*Y(Mg24)*state.rho) +
        (-screened_rates(k_p_p31_to_he4_si28)*Y(P31)*Y(H1)*state.rho + screened_rates(k_he4_si28_to_p_p31_derived)*Y(He4)*Y(Si28)*state.rho);

    ydot_nuc(He4) =
        (-screened_rates(k_he4_c12_to_o16)*Y(C12)*Y(He4)*state.rho + screened_rates(k_o16_to_he4_c12_derived)*Y(O16)) +
        (-screened_rates(k_he4_n14_to_f18)*Y(He4)*Y(N14)*state.rho + screened_rates(k_f18_to_he4_n14_derived)*Y(F18)) +
        (-screened_rates(k_he4_o16_to_ne20)*Y(He4)*Y(O16)*state.rho + screened_rates(k_ne20_to_he4_o16_derived)*Y(Ne20)) +
        (-screened_rates(k_he4_f18_to_na22)*Y(F18)*Y(He4)*state.rho + screened_rates(k_na22_to_he4_f18_derived)*Y(Na22)) +
        (-screened_rates(k_he4_ne20_to_mg24)*Y(He4)*Y(Ne20)*state.rho + screened_rates(k_mg24_to_he4_ne20_derived)*Y(Mg24)) +
        (-screened_rates(k_he4_mg24_to_si28)*Y(He4)*Y(Mg24)*state.rho + screened_rates(k_si28_to_he4_mg24_derived)*Y(Si28)) +
        (-screened_rates(k_he4_si28_to_s32)*Y(He4)*Y(Si28)*state.rho + screened_rates(k_s32_to_he4_si28_derived)*Y(S32)) +
        0.5*screened_rates(k_c12_c12_to_he4_ne20)*std::pow(Y(C12), 2)*state.rho +
        (-screened_rates(k_he4_n13_to_p_o16)*Y(He4)*Y(N13)*state.rho + screened_rates(k_p_o16_to_he4_n13_derived)*Y(O16)*Y(H1)*state.rho) +
        screened_rates(k_c12_o16_to_he4_mg24)*Y(C12)*Y(O16)*state.rho +
        0.5*screened_rates(k_o16_o16_to_he4_si28)*std::pow(Y(O16), 2)*state.rho +
        (-screened_rates(k_he4_f18_to_p_ne21)*Y(F18)*Y(He4)*state.rho + screened_rates(k_p_ne21_to_he4_f18_derived)*Y(Ne21)*Y(H1)*state.rho) +
        (-screened_rates(k_he4_ne21_to_n_mg24)*Y(He4)*Y(Ne21)*state.rho + screened_rates(k_n_mg24_to_he4_ne21_derived)*Y(Mg24)*Y(N)*state.rho) +
        (screened_rates(k_p_na23_to_he4_ne20)*Y(Na23)*Y(H1)*state.rho + -screened_rates(k_he4_ne20_to_p_na23_derived)*Y(He4)*Y(Ne20)*state.rho) +
        (screened_rates(k_p_al27_to_he4_mg24)*Y(Al27)*Y(H1)*state.rho + -screened_rates(k_he4_mg24_to_p_al27_derived)*Y(He4)*Y(Mg24)*state.rho) +
        (screened_rates(k_p_p31_to_he4_si28)*Y(P31)*Y(H1)*state.rho + -screened_rates(k_he4_si28_to_p_p31_derived)*Y(He4)*Y(Si28)*state.rho) +
        (-0.5*screened_rates(k_he4_he4_he4_to_c12)*std::pow(Y(He4), 3)*std::pow(state.rho, 2) + 3.0*screened_rates(k_c12_to_he4_he4_he4_derived)*Y(C12)) +
        (-screened_rates(k_s32_he4_to_ar36_approx)*Y(He4)*Y(S32)*state.rho + screened_rates(k_ar36_to_s32_he4_approx)*Y(Ar36)) +
        (-screened_rates(k_ar36_he4_to_ca40_approx)*Y(Ar36)*Y(He4)*state.rho + screened_rates(k_ca40_to_ar36_he4_approx)*Y(Ca40)) +
        (-screened_rates(k_ca40_he4_to_ti44_approx)*Y(Ca40)*Y(He4)*state.rho + screened_rates(k_ti44_to_ca40_he4_approx)*Y(Ti44)) +
        (-screened_rates(k_ti44_he4_to_cr48_approx)*Y(He4)*Y(Ti44)*state.rho + screened_rates(k_cr48_to_ti44_he4_approx)*Y(Cr48)) +
        (-screened_rates(k_cr48_he4_to_fe52_approx)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_fe52_to_cr48_he4_approx)*Y(Fe52)) +
        (-screened_rates(k_fe52_he4_to_ni56_approx)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_ni56_to_fe52_he4_approx)*Y(Ni56));

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
        (-screened_rates(k_n_ne20_to_ne21)*Y(N)*Y(Ne20)*state.rho + screened_rates(k_ne21_to_n_ne20_derived)*Y(Ne21)) +
        (-screened_rates(k_he4_ne20_to_mg24)*Y(He4)*Y(Ne20)*state.rho + screened_rates(k_mg24_to_he4_ne20_derived)*Y(Mg24)) +
        0.5*screened_rates(k_c12_c12_to_he4_ne20)*std::pow(Y(C12), 2)*state.rho +
        (screened_rates(k_p_na23_to_he4_ne20)*Y(Na23)*Y(H1)*state.rho + -screened_rates(k_he4_ne20_to_p_na23_derived)*Y(He4)*Y(Ne20)*state.rho);

    ydot_nuc(Ne21) =
        (screened_rates(k_n_ne20_to_ne21)*Y(N)*Y(Ne20)*state.rho + -screened_rates(k_ne21_to_n_ne20_derived)*Y(Ne21)) +
        (-screened_rates(k_p_ne21_to_na22)*Y(Ne21)*Y(H1)*state.rho + screened_rates(k_na22_to_p_ne21_derived)*Y(Na22)) +
        (screened_rates(k_he4_f18_to_p_ne21)*Y(F18)*Y(He4)*state.rho + -screened_rates(k_p_ne21_to_he4_f18_derived)*Y(Ne21)*Y(H1)*state.rho) +
        (-screened_rates(k_he4_ne21_to_n_mg24)*Y(He4)*Y(Ne21)*state.rho + screened_rates(k_n_mg24_to_he4_ne21_derived)*Y(Mg24)*Y(N)*state.rho);

    ydot_nuc(Na22) =
        (screened_rates(k_he4_f18_to_na22)*Y(F18)*Y(He4)*state.rho + -screened_rates(k_na22_to_he4_f18_derived)*Y(Na22)) +
        (screened_rates(k_p_ne21_to_na22)*Y(Ne21)*Y(H1)*state.rho + -screened_rates(k_na22_to_p_ne21_derived)*Y(Na22)) +
        (-screened_rates(k_n_na22_to_na23)*Y(N)*Y(Na22)*state.rho + screened_rates(k_na23_to_n_na22_derived)*Y(Na23));

    ydot_nuc(Na23) =
        (screened_rates(k_n_na22_to_na23)*Y(N)*Y(Na22)*state.rho + -screened_rates(k_na23_to_n_na22_derived)*Y(Na23)) +
        (-screened_rates(k_p_na23_to_mg24)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_mg24_to_p_na23_derived)*Y(Mg24)) +
        0.5*screened_rates(k_c12_c12_to_p_na23)*std::pow(Y(C12), 2)*state.rho +
        (-screened_rates(k_p_na23_to_he4_ne20)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_he4_ne20_to_p_na23_derived)*Y(He4)*Y(Ne20)*state.rho);

    ydot_nuc(Mg24) =
        (screened_rates(k_he4_ne20_to_mg24)*Y(He4)*Y(Ne20)*state.rho + -screened_rates(k_mg24_to_he4_ne20_derived)*Y(Mg24)) +
        (screened_rates(k_p_na23_to_mg24)*Y(Na23)*Y(H1)*state.rho + -screened_rates(k_mg24_to_p_na23_derived)*Y(Mg24)) +
        (-screened_rates(k_he4_mg24_to_si28)*Y(He4)*Y(Mg24)*state.rho + screened_rates(k_si28_to_he4_mg24_derived)*Y(Si28)) +
        screened_rates(k_c12_o16_to_he4_mg24)*Y(C12)*Y(O16)*state.rho +
        (screened_rates(k_he4_ne21_to_n_mg24)*Y(He4)*Y(Ne21)*state.rho + -screened_rates(k_n_mg24_to_he4_ne21_derived)*Y(Mg24)*Y(N)*state.rho) +
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
        (screened_rates(k_ti44_he4_to_cr48_approx)*Y(He4)*Y(Ti44)*state.rho + -screened_rates(k_cr48_to_ti44_he4_approx)*Y(Cr48)) +
        (-screened_rates(k_cr48_he4_to_fe52_approx)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_fe52_to_cr48_he4_approx)*Y(Fe52));

    ydot_nuc(Fe52) =
        (screened_rates(k_cr48_he4_to_fe52_approx)*Y(Cr48)*Y(He4)*state.rho + -screened_rates(k_fe52_to_cr48_he4_approx)*Y(Fe52)) +
        (-screened_rates(k_fe52_he4_to_ni56_approx)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_ni56_to_fe52_he4_approx)*Y(Ni56));

    ydot_nuc(Ni56) =
        (screened_rates(k_fe52_he4_to_ni56_approx)*Y(Fe52)*Y(He4)*state.rho + -screened_rates(k_ni56_to_fe52_he4_approx)*Y(Ni56));

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

    // Get the thermal neutrino losses

    Real sneut, dsneutdt, dsneutdd, snuda, snudz;

    sneut5(state.T, state.rho, state.abar, state.zbar, sneut, dsneutdt, dsneutdd, snuda, snudz);

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

    scratch = -screened_rates(k_n_mg24_to_he4_ne21_derived)*Y(Mg24)*state.rho - screened_rates(k_n_na22_to_na23)*Y(Na22)*state.rho - screened_rates(k_n_ne20_to_ne21)*Y(Ne20)*state.rho;
    jac.set(N, N, scratch);

    scratch = screened_rates(k_he4_ne21_to_n_mg24)*Y(Ne21)*state.rho;
    jac.set(N, He4, scratch);

    scratch = -screened_rates(k_n_ne20_to_ne21)*Y(N)*state.rho;
    jac.set(N, Ne20, scratch);

    scratch = screened_rates(k_he4_ne21_to_n_mg24)*Y(He4)*state.rho + screened_rates(k_ne21_to_n_ne20_derived);
    jac.set(N, Ne21, scratch);

    scratch = -screened_rates(k_n_na22_to_na23)*Y(N)*state.rho;
    jac.set(N, Na22, scratch);

    scratch = screened_rates(k_na23_to_n_na22_derived);
    jac.set(N, Na23, scratch);

    scratch = -screened_rates(k_n_mg24_to_he4_ne21_derived)*Y(N)*state.rho;
    jac.set(N, Mg24, scratch);

    scratch = -screened_rates(k_p_al27_to_he4_mg24)*Y(Al27)*state.rho - screened_rates(k_p_al27_to_si28)*Y(Al27)*state.rho - screened_rates(k_p_c12_to_n13)*Y(C12)*state.rho - screened_rates(k_p_na23_to_he4_ne20)*Y(Na23)*state.rho - screened_rates(k_p_na23_to_mg24)*Y(Na23)*state.rho - screened_rates(k_p_ne21_to_he4_f18_derived)*Y(Ne21)*state.rho - screened_rates(k_p_ne21_to_na22)*Y(Ne21)*state.rho - screened_rates(k_p_o16_to_he4_n13_derived)*Y(O16)*state.rho - screened_rates(k_p_p31_to_he4_si28)*Y(P31)*state.rho - screened_rates(k_p_p31_to_s32)*Y(P31)*state.rho;
    jac.set(H1, H1, scratch);

    scratch = screened_rates(k_he4_f18_to_p_ne21)*Y(F18)*state.rho + screened_rates(k_he4_mg24_to_p_al27_derived)*Y(Mg24)*state.rho + screened_rates(k_he4_n13_to_p_o16)*Y(N13)*state.rho + screened_rates(k_he4_ne20_to_p_na23_derived)*Y(Ne20)*state.rho + screened_rates(k_he4_si28_to_p_p31_derived)*Y(Si28)*state.rho;
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

    scratch = screened_rates(k_n_mg24_to_he4_ne21_derived)*Y(Mg24)*state.rho;
    jac.set(He4, N, scratch);

    scratch = screened_rates(k_p_al27_to_he4_mg24)*Y(Al27)*state.rho + screened_rates(k_p_na23_to_he4_ne20)*Y(Na23)*state.rho + screened_rates(k_p_ne21_to_he4_f18_derived)*Y(Ne21)*state.rho + screened_rates(k_p_o16_to_he4_n13_derived)*Y(O16)*state.rho + screened_rates(k_p_p31_to_he4_si28)*Y(P31)*state.rho;
    jac.set(He4, H1, scratch);

    scratch = -screened_rates(k_ar36_he4_to_ca40_approx)*Y(Ar36)*state.rho - screened_rates(k_ca40_he4_to_ti44_approx)*Y(Ca40)*state.rho - screened_rates(k_cr48_he4_to_fe52_approx)*Y(Cr48)*state.rho - screened_rates(k_fe52_he4_to_ni56_approx)*Y(Fe52)*state.rho - screened_rates(k_he4_c12_to_o16)*Y(C12)*state.rho - screened_rates(k_he4_f18_to_na22)*Y(F18)*state.rho - screened_rates(k_he4_f18_to_p_ne21)*Y(F18)*state.rho - 1.5*screened_rates(k_he4_he4_he4_to_c12)*std::pow(Y(He4), 2)*std::pow(state.rho, 2) - screened_rates(k_he4_mg24_to_p_al27_derived)*Y(Mg24)*state.rho - screened_rates(k_he4_mg24_to_si28)*Y(Mg24)*state.rho - screened_rates(k_he4_n13_to_p_o16)*Y(N13)*state.rho - screened_rates(k_he4_n14_to_f18)*Y(N14)*state.rho - screened_rates(k_he4_ne20_to_mg24)*Y(Ne20)*state.rho - screened_rates(k_he4_ne20_to_p_na23_derived)*Y(Ne20)*state.rho - screened_rates(k_he4_ne21_to_n_mg24)*Y(Ne21)*state.rho - screened_rates(k_he4_o16_to_ne20)*Y(O16)*state.rho - screened_rates(k_he4_si28_to_p_p31_derived)*Y(Si28)*state.rho - screened_rates(k_he4_si28_to_s32)*Y(Si28)*state.rho - screened_rates(k_s32_he4_to_ar36_approx)*Y(S32)*state.rho - screened_rates(k_ti44_he4_to_cr48_approx)*Y(Ti44)*state.rho;
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

    scratch = -screened_rates(k_he4_ne21_to_n_mg24)*Y(He4)*state.rho + screened_rates(k_p_ne21_to_he4_f18_derived)*Y(H1)*state.rho;
    jac.set(He4, Ne21, scratch);

    scratch = screened_rates(k_na22_to_he4_f18_derived);
    jac.set(He4, Na22, scratch);

    scratch = screened_rates(k_p_na23_to_he4_ne20)*Y(H1)*state.rho;
    jac.set(He4, Na23, scratch);

    scratch = -screened_rates(k_he4_mg24_to_p_al27_derived)*Y(He4)*state.rho - screened_rates(k_he4_mg24_to_si28)*Y(He4)*state.rho + screened_rates(k_mg24_to_he4_ne20_derived) + screened_rates(k_n_mg24_to_he4_ne21_derived)*Y(N)*state.rho;
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

    scratch = -screened_rates(k_cr48_he4_to_fe52_approx)*Y(He4)*state.rho + screened_rates(k_cr48_to_ti44_he4_approx);
    jac.set(He4, Cr48, scratch);

    scratch = -screened_rates(k_fe52_he4_to_ni56_approx)*Y(He4)*state.rho + screened_rates(k_fe52_to_cr48_he4_approx);
    jac.set(He4, Fe52, scratch);

    scratch = screened_rates(k_ni56_to_fe52_he4_approx);
    jac.set(He4, Ni56, scratch);

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

    scratch = -screened_rates(k_n_ne20_to_ne21)*Y(Ne20)*state.rho;
    jac.set(Ne20, N, scratch);

    scratch = screened_rates(k_p_na23_to_he4_ne20)*Y(Na23)*state.rho;
    jac.set(Ne20, H1, scratch);

    scratch = -screened_rates(k_he4_ne20_to_mg24)*Y(Ne20)*state.rho - screened_rates(k_he4_ne20_to_p_na23_derived)*Y(Ne20)*state.rho + screened_rates(k_he4_o16_to_ne20)*Y(O16)*state.rho;
    jac.set(Ne20, He4, scratch);

    scratch = 1.0*screened_rates(k_c12_c12_to_he4_ne20)*Y(C12)*state.rho;
    jac.set(Ne20, C12, scratch);

    scratch = screened_rates(k_he4_o16_to_ne20)*Y(He4)*state.rho;
    jac.set(Ne20, O16, scratch);

    scratch = -screened_rates(k_he4_ne20_to_mg24)*Y(He4)*state.rho - screened_rates(k_he4_ne20_to_p_na23_derived)*Y(He4)*state.rho - screened_rates(k_n_ne20_to_ne21)*Y(N)*state.rho - screened_rates(k_ne20_to_he4_o16_derived);
    jac.set(Ne20, Ne20, scratch);

    scratch = screened_rates(k_ne21_to_n_ne20_derived);
    jac.set(Ne20, Ne21, scratch);

    scratch = screened_rates(k_p_na23_to_he4_ne20)*Y(H1)*state.rho;
    jac.set(Ne20, Na23, scratch);

    scratch = screened_rates(k_mg24_to_he4_ne20_derived);
    jac.set(Ne20, Mg24, scratch);

    scratch = screened_rates(k_n_mg24_to_he4_ne21_derived)*Y(Mg24)*state.rho + screened_rates(k_n_ne20_to_ne21)*Y(Ne20)*state.rho;
    jac.set(Ne21, N, scratch);

    scratch = -screened_rates(k_p_ne21_to_he4_f18_derived)*Y(Ne21)*state.rho - screened_rates(k_p_ne21_to_na22)*Y(Ne21)*state.rho;
    jac.set(Ne21, H1, scratch);

    scratch = screened_rates(k_he4_f18_to_p_ne21)*Y(F18)*state.rho - screened_rates(k_he4_ne21_to_n_mg24)*Y(Ne21)*state.rho;
    jac.set(Ne21, He4, scratch);

    scratch = screened_rates(k_he4_f18_to_p_ne21)*Y(He4)*state.rho;
    jac.set(Ne21, F18, scratch);

    scratch = screened_rates(k_n_ne20_to_ne21)*Y(N)*state.rho;
    jac.set(Ne21, Ne20, scratch);

    scratch = -screened_rates(k_he4_ne21_to_n_mg24)*Y(He4)*state.rho - screened_rates(k_ne21_to_n_ne20_derived) - screened_rates(k_p_ne21_to_he4_f18_derived)*Y(H1)*state.rho - screened_rates(k_p_ne21_to_na22)*Y(H1)*state.rho;
    jac.set(Ne21, Ne21, scratch);

    scratch = screened_rates(k_na22_to_p_ne21_derived);
    jac.set(Ne21, Na22, scratch);

    scratch = screened_rates(k_n_mg24_to_he4_ne21_derived)*Y(N)*state.rho;
    jac.set(Ne21, Mg24, scratch);

    scratch = -screened_rates(k_n_na22_to_na23)*Y(Na22)*state.rho;
    jac.set(Na22, N, scratch);

    scratch = screened_rates(k_p_ne21_to_na22)*Y(Ne21)*state.rho;
    jac.set(Na22, H1, scratch);

    scratch = screened_rates(k_he4_f18_to_na22)*Y(F18)*state.rho;
    jac.set(Na22, He4, scratch);

    scratch = screened_rates(k_he4_f18_to_na22)*Y(He4)*state.rho;
    jac.set(Na22, F18, scratch);

    scratch = screened_rates(k_p_ne21_to_na22)*Y(H1)*state.rho;
    jac.set(Na22, Ne21, scratch);

    scratch = -screened_rates(k_n_na22_to_na23)*Y(N)*state.rho - screened_rates(k_na22_to_he4_f18_derived) - screened_rates(k_na22_to_p_ne21_derived);
    jac.set(Na22, Na22, scratch);

    scratch = screened_rates(k_na23_to_n_na22_derived);
    jac.set(Na22, Na23, scratch);

    scratch = screened_rates(k_n_na22_to_na23)*Y(Na22)*state.rho;
    jac.set(Na23, N, scratch);

    scratch = -screened_rates(k_p_na23_to_he4_ne20)*Y(Na23)*state.rho - screened_rates(k_p_na23_to_mg24)*Y(Na23)*state.rho;
    jac.set(Na23, H1, scratch);

    scratch = screened_rates(k_he4_ne20_to_p_na23_derived)*Y(Ne20)*state.rho;
    jac.set(Na23, He4, scratch);

    scratch = 1.0*screened_rates(k_c12_c12_to_p_na23)*Y(C12)*state.rho;
    jac.set(Na23, C12, scratch);

    scratch = screened_rates(k_he4_ne20_to_p_na23_derived)*Y(He4)*state.rho;
    jac.set(Na23, Ne20, scratch);

    scratch = screened_rates(k_n_na22_to_na23)*Y(N)*state.rho;
    jac.set(Na23, Na22, scratch);

    scratch = -screened_rates(k_na23_to_n_na22_derived) - screened_rates(k_p_na23_to_he4_ne20)*Y(H1)*state.rho - screened_rates(k_p_na23_to_mg24)*Y(H1)*state.rho;
    jac.set(Na23, Na23, scratch);

    scratch = screened_rates(k_mg24_to_p_na23_derived);
    jac.set(Na23, Mg24, scratch);

    scratch = -screened_rates(k_n_mg24_to_he4_ne21_derived)*Y(Mg24)*state.rho;
    jac.set(Mg24, N, scratch);

    scratch = screened_rates(k_p_al27_to_he4_mg24)*Y(Al27)*state.rho + screened_rates(k_p_na23_to_mg24)*Y(Na23)*state.rho;
    jac.set(Mg24, H1, scratch);

    scratch = -screened_rates(k_he4_mg24_to_p_al27_derived)*Y(Mg24)*state.rho - screened_rates(k_he4_mg24_to_si28)*Y(Mg24)*state.rho + screened_rates(k_he4_ne20_to_mg24)*Y(Ne20)*state.rho + screened_rates(k_he4_ne21_to_n_mg24)*Y(Ne21)*state.rho;
    jac.set(Mg24, He4, scratch);

    scratch = 1.0*screened_rates(k_c12_c12_to_mg24_modified)*Y(C12)*state.rho + screened_rates(k_c12_o16_to_he4_mg24)*Y(O16)*state.rho;
    jac.set(Mg24, C12, scratch);

    scratch = screened_rates(k_c12_o16_to_he4_mg24)*Y(C12)*state.rho;
    jac.set(Mg24, O16, scratch);

    scratch = screened_rates(k_he4_ne20_to_mg24)*Y(He4)*state.rho;
    jac.set(Mg24, Ne20, scratch);

    scratch = screened_rates(k_he4_ne21_to_n_mg24)*Y(He4)*state.rho;
    jac.set(Mg24, Ne21, scratch);

    scratch = screened_rates(k_p_na23_to_mg24)*Y(H1)*state.rho;
    jac.set(Mg24, Na23, scratch);

    scratch = -screened_rates(k_he4_mg24_to_p_al27_derived)*Y(He4)*state.rho - screened_rates(k_he4_mg24_to_si28)*Y(He4)*state.rho - screened_rates(k_mg24_to_he4_ne20_derived) - screened_rates(k_mg24_to_p_na23_derived) - screened_rates(k_n_mg24_to_he4_ne21_derived)*Y(N)*state.rho;
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

    scratch = -screened_rates(k_cr48_he4_to_fe52_approx)*Y(Cr48)*state.rho + screened_rates(k_ti44_he4_to_cr48_approx)*Y(Ti44)*state.rho;
    jac.set(Cr48, He4, scratch);

    scratch = screened_rates(k_ti44_he4_to_cr48_approx)*Y(He4)*state.rho;
    jac.set(Cr48, Ti44, scratch);

    scratch = -screened_rates(k_cr48_he4_to_fe52_approx)*Y(He4)*state.rho - screened_rates(k_cr48_to_ti44_he4_approx);
    jac.set(Cr48, Cr48, scratch);

    scratch = screened_rates(k_fe52_to_cr48_he4_approx);
    jac.set(Cr48, Fe52, scratch);

    scratch = screened_rates(k_cr48_he4_to_fe52_approx)*Y(Cr48)*state.rho - screened_rates(k_fe52_he4_to_ni56_approx)*Y(Fe52)*state.rho;
    jac.set(Fe52, He4, scratch);

    scratch = screened_rates(k_cr48_he4_to_fe52_approx)*Y(He4)*state.rho;
    jac.set(Fe52, Cr48, scratch);

    scratch = -screened_rates(k_fe52_he4_to_ni56_approx)*Y(He4)*state.rho - screened_rates(k_fe52_to_cr48_he4_approx);
    jac.set(Fe52, Fe52, scratch);

    scratch = screened_rates(k_ni56_to_fe52_he4_approx);
    jac.set(Fe52, Ni56, scratch);

    scratch = screened_rates(k_fe52_he4_to_ni56_approx)*Y(Fe52)*state.rho;
    jac.set(Ni56, He4, scratch);

    scratch = screened_rates(k_fe52_he4_to_ni56_approx)*Y(He4)*state.rho;
    jac.set(Ni56, Fe52, scratch);

    scratch = -screened_rates(k_ni56_to_fe52_he4_approx);
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

    Real sneut, dsneutdt, dsneutdd, snuda, snudz;
    sneut5(state.T, state.rho, state.abar, state.zbar, sneut, dsneutdt, dsneutdd, snuda, snudz);

    for (int j = 1; j <= NumSpec; ++j) {
       Real b1 = (-state.abar * state.abar * snuda + (zion[j-1] - state.zbar) * state.abar * snudz);
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
