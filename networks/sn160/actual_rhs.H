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

    for (int n = 1; n <= NumSpec; ++n) {
        enuc += dydt(n) * network::mion(n);
    }

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


#ifdef SCREENING
    // Evaluate screening factors

    amrex::Real ratraw, dratraw_dT;
    amrex::Real scor, dscor_dt;
    [[maybe_unused]] amrex::Real scor2, dscor2_dt;


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 1.0_rt, 1.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_p_to_d_weak_bet_pos_);
    rate_eval.screened_rates(k_p_p_to_d_weak_bet_pos_) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_p_to_d_weak_bet_pos_);
        rate_eval.dscreened_rates_dT(k_p_p_to_d_weak_bet_pos_) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_p_to_d_weak_electron_capture);
    rate_eval.screened_rates(k_p_p_to_d_weak_electron_capture) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_p_to_d_weak_electron_capture);
        rate_eval.dscreened_rates_dT(k_p_p_to_d_weak_electron_capture) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_n_p_p_to_p_d);
    rate_eval.screened_rates(k_n_p_p_to_p_d) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_n_p_p_to_p_d);
        rate_eval.dscreened_rates_dT(k_n_p_p_to_p_d) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_p_He4_He4_to_He3_Be7);
    rate_eval.screened_rates(k_p_p_He4_He4_to_He3_Be7) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_p_He4_He4_to_He3_Be7);
        rate_eval.dscreened_rates_dT(k_p_p_He4_He4_to_He3_Be7) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 1.0_rt, 2.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_d_to_He3);
    rate_eval.screened_rates(k_p_d_to_He3) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_d_to_He3);
        rate_eval.dscreened_rates_dT(k_p_d_to_He3) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_d_to_n_p_p);
    rate_eval.screened_rates(k_p_d_to_n_p_p) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_d_to_n_p_p);
        rate_eval.dscreened_rates_dT(k_p_d_to_n_p_p) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 2.0_rt, 1.0_rt, 2.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_d_d_to_He4);
    rate_eval.screened_rates(k_d_d_to_He4) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_d_d_to_He4);
        rate_eval.dscreened_rates_dT(k_d_d_to_He4) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_d_d_to_n_He3);
    rate_eval.screened_rates(k_d_d_to_n_He3) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_d_d_to_n_He3);
        rate_eval.dscreened_rates_dT(k_d_d_to_n_He3) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 2.0_rt, 2.0_rt, 4.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_d_to_Li6);
    rate_eval.screened_rates(k_He4_d_to_Li6) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_d_to_Li6);
        rate_eval.dscreened_rates_dT(k_He4_d_to_Li6) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 2.0_rt, 3.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_He3_to_He4_weak_bet_pos_);
    rate_eval.screened_rates(k_p_He3_to_He4_weak_bet_pos_) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_He3_to_He4_weak_bet_pos_);
        rate_eval.dscreened_rates_dT(k_p_He3_to_He4_weak_bet_pos_) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 2.0_rt, 3.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_He3_to_Be7);
    rate_eval.screened_rates(k_He4_He3_to_Be7) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_He3_to_Be7);
        rate_eval.dscreened_rates_dT(k_He4_He3_to_Be7) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_He3_to_p_Li6);
    rate_eval.screened_rates(k_He4_He3_to_p_Li6) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_He3_to_p_Li6);
        rate_eval.dscreened_rates_dT(k_He4_He3_to_p_Li6) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 3.0_rt, 6.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Li6_to_Be7);
    rate_eval.screened_rates(k_p_Li6_to_Be7) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Li6_to_Be7);
        rate_eval.dscreened_rates_dT(k_p_Li6_to_Be7) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Li6_to_He4_He3);
    rate_eval.screened_rates(k_p_Li6_to_He4_He3) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Li6_to_He4_He3);
        rate_eval.dscreened_rates_dT(k_p_Li6_to_He4_He3) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 3.0_rt, 6.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Li6_to_B10);
    rate_eval.screened_rates(k_He4_Li6_to_B10) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Li6_to_B10);
        rate_eval.dscreened_rates_dT(k_He4_Li6_to_B10) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Li6_to_p_Be9);
    rate_eval.screened_rates(k_He4_Li6_to_p_Be9) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Li6_to_p_Be9);
        rate_eval.dscreened_rates_dT(k_He4_Li6_to_p_Be9) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 3.0_rt, 7.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Li7_to_B11);
    rate_eval.screened_rates(k_He4_Li7_to_B11) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Li7_to_B11);
        rate_eval.dscreened_rates_dT(k_He4_Li7_to_B11) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Li7_to_n_B10);
    rate_eval.screened_rates(k_He4_Li7_to_n_B10) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Li7_to_n_B10);
        rate_eval.dscreened_rates_dT(k_He4_Li7_to_n_B10) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 4.0_rt, 7.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Be7_to_B8);
    rate_eval.screened_rates(k_p_Be7_to_B8) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Be7_to_B8);
        rate_eval.dscreened_rates_dT(k_p_Be7_to_B8) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 4.0_rt, 9.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Be9_to_B10);
    rate_eval.screened_rates(k_p_Be9_to_B10) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Be9_to_B10);
        rate_eval.dscreened_rates_dT(k_p_Be9_to_B10) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Be9_to_He4_Li6);
    rate_eval.screened_rates(k_p_Be9_to_He4_Li6) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Be9_to_He4_Li6);
        rate_eval.dscreened_rates_dT(k_p_Be9_to_He4_Li6) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Be9_to_d_He4_He4);
    rate_eval.screened_rates(k_p_Be9_to_d_He4_He4) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Be9_to_d_He4_He4);
        rate_eval.dscreened_rates_dT(k_p_Be9_to_d_He4_He4) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Be9_to_n_p_He4_He4);
    rate_eval.screened_rates(k_p_Be9_to_n_p_He4_He4) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Be9_to_n_p_He4_He4);
        rate_eval.dscreened_rates_dT(k_p_Be9_to_n_p_He4_He4) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 5.0_rt, 11.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_B11_to_C12);
    rate_eval.screened_rates(k_p_B11_to_C12) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_B11_to_C12);
        rate_eval.dscreened_rates_dT(k_p_B11_to_C12) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_B11_to_He4_He4_He4);
    rate_eval.screened_rates(k_p_B11_to_He4_He4_He4) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_B11_to_He4_He4_He4);
        rate_eval.dscreened_rates_dT(k_p_B11_to_He4_He4_He4) = ratraw * dscor_dt + dratraw_dT * scor;
    }


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

    ratraw = rate_eval.screened_rates(k_He4_C12_to_n_O15);
    rate_eval.screened_rates(k_He4_C12_to_n_O15) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_C12_to_n_O15);
        rate_eval.dscreened_rates_dT(k_He4_C12_to_n_O15) = ratraw * dscor_dt + dratraw_dT * scor;
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


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_C13_to_N14);
    rate_eval.screened_rates(k_p_C13_to_N14) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_C13_to_N14);
        rate_eval.dscreened_rates_dT(k_p_C13_to_N14) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_C13_to_n_N13);
    rate_eval.screened_rates(k_p_C13_to_n_N13) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_C13_to_n_N13);
        rate_eval.dscreened_rates_dT(k_p_C13_to_n_N13) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_C13_to_He4_B10);
    rate_eval.screened_rates(k_p_C13_to_He4_B10) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_C13_to_He4_B10);
        rate_eval.dscreened_rates_dT(k_p_C13_to_He4_B10) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 6.0_rt, 14.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_C14_to_N15);
    rate_eval.screened_rates(k_p_C14_to_N15) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_C14_to_N15);
        rate_eval.dscreened_rates_dT(k_p_C14_to_N15) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_C14_to_n_N14);
    rate_eval.screened_rates(k_p_C14_to_n_N14) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_C14_to_n_N14);
        rate_eval.dscreened_rates_dT(k_p_C14_to_n_N14) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_C14_to_He4_B11);
    rate_eval.screened_rates(k_p_C14_to_He4_B11) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_C14_to_He4_B11);
        rate_eval.dscreened_rates_dT(k_p_C14_to_He4_B11) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 6.0_rt, 14.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_C14_to_O18);
    rate_eval.screened_rates(k_He4_C14_to_O18) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_C14_to_O18);
        rate_eval.dscreened_rates_dT(k_He4_C14_to_O18) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_C14_to_n_O17);
    rate_eval.screened_rates(k_He4_C14_to_n_O17) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_C14_to_n_O17);
        rate_eval.dscreened_rates_dT(k_He4_C14_to_n_O17) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 7.0_rt, 13.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
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


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_N14_to_O15);
    rate_eval.screened_rates(k_p_N14_to_O15) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_N14_to_O15);
        rate_eval.dscreened_rates_dT(k_p_N14_to_O15) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_N14_to_n_O14);
    rate_eval.screened_rates(k_p_N14_to_n_O14) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_N14_to_n_O14);
        rate_eval.dscreened_rates_dT(k_p_N14_to_n_O14) = ratraw * dscor_dt + dratraw_dT * scor;
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

    ratraw = rate_eval.screened_rates(k_He4_N14_to_n_F17);
    rate_eval.screened_rates(k_He4_N14_to_n_F17) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_N14_to_n_F17);
        rate_eval.dscreened_rates_dT(k_He4_N14_to_n_F17) = ratraw * dscor_dt + dratraw_dT * scor;
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


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_N15_to_O16);
    rate_eval.screened_rates(k_p_N15_to_O16) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_N15_to_O16);
        rate_eval.dscreened_rates_dT(k_p_N15_to_O16) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_N15_to_n_O15);
    rate_eval.screened_rates(k_p_N15_to_n_O15) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_N15_to_n_O15);
        rate_eval.dscreened_rates_dT(k_p_N15_to_n_O15) = ratraw * dscor_dt + dratraw_dT * scor;
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


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_N15_to_F19);
    rate_eval.screened_rates(k_He4_N15_to_F19) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_N15_to_F19);
        rate_eval.dscreened_rates_dT(k_He4_N15_to_F19) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_N15_to_n_F18);
    rate_eval.screened_rates(k_He4_N15_to_n_F18) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_N15_to_n_F18);
        rate_eval.dscreened_rates_dT(k_He4_N15_to_n_F18) = ratraw * dscor_dt + dratraw_dT * scor;
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


        actual_screen(pstate, scn_fac, scor, dscor_dt);
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


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_O15_to_Ne19);
    rate_eval.screened_rates(k_He4_O15_to_Ne19) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_O15_to_Ne19);
        rate_eval.dscreened_rates_dT(k_He4_O15_to_Ne19) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_O15_to_n_Ne18);
    rate_eval.screened_rates(k_He4_O15_to_n_Ne18) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_O15_to_n_Ne18);
        rate_eval.dscreened_rates_dT(k_He4_O15_to_n_Ne18) = ratraw * dscor_dt + dratraw_dT * scor;
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


        actual_screen(pstate, scn_fac, scor, dscor_dt);
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


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_O16_to_Ne20);
    rate_eval.screened_rates(k_He4_O16_to_Ne20) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_O16_to_Ne20);
        rate_eval.dscreened_rates_dT(k_He4_O16_to_Ne20) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_O16_to_n_Ne19);
    rate_eval.screened_rates(k_He4_O16_to_n_Ne19) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_O16_to_n_Ne19);
        rate_eval.dscreened_rates_dT(k_He4_O16_to_n_Ne19) = ratraw * dscor_dt + dratraw_dT * scor;
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


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_O17_to_F18);
    rate_eval.screened_rates(k_p_O17_to_F18) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_O17_to_F18);
        rate_eval.dscreened_rates_dT(k_p_O17_to_F18) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_O17_to_n_F17);
    rate_eval.screened_rates(k_p_O17_to_n_F17) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_O17_to_n_F17);
        rate_eval.dscreened_rates_dT(k_p_O17_to_n_F17) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_O17_to_He4_N14);
    rate_eval.screened_rates(k_p_O17_to_He4_N14) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_O17_to_He4_N14);
        rate_eval.dscreened_rates_dT(k_p_O17_to_He4_N14) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 8.0_rt, 17.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_O17_to_Ne21);
    rate_eval.screened_rates(k_He4_O17_to_Ne21) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_O17_to_Ne21);
        rate_eval.dscreened_rates_dT(k_He4_O17_to_Ne21) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_O17_to_n_Ne20);
    rate_eval.screened_rates(k_He4_O17_to_n_Ne20) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_O17_to_n_Ne20);
        rate_eval.dscreened_rates_dT(k_He4_O17_to_n_Ne20) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 8.0_rt, 18.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_O18_to_F19);
    rate_eval.screened_rates(k_p_O18_to_F19) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_O18_to_F19);
        rate_eval.dscreened_rates_dT(k_p_O18_to_F19) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_O18_to_n_F18);
    rate_eval.screened_rates(k_p_O18_to_n_F18) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_O18_to_n_F18);
        rate_eval.dscreened_rates_dT(k_p_O18_to_n_F18) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_O18_to_He4_N15);
    rate_eval.screened_rates(k_p_O18_to_He4_N15) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_O18_to_He4_N15);
        rate_eval.dscreened_rates_dT(k_p_O18_to_He4_N15) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 8.0_rt, 18.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_O18_to_Ne22);
    rate_eval.screened_rates(k_He4_O18_to_Ne22) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_O18_to_Ne22);
        rate_eval.dscreened_rates_dT(k_He4_O18_to_Ne22) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_O18_to_n_Ne21);
    rate_eval.screened_rates(k_He4_O18_to_n_Ne21) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_O18_to_n_Ne21);
        rate_eval.dscreened_rates_dT(k_He4_O18_to_n_Ne21) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 9.0_rt, 17.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
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
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 9.0_rt, 17.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_F17_to_Na21);
    rate_eval.screened_rates(k_He4_F17_to_Na21) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_F17_to_Na21);
        rate_eval.dscreened_rates_dT(k_He4_F17_to_Na21) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_F17_to_p_Ne20);
    rate_eval.screened_rates(k_He4_F17_to_p_Ne20) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_F17_to_p_Ne20);
        rate_eval.dscreened_rates_dT(k_He4_F17_to_p_Ne20) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 9.0_rt, 18.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_F18_to_Ne19);
    rate_eval.screened_rates(k_p_F18_to_Ne19) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_F18_to_Ne19);
        rate_eval.dscreened_rates_dT(k_p_F18_to_Ne19) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_F18_to_n_Ne18);
    rate_eval.screened_rates(k_p_F18_to_n_Ne18) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_F18_to_n_Ne18);
        rate_eval.dscreened_rates_dT(k_p_F18_to_n_Ne18) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_F18_to_He4_O15);
    rate_eval.screened_rates(k_p_F18_to_He4_O15) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_F18_to_He4_O15);
        rate_eval.dscreened_rates_dT(k_p_F18_to_He4_O15) = ratraw * dscor_dt + dratraw_dT * scor;
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

    ratraw = rate_eval.screened_rates(k_He4_F18_to_n_Na21);
    rate_eval.screened_rates(k_He4_F18_to_n_Na21) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_F18_to_n_Na21);
        rate_eval.dscreened_rates_dT(k_He4_F18_to_n_Na21) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_F18_to_p_Ne21);
    rate_eval.screened_rates(k_He4_F18_to_p_Ne21) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_F18_to_p_Ne21);
        rate_eval.dscreened_rates_dT(k_He4_F18_to_p_Ne21) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 9.0_rt, 19.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_F19_to_Ne20);
    rate_eval.screened_rates(k_p_F19_to_Ne20) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_F19_to_Ne20);
        rate_eval.dscreened_rates_dT(k_p_F19_to_Ne20) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_F19_to_n_Ne19);
    rate_eval.screened_rates(k_p_F19_to_n_Ne19) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_F19_to_n_Ne19);
        rate_eval.dscreened_rates_dT(k_p_F19_to_n_Ne19) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_F19_to_He4_O16);
    rate_eval.screened_rates(k_p_F19_to_He4_O16) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_F19_to_He4_O16);
        rate_eval.dscreened_rates_dT(k_p_F19_to_He4_O16) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 9.0_rt, 19.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_F19_to_Na23);
    rate_eval.screened_rates(k_He4_F19_to_Na23) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_F19_to_Na23);
        rate_eval.dscreened_rates_dT(k_He4_F19_to_Na23) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_F19_to_n_Na22);
    rate_eval.screened_rates(k_He4_F19_to_n_Na22) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_F19_to_n_Na22);
        rate_eval.dscreened_rates_dT(k_He4_F19_to_n_Na22) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_F19_to_p_Ne22);
    rate_eval.screened_rates(k_He4_F19_to_p_Ne22) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_F19_to_p_Ne22);
        rate_eval.dscreened_rates_dT(k_He4_F19_to_p_Ne22) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 10.0_rt, 19.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ne19_to_Mg23);
    rate_eval.screened_rates(k_He4_Ne19_to_Mg23) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ne19_to_Mg23);
        rate_eval.dscreened_rates_dT(k_He4_Ne19_to_Mg23) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ne19_to_p_Na22);
    rate_eval.screened_rates(k_He4_Ne19_to_p_Na22) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ne19_to_p_Na22);
        rate_eval.dscreened_rates_dT(k_He4_Ne19_to_p_Na22) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 10.0_rt, 20.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ne20_to_Na21);
    rate_eval.screened_rates(k_p_Ne20_to_Na21) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ne20_to_Na21);
        rate_eval.dscreened_rates_dT(k_p_Ne20_to_Na21) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ne20_to_He4_F17);
    rate_eval.screened_rates(k_p_Ne20_to_He4_F17) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ne20_to_He4_F17);
        rate_eval.dscreened_rates_dT(k_p_Ne20_to_He4_F17) = ratraw * dscor_dt + dratraw_dT * scor;
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

    ratraw = rate_eval.screened_rates(k_He4_Ne20_to_n_Mg23);
    rate_eval.screened_rates(k_He4_Ne20_to_n_Mg23) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ne20_to_n_Mg23);
        rate_eval.dscreened_rates_dT(k_He4_Ne20_to_n_Mg23) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ne20_to_p_Na23);
    rate_eval.screened_rates(k_He4_Ne20_to_p_Na23) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ne20_to_p_Na23);
        rate_eval.dscreened_rates_dT(k_He4_Ne20_to_p_Na23) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ne20_to_C12_C12);
    rate_eval.screened_rates(k_He4_Ne20_to_C12_C12) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ne20_to_C12_C12);
        rate_eval.dscreened_rates_dT(k_He4_Ne20_to_C12_C12) = ratraw * dscor_dt + dratraw_dT * scor;
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

    ratraw = rate_eval.screened_rates(k_p_Ne21_to_n_Na21);
    rate_eval.screened_rates(k_p_Ne21_to_n_Na21) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ne21_to_n_Na21);
        rate_eval.dscreened_rates_dT(k_p_Ne21_to_n_Na21) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ne21_to_He4_F18);
    rate_eval.screened_rates(k_p_Ne21_to_He4_F18) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ne21_to_He4_F18);
        rate_eval.dscreened_rates_dT(k_p_Ne21_to_He4_F18) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 10.0_rt, 21.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ne21_to_Mg25);
    rate_eval.screened_rates(k_He4_Ne21_to_Mg25) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ne21_to_Mg25);
        rate_eval.dscreened_rates_dT(k_He4_Ne21_to_Mg25) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ne21_to_n_Mg24);
    rate_eval.screened_rates(k_He4_Ne21_to_n_Mg24) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ne21_to_n_Mg24);
        rate_eval.dscreened_rates_dT(k_He4_Ne21_to_n_Mg24) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 10.0_rt, 22.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ne22_to_Na23);
    rate_eval.screened_rates(k_p_Ne22_to_Na23) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ne22_to_Na23);
        rate_eval.dscreened_rates_dT(k_p_Ne22_to_Na23) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ne22_to_n_Na22);
    rate_eval.screened_rates(k_p_Ne22_to_n_Na22) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ne22_to_n_Na22);
        rate_eval.dscreened_rates_dT(k_p_Ne22_to_n_Na22) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ne22_to_He4_F19);
    rate_eval.screened_rates(k_p_Ne22_to_He4_F19) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ne22_to_He4_F19);
        rate_eval.dscreened_rates_dT(k_p_Ne22_to_He4_F19) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 10.0_rt, 22.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ne22_to_Mg26);
    rate_eval.screened_rates(k_He4_Ne22_to_Mg26) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ne22_to_Mg26);
        rate_eval.dscreened_rates_dT(k_He4_Ne22_to_Mg26) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ne22_to_n_Mg25);
    rate_eval.screened_rates(k_He4_Ne22_to_n_Mg25) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ne22_to_n_Mg25);
        rate_eval.dscreened_rates_dT(k_He4_Ne22_to_n_Mg25) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 11.0_rt, 21.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Na21_to_Al25);
    rate_eval.screened_rates(k_He4_Na21_to_Al25) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Na21_to_Al25);
        rate_eval.dscreened_rates_dT(k_He4_Na21_to_Al25) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Na21_to_p_Mg24);
    rate_eval.screened_rates(k_He4_Na21_to_p_Mg24) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Na21_to_p_Mg24);
        rate_eval.dscreened_rates_dT(k_He4_Na21_to_p_Mg24) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 11.0_rt, 22.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Na22_to_Mg23);
    rate_eval.screened_rates(k_p_Na22_to_Mg23) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Na22_to_Mg23);
        rate_eval.dscreened_rates_dT(k_p_Na22_to_Mg23) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Na22_to_He4_Ne19);
    rate_eval.screened_rates(k_p_Na22_to_He4_Ne19) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Na22_to_He4_Ne19);
        rate_eval.dscreened_rates_dT(k_p_Na22_to_He4_Ne19) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 11.0_rt, 22.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Na22_to_Al26);
    rate_eval.screened_rates(k_He4_Na22_to_Al26) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Na22_to_Al26);
        rate_eval.dscreened_rates_dT(k_He4_Na22_to_Al26) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Na22_to_n_Al25);
    rate_eval.screened_rates(k_He4_Na22_to_n_Al25) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Na22_to_n_Al25);
        rate_eval.dscreened_rates_dT(k_He4_Na22_to_n_Al25) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Na22_to_p_Mg25);
    rate_eval.screened_rates(k_He4_Na22_to_p_Mg25) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Na22_to_p_Mg25);
        rate_eval.dscreened_rates_dT(k_He4_Na22_to_p_Mg25) = ratraw * dscor_dt + dratraw_dT * scor;
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

    ratraw = rate_eval.screened_rates(k_p_Na23_to_n_Mg23);
    rate_eval.screened_rates(k_p_Na23_to_n_Mg23) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Na23_to_n_Mg23);
        rate_eval.dscreened_rates_dT(k_p_Na23_to_n_Mg23) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Na23_to_He4_Ne20);
    rate_eval.screened_rates(k_p_Na23_to_He4_Ne20) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Na23_to_He4_Ne20);
        rate_eval.dscreened_rates_dT(k_p_Na23_to_He4_Ne20) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Na23_to_C12_C12);
    rate_eval.screened_rates(k_p_Na23_to_C12_C12) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Na23_to_C12_C12);
        rate_eval.dscreened_rates_dT(k_p_Na23_to_C12_C12) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 11.0_rt, 23.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Na23_to_Al27);
    rate_eval.screened_rates(k_He4_Na23_to_Al27) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Na23_to_Al27);
        rate_eval.dscreened_rates_dT(k_He4_Na23_to_Al27) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Na23_to_n_Al26);
    rate_eval.screened_rates(k_He4_Na23_to_n_Al26) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Na23_to_n_Al26);
        rate_eval.dscreened_rates_dT(k_He4_Na23_to_n_Al26) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Na23_to_p_Mg26);
    rate_eval.screened_rates(k_He4_Na23_to_p_Mg26) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Na23_to_p_Mg26);
        rate_eval.dscreened_rates_dT(k_He4_Na23_to_p_Mg26) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 12.0_rt, 24.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Mg24_to_Al25);
    rate_eval.screened_rates(k_p_Mg24_to_Al25) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mg24_to_Al25);
        rate_eval.dscreened_rates_dT(k_p_Mg24_to_Al25) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Mg24_to_He4_Na21);
    rate_eval.screened_rates(k_p_Mg24_to_He4_Na21) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mg24_to_He4_Na21);
        rate_eval.dscreened_rates_dT(k_p_Mg24_to_He4_Na21) = ratraw * dscor_dt + dratraw_dT * scor;
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

    ratraw = rate_eval.screened_rates(k_He4_Mg24_to_p_Al27);
    rate_eval.screened_rates(k_He4_Mg24_to_p_Al27) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mg24_to_p_Al27);
        rate_eval.dscreened_rates_dT(k_He4_Mg24_to_p_Al27) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Mg24_to_C12_O16);
    rate_eval.screened_rates(k_He4_Mg24_to_C12_O16) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mg24_to_C12_O16);
        rate_eval.dscreened_rates_dT(k_He4_Mg24_to_C12_O16) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 12.0_rt, 25.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Mg25_to_Al26);
    rate_eval.screened_rates(k_p_Mg25_to_Al26) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mg25_to_Al26);
        rate_eval.dscreened_rates_dT(k_p_Mg25_to_Al26) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Mg25_to_n_Al25);
    rate_eval.screened_rates(k_p_Mg25_to_n_Al25) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mg25_to_n_Al25);
        rate_eval.dscreened_rates_dT(k_p_Mg25_to_n_Al25) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Mg25_to_He4_Na22);
    rate_eval.screened_rates(k_p_Mg25_to_He4_Na22) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mg25_to_He4_Na22);
        rate_eval.dscreened_rates_dT(k_p_Mg25_to_He4_Na22) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 12.0_rt, 25.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Mg25_to_Si29);
    rate_eval.screened_rates(k_He4_Mg25_to_Si29) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mg25_to_Si29);
        rate_eval.dscreened_rates_dT(k_He4_Mg25_to_Si29) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Mg25_to_n_Si28);
    rate_eval.screened_rates(k_He4_Mg25_to_n_Si28) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mg25_to_n_Si28);
        rate_eval.dscreened_rates_dT(k_He4_Mg25_to_n_Si28) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 12.0_rt, 26.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Mg26_to_Al27);
    rate_eval.screened_rates(k_p_Mg26_to_Al27) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mg26_to_Al27);
        rate_eval.dscreened_rates_dT(k_p_Mg26_to_Al27) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Mg26_to_n_Al26);
    rate_eval.screened_rates(k_p_Mg26_to_n_Al26) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mg26_to_n_Al26);
        rate_eval.dscreened_rates_dT(k_p_Mg26_to_n_Al26) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Mg26_to_He4_Na23);
    rate_eval.screened_rates(k_p_Mg26_to_He4_Na23) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mg26_to_He4_Na23);
        rate_eval.dscreened_rates_dT(k_p_Mg26_to_He4_Na23) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 12.0_rt, 26.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Mg26_to_Si30);
    rate_eval.screened_rates(k_He4_Mg26_to_Si30) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mg26_to_Si30);
        rate_eval.dscreened_rates_dT(k_He4_Mg26_to_Si30) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Mg26_to_n_Si29);
    rate_eval.screened_rates(k_He4_Mg26_to_n_Si29) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mg26_to_n_Si29);
        rate_eval.dscreened_rates_dT(k_He4_Mg26_to_n_Si29) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 13.0_rt, 25.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Al25_to_P29);
    rate_eval.screened_rates(k_He4_Al25_to_P29) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Al25_to_P29);
        rate_eval.dscreened_rates_dT(k_He4_Al25_to_P29) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Al25_to_p_Si28);
    rate_eval.screened_rates(k_He4_Al25_to_p_Si28) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Al25_to_p_Si28);
        rate_eval.dscreened_rates_dT(k_He4_Al25_to_p_Si28) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 13.0_rt, 26.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Al26_to_P30);
    rate_eval.screened_rates(k_He4_Al26_to_P30) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Al26_to_P30);
        rate_eval.dscreened_rates_dT(k_He4_Al26_to_P30) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Al26_to_n_P29);
    rate_eval.screened_rates(k_He4_Al26_to_n_P29) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Al26_to_n_P29);
        rate_eval.dscreened_rates_dT(k_He4_Al26_to_n_P29) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Al26_to_p_Si29);
    rate_eval.screened_rates(k_He4_Al26_to_p_Si29) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Al26_to_p_Si29);
        rate_eval.dscreened_rates_dT(k_He4_Al26_to_p_Si29) = ratraw * dscor_dt + dratraw_dT * scor;
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

    ratraw = rate_eval.screened_rates(k_p_Al27_to_C12_O16);
    rate_eval.screened_rates(k_p_Al27_to_C12_O16) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Al27_to_C12_O16);
        rate_eval.dscreened_rates_dT(k_p_Al27_to_C12_O16) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 13.0_rt, 27.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Al27_to_P31);
    rate_eval.screened_rates(k_He4_Al27_to_P31) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Al27_to_P31);
        rate_eval.dscreened_rates_dT(k_He4_Al27_to_P31) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Al27_to_n_P30);
    rate_eval.screened_rates(k_He4_Al27_to_n_P30) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Al27_to_n_P30);
        rate_eval.dscreened_rates_dT(k_He4_Al27_to_n_P30) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Al27_to_p_Si30);
    rate_eval.screened_rates(k_He4_Al27_to_p_Si30) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Al27_to_p_Si30);
        rate_eval.dscreened_rates_dT(k_He4_Al27_to_p_Si30) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 14.0_rt, 28.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Si28_to_P29);
    rate_eval.screened_rates(k_p_Si28_to_P29) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Si28_to_P29);
        rate_eval.dscreened_rates_dT(k_p_Si28_to_P29) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Si28_to_He4_Al25);
    rate_eval.screened_rates(k_p_Si28_to_He4_Al25) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Si28_to_He4_Al25);
        rate_eval.dscreened_rates_dT(k_p_Si28_to_He4_Al25) = ratraw * dscor_dt + dratraw_dT * scor;
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

    ratraw = rate_eval.screened_rates(k_He4_Si28_to_p_P31);
    rate_eval.screened_rates(k_He4_Si28_to_p_P31) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Si28_to_p_P31);
        rate_eval.dscreened_rates_dT(k_He4_Si28_to_p_P31) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Si28_to_C12_Ne20);
    rate_eval.screened_rates(k_He4_Si28_to_C12_Ne20) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Si28_to_C12_Ne20);
        rate_eval.dscreened_rates_dT(k_He4_Si28_to_C12_Ne20) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Si28_to_O16_O16);
    rate_eval.screened_rates(k_He4_Si28_to_O16_O16) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Si28_to_O16_O16);
        rate_eval.dscreened_rates_dT(k_He4_Si28_to_O16_O16) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 14.0_rt, 29.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Si29_to_P30);
    rate_eval.screened_rates(k_p_Si29_to_P30) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Si29_to_P30);
        rate_eval.dscreened_rates_dT(k_p_Si29_to_P30) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Si29_to_n_P29);
    rate_eval.screened_rates(k_p_Si29_to_n_P29) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Si29_to_n_P29);
        rate_eval.dscreened_rates_dT(k_p_Si29_to_n_P29) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Si29_to_He4_Al26);
    rate_eval.screened_rates(k_p_Si29_to_He4_Al26) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Si29_to_He4_Al26);
        rate_eval.dscreened_rates_dT(k_p_Si29_to_He4_Al26) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 14.0_rt, 29.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Si29_to_S33);
    rate_eval.screened_rates(k_He4_Si29_to_S33) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Si29_to_S33);
        rate_eval.dscreened_rates_dT(k_He4_Si29_to_S33) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Si29_to_n_S32);
    rate_eval.screened_rates(k_He4_Si29_to_n_S32) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Si29_to_n_S32);
        rate_eval.dscreened_rates_dT(k_He4_Si29_to_n_S32) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Si29_to_p_P32);
    rate_eval.screened_rates(k_He4_Si29_to_p_P32) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Si29_to_p_P32);
        rate_eval.dscreened_rates_dT(k_He4_Si29_to_p_P32) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 14.0_rt, 30.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Si30_to_P31);
    rate_eval.screened_rates(k_p_Si30_to_P31) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Si30_to_P31);
        rate_eval.dscreened_rates_dT(k_p_Si30_to_P31) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Si30_to_n_P30);
    rate_eval.screened_rates(k_p_Si30_to_n_P30) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Si30_to_n_P30);
        rate_eval.dscreened_rates_dT(k_p_Si30_to_n_P30) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Si30_to_He4_Al27);
    rate_eval.screened_rates(k_p_Si30_to_He4_Al27) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Si30_to_He4_Al27);
        rate_eval.dscreened_rates_dT(k_p_Si30_to_He4_Al27) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 14.0_rt, 30.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Si30_to_S34);
    rate_eval.screened_rates(k_He4_Si30_to_S34) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Si30_to_S34);
        rate_eval.dscreened_rates_dT(k_He4_Si30_to_S34) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Si30_to_n_S33);
    rate_eval.screened_rates(k_He4_Si30_to_n_S33) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Si30_to_n_S33);
        rate_eval.dscreened_rates_dT(k_He4_Si30_to_n_S33) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Si30_to_p_P33);
    rate_eval.screened_rates(k_He4_Si30_to_p_P33) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Si30_to_p_P33);
        rate_eval.dscreened_rates_dT(k_He4_Si30_to_p_P33) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 14.0_rt, 31.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Si31_to_P32);
    rate_eval.screened_rates(k_p_Si31_to_P32) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Si31_to_P32);
        rate_eval.dscreened_rates_dT(k_p_Si31_to_P32) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Si31_to_n_P31);
    rate_eval.screened_rates(k_p_Si31_to_n_P31) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Si31_to_n_P31);
        rate_eval.dscreened_rates_dT(k_p_Si31_to_n_P31) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 14.0_rt, 31.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Si31_to_S35);
    rate_eval.screened_rates(k_He4_Si31_to_S35) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Si31_to_S35);
        rate_eval.dscreened_rates_dT(k_He4_Si31_to_S35) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Si31_to_n_S34);
    rate_eval.screened_rates(k_He4_Si31_to_n_S34) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Si31_to_n_S34);
        rate_eval.dscreened_rates_dT(k_He4_Si31_to_n_S34) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 14.0_rt, 32.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Si32_to_P33);
    rate_eval.screened_rates(k_p_Si32_to_P33) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Si32_to_P33);
        rate_eval.dscreened_rates_dT(k_p_Si32_to_P33) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Si32_to_n_P32);
    rate_eval.screened_rates(k_p_Si32_to_n_P32) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Si32_to_n_P32);
        rate_eval.dscreened_rates_dT(k_p_Si32_to_n_P32) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 14.0_rt, 32.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Si32_to_S36);
    rate_eval.screened_rates(k_He4_Si32_to_S36) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Si32_to_S36);
        rate_eval.dscreened_rates_dT(k_He4_Si32_to_S36) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Si32_to_n_S35);
    rate_eval.screened_rates(k_He4_Si32_to_n_S35) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Si32_to_n_S35);
        rate_eval.dscreened_rates_dT(k_He4_Si32_to_n_S35) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 15.0_rt, 29.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_P29_to_Cl33);
    rate_eval.screened_rates(k_He4_P29_to_Cl33) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_P29_to_Cl33);
        rate_eval.dscreened_rates_dT(k_He4_P29_to_Cl33) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_P29_to_p_S32);
    rate_eval.screened_rates(k_He4_P29_to_p_S32) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_P29_to_p_S32);
        rate_eval.dscreened_rates_dT(k_He4_P29_to_p_S32) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 15.0_rt, 30.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_P30_to_Cl34);
    rate_eval.screened_rates(k_He4_P30_to_Cl34) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_P30_to_Cl34);
        rate_eval.dscreened_rates_dT(k_He4_P30_to_Cl34) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_P30_to_n_Cl33);
    rate_eval.screened_rates(k_He4_P30_to_n_Cl33) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_P30_to_n_Cl33);
        rate_eval.dscreened_rates_dT(k_He4_P30_to_n_Cl33) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_P30_to_p_S33);
    rate_eval.screened_rates(k_He4_P30_to_p_S33) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_P30_to_p_S33);
        rate_eval.dscreened_rates_dT(k_He4_P30_to_p_S33) = ratraw * dscor_dt + dratraw_dT * scor;
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

    ratraw = rate_eval.screened_rates(k_p_P31_to_C12_Ne20);
    rate_eval.screened_rates(k_p_P31_to_C12_Ne20) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_P31_to_C12_Ne20);
        rate_eval.dscreened_rates_dT(k_p_P31_to_C12_Ne20) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_P31_to_O16_O16);
    rate_eval.screened_rates(k_p_P31_to_O16_O16) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_P31_to_O16_O16);
        rate_eval.dscreened_rates_dT(k_p_P31_to_O16_O16) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 15.0_rt, 31.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_P31_to_Cl35);
    rate_eval.screened_rates(k_He4_P31_to_Cl35) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_P31_to_Cl35);
        rate_eval.dscreened_rates_dT(k_He4_P31_to_Cl35) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_P31_to_n_Cl34);
    rate_eval.screened_rates(k_He4_P31_to_n_Cl34) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_P31_to_n_Cl34);
        rate_eval.dscreened_rates_dT(k_He4_P31_to_n_Cl34) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_P31_to_p_S34);
    rate_eval.screened_rates(k_He4_P31_to_p_S34) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_P31_to_p_S34);
        rate_eval.dscreened_rates_dT(k_He4_P31_to_p_S34) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 15.0_rt, 32.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_P32_to_S33);
    rate_eval.screened_rates(k_p_P32_to_S33) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_P32_to_S33);
        rate_eval.dscreened_rates_dT(k_p_P32_to_S33) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_P32_to_n_S32);
    rate_eval.screened_rates(k_p_P32_to_n_S32) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_P32_to_n_S32);
        rate_eval.dscreened_rates_dT(k_p_P32_to_n_S32) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_P32_to_He4_Si29);
    rate_eval.screened_rates(k_p_P32_to_He4_Si29) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_P32_to_He4_Si29);
        rate_eval.dscreened_rates_dT(k_p_P32_to_He4_Si29) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 15.0_rt, 32.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_P32_to_Cl36);
    rate_eval.screened_rates(k_He4_P32_to_Cl36) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_P32_to_Cl36);
        rate_eval.dscreened_rates_dT(k_He4_P32_to_Cl36) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_P32_to_n_Cl35);
    rate_eval.screened_rates(k_He4_P32_to_n_Cl35) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_P32_to_n_Cl35);
        rate_eval.dscreened_rates_dT(k_He4_P32_to_n_Cl35) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_P32_to_p_S35);
    rate_eval.screened_rates(k_He4_P32_to_p_S35) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_P32_to_p_S35);
        rate_eval.dscreened_rates_dT(k_He4_P32_to_p_S35) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 15.0_rt, 33.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_P33_to_S34);
    rate_eval.screened_rates(k_p_P33_to_S34) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_P33_to_S34);
        rate_eval.dscreened_rates_dT(k_p_P33_to_S34) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_P33_to_n_S33);
    rate_eval.screened_rates(k_p_P33_to_n_S33) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_P33_to_n_S33);
        rate_eval.dscreened_rates_dT(k_p_P33_to_n_S33) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_P33_to_He4_Si30);
    rate_eval.screened_rates(k_p_P33_to_He4_Si30) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_P33_to_He4_Si30);
        rate_eval.dscreened_rates_dT(k_p_P33_to_He4_Si30) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 15.0_rt, 33.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_P33_to_Cl37);
    rate_eval.screened_rates(k_He4_P33_to_Cl37) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_P33_to_Cl37);
        rate_eval.dscreened_rates_dT(k_He4_P33_to_Cl37) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_P33_to_n_Cl36);
    rate_eval.screened_rates(k_He4_P33_to_n_Cl36) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_P33_to_n_Cl36);
        rate_eval.dscreened_rates_dT(k_He4_P33_to_n_Cl36) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_P33_to_p_S36);
    rate_eval.screened_rates(k_He4_P33_to_p_S36) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_P33_to_p_S36);
        rate_eval.dscreened_rates_dT(k_He4_P33_to_p_S36) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 16.0_rt, 32.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_S32_to_Cl33);
    rate_eval.screened_rates(k_p_S32_to_Cl33) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_S32_to_Cl33);
        rate_eval.dscreened_rates_dT(k_p_S32_to_Cl33) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_S32_to_He4_P29);
    rate_eval.screened_rates(k_p_S32_to_He4_P29) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_S32_to_He4_P29);
        rate_eval.dscreened_rates_dT(k_p_S32_to_He4_P29) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 16.0_rt, 32.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_S32_to_Ar36);
    rate_eval.screened_rates(k_He4_S32_to_Ar36) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_S32_to_Ar36);
        rate_eval.dscreened_rates_dT(k_He4_S32_to_Ar36) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_S32_to_p_Cl35);
    rate_eval.screened_rates(k_He4_S32_to_p_Cl35) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_S32_to_p_Cl35);
        rate_eval.dscreened_rates_dT(k_He4_S32_to_p_Cl35) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 16.0_rt, 33.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_S33_to_Cl34);
    rate_eval.screened_rates(k_p_S33_to_Cl34) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_S33_to_Cl34);
        rate_eval.dscreened_rates_dT(k_p_S33_to_Cl34) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_S33_to_n_Cl33);
    rate_eval.screened_rates(k_p_S33_to_n_Cl33) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_S33_to_n_Cl33);
        rate_eval.dscreened_rates_dT(k_p_S33_to_n_Cl33) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_S33_to_He4_P30);
    rate_eval.screened_rates(k_p_S33_to_He4_P30) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_S33_to_He4_P30);
        rate_eval.dscreened_rates_dT(k_p_S33_to_He4_P30) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 16.0_rt, 33.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_S33_to_Ar37);
    rate_eval.screened_rates(k_He4_S33_to_Ar37) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_S33_to_Ar37);
        rate_eval.dscreened_rates_dT(k_He4_S33_to_Ar37) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_S33_to_n_Ar36);
    rate_eval.screened_rates(k_He4_S33_to_n_Ar36) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_S33_to_n_Ar36);
        rate_eval.dscreened_rates_dT(k_He4_S33_to_n_Ar36) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_S33_to_p_Cl36);
    rate_eval.screened_rates(k_He4_S33_to_p_Cl36) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_S33_to_p_Cl36);
        rate_eval.dscreened_rates_dT(k_He4_S33_to_p_Cl36) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 16.0_rt, 34.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_S34_to_Cl35);
    rate_eval.screened_rates(k_p_S34_to_Cl35) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_S34_to_Cl35);
        rate_eval.dscreened_rates_dT(k_p_S34_to_Cl35) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_S34_to_n_Cl34);
    rate_eval.screened_rates(k_p_S34_to_n_Cl34) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_S34_to_n_Cl34);
        rate_eval.dscreened_rates_dT(k_p_S34_to_n_Cl34) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_S34_to_He4_P31);
    rate_eval.screened_rates(k_p_S34_to_He4_P31) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_S34_to_He4_P31);
        rate_eval.dscreened_rates_dT(k_p_S34_to_He4_P31) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 16.0_rt, 34.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_S34_to_Ar38);
    rate_eval.screened_rates(k_He4_S34_to_Ar38) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_S34_to_Ar38);
        rate_eval.dscreened_rates_dT(k_He4_S34_to_Ar38) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_S34_to_n_Ar37);
    rate_eval.screened_rates(k_He4_S34_to_n_Ar37) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_S34_to_n_Ar37);
        rate_eval.dscreened_rates_dT(k_He4_S34_to_n_Ar37) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_S34_to_p_Cl37);
    rate_eval.screened_rates(k_He4_S34_to_p_Cl37) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_S34_to_p_Cl37);
        rate_eval.dscreened_rates_dT(k_He4_S34_to_p_Cl37) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 16.0_rt, 35.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_S35_to_Cl36);
    rate_eval.screened_rates(k_p_S35_to_Cl36) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_S35_to_Cl36);
        rate_eval.dscreened_rates_dT(k_p_S35_to_Cl36) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_S35_to_n_Cl35);
    rate_eval.screened_rates(k_p_S35_to_n_Cl35) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_S35_to_n_Cl35);
        rate_eval.dscreened_rates_dT(k_p_S35_to_n_Cl35) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_S35_to_He4_P32);
    rate_eval.screened_rates(k_p_S35_to_He4_P32) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_S35_to_He4_P32);
        rate_eval.dscreened_rates_dT(k_p_S35_to_He4_P32) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 16.0_rt, 35.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_S35_to_Ar39);
    rate_eval.screened_rates(k_He4_S35_to_Ar39) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_S35_to_Ar39);
        rate_eval.dscreened_rates_dT(k_He4_S35_to_Ar39) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_S35_to_n_Ar38);
    rate_eval.screened_rates(k_He4_S35_to_n_Ar38) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_S35_to_n_Ar38);
        rate_eval.dscreened_rates_dT(k_He4_S35_to_n_Ar38) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 16.0_rt, 36.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_S36_to_Cl37);
    rate_eval.screened_rates(k_p_S36_to_Cl37) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_S36_to_Cl37);
        rate_eval.dscreened_rates_dT(k_p_S36_to_Cl37) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_S36_to_n_Cl36);
    rate_eval.screened_rates(k_p_S36_to_n_Cl36) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_S36_to_n_Cl36);
        rate_eval.dscreened_rates_dT(k_p_S36_to_n_Cl36) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_S36_to_He4_P33);
    rate_eval.screened_rates(k_p_S36_to_He4_P33) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_S36_to_He4_P33);
        rate_eval.dscreened_rates_dT(k_p_S36_to_He4_P33) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 16.0_rt, 36.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_S36_to_Ar40);
    rate_eval.screened_rates(k_He4_S36_to_Ar40) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_S36_to_Ar40);
        rate_eval.dscreened_rates_dT(k_He4_S36_to_Ar40) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_S36_to_n_Ar39);
    rate_eval.screened_rates(k_He4_S36_to_n_Ar39) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_S36_to_n_Ar39);
        rate_eval.dscreened_rates_dT(k_He4_S36_to_n_Ar39) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 17.0_rt, 33.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Cl33_to_K37);
    rate_eval.screened_rates(k_He4_Cl33_to_K37) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cl33_to_K37);
        rate_eval.dscreened_rates_dT(k_He4_Cl33_to_K37) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cl33_to_p_Ar36);
    rate_eval.screened_rates(k_He4_Cl33_to_p_Ar36) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cl33_to_p_Ar36);
        rate_eval.dscreened_rates_dT(k_He4_Cl33_to_p_Ar36) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 17.0_rt, 34.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Cl34_to_K38);
    rate_eval.screened_rates(k_He4_Cl34_to_K38) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cl34_to_K38);
        rate_eval.dscreened_rates_dT(k_He4_Cl34_to_K38) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cl34_to_n_K37);
    rate_eval.screened_rates(k_He4_Cl34_to_n_K37) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cl34_to_n_K37);
        rate_eval.dscreened_rates_dT(k_He4_Cl34_to_n_K37) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cl34_to_p_Ar37);
    rate_eval.screened_rates(k_He4_Cl34_to_p_Ar37) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cl34_to_p_Ar37);
        rate_eval.dscreened_rates_dT(k_He4_Cl34_to_p_Ar37) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 17.0_rt, 35.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Cl35_to_Ar36);
    rate_eval.screened_rates(k_p_Cl35_to_Ar36) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cl35_to_Ar36);
        rate_eval.dscreened_rates_dT(k_p_Cl35_to_Ar36) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cl35_to_He4_S32);
    rate_eval.screened_rates(k_p_Cl35_to_He4_S32) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cl35_to_He4_S32);
        rate_eval.dscreened_rates_dT(k_p_Cl35_to_He4_S32) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 17.0_rt, 35.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Cl35_to_K39);
    rate_eval.screened_rates(k_He4_Cl35_to_K39) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cl35_to_K39);
        rate_eval.dscreened_rates_dT(k_He4_Cl35_to_K39) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cl35_to_n_K38);
    rate_eval.screened_rates(k_He4_Cl35_to_n_K38) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cl35_to_n_K38);
        rate_eval.dscreened_rates_dT(k_He4_Cl35_to_n_K38) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cl35_to_p_Ar38);
    rate_eval.screened_rates(k_He4_Cl35_to_p_Ar38) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cl35_to_p_Ar38);
        rate_eval.dscreened_rates_dT(k_He4_Cl35_to_p_Ar38) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 17.0_rt, 36.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Cl36_to_Ar37);
    rate_eval.screened_rates(k_p_Cl36_to_Ar37) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cl36_to_Ar37);
        rate_eval.dscreened_rates_dT(k_p_Cl36_to_Ar37) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cl36_to_n_Ar36);
    rate_eval.screened_rates(k_p_Cl36_to_n_Ar36) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cl36_to_n_Ar36);
        rate_eval.dscreened_rates_dT(k_p_Cl36_to_n_Ar36) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cl36_to_He4_S33);
    rate_eval.screened_rates(k_p_Cl36_to_He4_S33) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cl36_to_He4_S33);
        rate_eval.dscreened_rates_dT(k_p_Cl36_to_He4_S33) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 17.0_rt, 36.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Cl36_to_K40);
    rate_eval.screened_rates(k_He4_Cl36_to_K40) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cl36_to_K40);
        rate_eval.dscreened_rates_dT(k_He4_Cl36_to_K40) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cl36_to_n_K39);
    rate_eval.screened_rates(k_He4_Cl36_to_n_K39) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cl36_to_n_K39);
        rate_eval.dscreened_rates_dT(k_He4_Cl36_to_n_K39) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cl36_to_p_Ar39);
    rate_eval.screened_rates(k_He4_Cl36_to_p_Ar39) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cl36_to_p_Ar39);
        rate_eval.dscreened_rates_dT(k_He4_Cl36_to_p_Ar39) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 17.0_rt, 37.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Cl37_to_Ar38);
    rate_eval.screened_rates(k_p_Cl37_to_Ar38) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cl37_to_Ar38);
        rate_eval.dscreened_rates_dT(k_p_Cl37_to_Ar38) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cl37_to_n_Ar37);
    rate_eval.screened_rates(k_p_Cl37_to_n_Ar37) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cl37_to_n_Ar37);
        rate_eval.dscreened_rates_dT(k_p_Cl37_to_n_Ar37) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cl37_to_He4_S34);
    rate_eval.screened_rates(k_p_Cl37_to_He4_S34) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cl37_to_He4_S34);
        rate_eval.dscreened_rates_dT(k_p_Cl37_to_He4_S34) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 17.0_rt, 37.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Cl37_to_K41);
    rate_eval.screened_rates(k_He4_Cl37_to_K41) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cl37_to_K41);
        rate_eval.dscreened_rates_dT(k_He4_Cl37_to_K41) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cl37_to_n_K40);
    rate_eval.screened_rates(k_He4_Cl37_to_n_K40) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cl37_to_n_K40);
        rate_eval.dscreened_rates_dT(k_He4_Cl37_to_n_K40) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cl37_to_p_Ar40);
    rate_eval.screened_rates(k_He4_Cl37_to_p_Ar40) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cl37_to_p_Ar40);
        rate_eval.dscreened_rates_dT(k_He4_Cl37_to_p_Ar40) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 18.0_rt, 36.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ar36_to_K37);
    rate_eval.screened_rates(k_p_Ar36_to_K37) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ar36_to_K37);
        rate_eval.dscreened_rates_dT(k_p_Ar36_to_K37) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ar36_to_He4_Cl33);
    rate_eval.screened_rates(k_p_Ar36_to_He4_Cl33) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ar36_to_He4_Cl33);
        rate_eval.dscreened_rates_dT(k_p_Ar36_to_He4_Cl33) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 18.0_rt, 36.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ar36_to_Ca40);
    rate_eval.screened_rates(k_He4_Ar36_to_Ca40) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ar36_to_Ca40);
        rate_eval.dscreened_rates_dT(k_He4_Ar36_to_Ca40) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ar36_to_p_K39);
    rate_eval.screened_rates(k_He4_Ar36_to_p_K39) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ar36_to_p_K39);
        rate_eval.dscreened_rates_dT(k_He4_Ar36_to_p_K39) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 18.0_rt, 37.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ar37_to_K38);
    rate_eval.screened_rates(k_p_Ar37_to_K38) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ar37_to_K38);
        rate_eval.dscreened_rates_dT(k_p_Ar37_to_K38) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ar37_to_n_K37);
    rate_eval.screened_rates(k_p_Ar37_to_n_K37) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ar37_to_n_K37);
        rate_eval.dscreened_rates_dT(k_p_Ar37_to_n_K37) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ar37_to_He4_Cl34);
    rate_eval.screened_rates(k_p_Ar37_to_He4_Cl34) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ar37_to_He4_Cl34);
        rate_eval.dscreened_rates_dT(k_p_Ar37_to_He4_Cl34) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 18.0_rt, 37.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ar37_to_Ca41);
    rate_eval.screened_rates(k_He4_Ar37_to_Ca41) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ar37_to_Ca41);
        rate_eval.dscreened_rates_dT(k_He4_Ar37_to_Ca41) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ar37_to_n_Ca40);
    rate_eval.screened_rates(k_He4_Ar37_to_n_Ca40) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ar37_to_n_Ca40);
        rate_eval.dscreened_rates_dT(k_He4_Ar37_to_n_Ca40) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ar37_to_p_K40);
    rate_eval.screened_rates(k_He4_Ar37_to_p_K40) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ar37_to_p_K40);
        rate_eval.dscreened_rates_dT(k_He4_Ar37_to_p_K40) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 18.0_rt, 38.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ar38_to_K39);
    rate_eval.screened_rates(k_p_Ar38_to_K39) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ar38_to_K39);
        rate_eval.dscreened_rates_dT(k_p_Ar38_to_K39) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ar38_to_n_K38);
    rate_eval.screened_rates(k_p_Ar38_to_n_K38) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ar38_to_n_K38);
        rate_eval.dscreened_rates_dT(k_p_Ar38_to_n_K38) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ar38_to_He4_Cl35);
    rate_eval.screened_rates(k_p_Ar38_to_He4_Cl35) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ar38_to_He4_Cl35);
        rate_eval.dscreened_rates_dT(k_p_Ar38_to_He4_Cl35) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 18.0_rt, 38.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ar38_to_Ca42);
    rate_eval.screened_rates(k_He4_Ar38_to_Ca42) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ar38_to_Ca42);
        rate_eval.dscreened_rates_dT(k_He4_Ar38_to_Ca42) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ar38_to_n_Ca41);
    rate_eval.screened_rates(k_He4_Ar38_to_n_Ca41) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ar38_to_n_Ca41);
        rate_eval.dscreened_rates_dT(k_He4_Ar38_to_n_Ca41) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ar38_to_p_K41);
    rate_eval.screened_rates(k_He4_Ar38_to_p_K41) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ar38_to_p_K41);
        rate_eval.dscreened_rates_dT(k_He4_Ar38_to_p_K41) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 18.0_rt, 39.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ar39_to_K40);
    rate_eval.screened_rates(k_p_Ar39_to_K40) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ar39_to_K40);
        rate_eval.dscreened_rates_dT(k_p_Ar39_to_K40) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ar39_to_n_K39);
    rate_eval.screened_rates(k_p_Ar39_to_n_K39) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ar39_to_n_K39);
        rate_eval.dscreened_rates_dT(k_p_Ar39_to_n_K39) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ar39_to_He4_Cl36);
    rate_eval.screened_rates(k_p_Ar39_to_He4_Cl36) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ar39_to_He4_Cl36);
        rate_eval.dscreened_rates_dT(k_p_Ar39_to_He4_Cl36) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 18.0_rt, 39.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ar39_to_Ca43);
    rate_eval.screened_rates(k_He4_Ar39_to_Ca43) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ar39_to_Ca43);
        rate_eval.dscreened_rates_dT(k_He4_Ar39_to_Ca43) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ar39_to_n_Ca42);
    rate_eval.screened_rates(k_He4_Ar39_to_n_Ca42) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ar39_to_n_Ca42);
        rate_eval.dscreened_rates_dT(k_He4_Ar39_to_n_Ca42) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 18.0_rt, 40.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ar40_to_K41);
    rate_eval.screened_rates(k_p_Ar40_to_K41) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ar40_to_K41);
        rate_eval.dscreened_rates_dT(k_p_Ar40_to_K41) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ar40_to_n_K40);
    rate_eval.screened_rates(k_p_Ar40_to_n_K40) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ar40_to_n_K40);
        rate_eval.dscreened_rates_dT(k_p_Ar40_to_n_K40) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ar40_to_He4_Cl37);
    rate_eval.screened_rates(k_p_Ar40_to_He4_Cl37) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ar40_to_He4_Cl37);
        rate_eval.dscreened_rates_dT(k_p_Ar40_to_He4_Cl37) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 18.0_rt, 40.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ar40_to_Ca44);
    rate_eval.screened_rates(k_He4_Ar40_to_Ca44) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ar40_to_Ca44);
        rate_eval.dscreened_rates_dT(k_He4_Ar40_to_Ca44) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ar40_to_n_Ca43);
    rate_eval.screened_rates(k_He4_Ar40_to_n_Ca43) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ar40_to_n_Ca43);
        rate_eval.dscreened_rates_dT(k_He4_Ar40_to_n_Ca43) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 19.0_rt, 39.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_K39_to_Ca40);
    rate_eval.screened_rates(k_p_K39_to_Ca40) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_K39_to_Ca40);
        rate_eval.dscreened_rates_dT(k_p_K39_to_Ca40) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_K39_to_He4_Ar36);
    rate_eval.screened_rates(k_p_K39_to_He4_Ar36) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_K39_to_He4_Ar36);
        rate_eval.dscreened_rates_dT(k_p_K39_to_He4_Ar36) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 19.0_rt, 39.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_K39_to_Sc43);
    rate_eval.screened_rates(k_He4_K39_to_Sc43) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_K39_to_Sc43);
        rate_eval.dscreened_rates_dT(k_He4_K39_to_Sc43) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_K39_to_p_Ca42);
    rate_eval.screened_rates(k_He4_K39_to_p_Ca42) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_K39_to_p_Ca42);
        rate_eval.dscreened_rates_dT(k_He4_K39_to_p_Ca42) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 19.0_rt, 40.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_K40_to_Ca41);
    rate_eval.screened_rates(k_p_K40_to_Ca41) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_K40_to_Ca41);
        rate_eval.dscreened_rates_dT(k_p_K40_to_Ca41) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_K40_to_n_Ca40);
    rate_eval.screened_rates(k_p_K40_to_n_Ca40) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_K40_to_n_Ca40);
        rate_eval.dscreened_rates_dT(k_p_K40_to_n_Ca40) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_K40_to_He4_Ar37);
    rate_eval.screened_rates(k_p_K40_to_He4_Ar37) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_K40_to_He4_Ar37);
        rate_eval.dscreened_rates_dT(k_p_K40_to_He4_Ar37) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 19.0_rt, 40.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_K40_to_Sc44);
    rate_eval.screened_rates(k_He4_K40_to_Sc44) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_K40_to_Sc44);
        rate_eval.dscreened_rates_dT(k_He4_K40_to_Sc44) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_K40_to_n_Sc43);
    rate_eval.screened_rates(k_He4_K40_to_n_Sc43) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_K40_to_n_Sc43);
        rate_eval.dscreened_rates_dT(k_He4_K40_to_n_Sc43) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_K40_to_p_Ca43);
    rate_eval.screened_rates(k_He4_K40_to_p_Ca43) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_K40_to_p_Ca43);
        rate_eval.dscreened_rates_dT(k_He4_K40_to_p_Ca43) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 19.0_rt, 41.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_K41_to_Ca42);
    rate_eval.screened_rates(k_p_K41_to_Ca42) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_K41_to_Ca42);
        rate_eval.dscreened_rates_dT(k_p_K41_to_Ca42) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_K41_to_n_Ca41);
    rate_eval.screened_rates(k_p_K41_to_n_Ca41) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_K41_to_n_Ca41);
        rate_eval.dscreened_rates_dT(k_p_K41_to_n_Ca41) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_K41_to_He4_Ar38);
    rate_eval.screened_rates(k_p_K41_to_He4_Ar38) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_K41_to_He4_Ar38);
        rate_eval.dscreened_rates_dT(k_p_K41_to_He4_Ar38) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 19.0_rt, 41.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_K41_to_Sc45);
    rate_eval.screened_rates(k_He4_K41_to_Sc45) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_K41_to_Sc45);
        rate_eval.dscreened_rates_dT(k_He4_K41_to_Sc45) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_K41_to_n_Sc44);
    rate_eval.screened_rates(k_He4_K41_to_n_Sc44) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_K41_to_n_Sc44);
        rate_eval.dscreened_rates_dT(k_He4_K41_to_n_Sc44) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_K41_to_p_Ca44);
    rate_eval.screened_rates(k_He4_K41_to_p_Ca44) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_K41_to_p_Ca44);
        rate_eval.dscreened_rates_dT(k_He4_K41_to_p_Ca44) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 20.0_rt, 40.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ca40_to_Ti44);
    rate_eval.screened_rates(k_He4_Ca40_to_Ti44) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca40_to_Ti44);
        rate_eval.dscreened_rates_dT(k_He4_Ca40_to_Ti44) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ca40_to_p_Sc43);
    rate_eval.screened_rates(k_He4_Ca40_to_p_Sc43) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca40_to_p_Sc43);
        rate_eval.dscreened_rates_dT(k_He4_Ca40_to_p_Sc43) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 20.0_rt, 41.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ca41_to_Ti45);
    rate_eval.screened_rates(k_He4_Ca41_to_Ti45) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca41_to_Ti45);
        rate_eval.dscreened_rates_dT(k_He4_Ca41_to_Ti45) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ca41_to_n_Ti44);
    rate_eval.screened_rates(k_He4_Ca41_to_n_Ti44) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca41_to_n_Ti44);
        rate_eval.dscreened_rates_dT(k_He4_Ca41_to_n_Ti44) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ca41_to_p_Sc44);
    rate_eval.screened_rates(k_He4_Ca41_to_p_Sc44) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca41_to_p_Sc44);
        rate_eval.dscreened_rates_dT(k_He4_Ca41_to_p_Sc44) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 20.0_rt, 42.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ca42_to_Sc43);
    rate_eval.screened_rates(k_p_Ca42_to_Sc43) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ca42_to_Sc43);
        rate_eval.dscreened_rates_dT(k_p_Ca42_to_Sc43) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ca42_to_He4_K39);
    rate_eval.screened_rates(k_p_Ca42_to_He4_K39) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ca42_to_He4_K39);
        rate_eval.dscreened_rates_dT(k_p_Ca42_to_He4_K39) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 20.0_rt, 42.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ca42_to_Ti46);
    rate_eval.screened_rates(k_He4_Ca42_to_Ti46) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca42_to_Ti46);
        rate_eval.dscreened_rates_dT(k_He4_Ca42_to_Ti46) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ca42_to_n_Ti45);
    rate_eval.screened_rates(k_He4_Ca42_to_n_Ti45) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca42_to_n_Ti45);
        rate_eval.dscreened_rates_dT(k_He4_Ca42_to_n_Ti45) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ca42_to_p_Sc45);
    rate_eval.screened_rates(k_He4_Ca42_to_p_Sc45) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca42_to_p_Sc45);
        rate_eval.dscreened_rates_dT(k_He4_Ca42_to_p_Sc45) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 20.0_rt, 43.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ca43_to_Sc44);
    rate_eval.screened_rates(k_p_Ca43_to_Sc44) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ca43_to_Sc44);
        rate_eval.dscreened_rates_dT(k_p_Ca43_to_Sc44) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ca43_to_n_Sc43);
    rate_eval.screened_rates(k_p_Ca43_to_n_Sc43) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ca43_to_n_Sc43);
        rate_eval.dscreened_rates_dT(k_p_Ca43_to_n_Sc43) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ca43_to_He4_K40);
    rate_eval.screened_rates(k_p_Ca43_to_He4_K40) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ca43_to_He4_K40);
        rate_eval.dscreened_rates_dT(k_p_Ca43_to_He4_K40) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 20.0_rt, 43.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ca43_to_Ti47);
    rate_eval.screened_rates(k_He4_Ca43_to_Ti47) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca43_to_Ti47);
        rate_eval.dscreened_rates_dT(k_He4_Ca43_to_Ti47) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ca43_to_n_Ti46);
    rate_eval.screened_rates(k_He4_Ca43_to_n_Ti46) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca43_to_n_Ti46);
        rate_eval.dscreened_rates_dT(k_He4_Ca43_to_n_Ti46) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ca43_to_p_Sc46);
    rate_eval.screened_rates(k_He4_Ca43_to_p_Sc46) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca43_to_p_Sc46);
        rate_eval.dscreened_rates_dT(k_He4_Ca43_to_p_Sc46) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 20.0_rt, 44.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ca44_to_Sc45);
    rate_eval.screened_rates(k_p_Ca44_to_Sc45) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ca44_to_Sc45);
        rate_eval.dscreened_rates_dT(k_p_Ca44_to_Sc45) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ca44_to_n_Sc44);
    rate_eval.screened_rates(k_p_Ca44_to_n_Sc44) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ca44_to_n_Sc44);
        rate_eval.dscreened_rates_dT(k_p_Ca44_to_n_Sc44) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ca44_to_He4_K41);
    rate_eval.screened_rates(k_p_Ca44_to_He4_K41) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ca44_to_He4_K41);
        rate_eval.dscreened_rates_dT(k_p_Ca44_to_He4_K41) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 20.0_rt, 44.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ca44_to_Ti48);
    rate_eval.screened_rates(k_He4_Ca44_to_Ti48) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca44_to_Ti48);
        rate_eval.dscreened_rates_dT(k_He4_Ca44_to_Ti48) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ca44_to_n_Ti47);
    rate_eval.screened_rates(k_He4_Ca44_to_n_Ti47) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca44_to_n_Ti47);
        rate_eval.dscreened_rates_dT(k_He4_Ca44_to_n_Ti47) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ca44_to_p_Sc47);
    rate_eval.screened_rates(k_He4_Ca44_to_p_Sc47) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca44_to_p_Sc47);
        rate_eval.dscreened_rates_dT(k_He4_Ca44_to_p_Sc47) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 20.0_rt, 45.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ca45_to_Sc46);
    rate_eval.screened_rates(k_p_Ca45_to_Sc46) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ca45_to_Sc46);
        rate_eval.dscreened_rates_dT(k_p_Ca45_to_Sc46) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ca45_to_n_Sc45);
    rate_eval.screened_rates(k_p_Ca45_to_n_Sc45) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ca45_to_n_Sc45);
        rate_eval.dscreened_rates_dT(k_p_Ca45_to_n_Sc45) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 20.0_rt, 45.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ca45_to_Ti49);
    rate_eval.screened_rates(k_He4_Ca45_to_Ti49) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca45_to_Ti49);
        rate_eval.dscreened_rates_dT(k_He4_Ca45_to_Ti49) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ca45_to_n_Ti48);
    rate_eval.screened_rates(k_He4_Ca45_to_n_Ti48) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca45_to_n_Ti48);
        rate_eval.dscreened_rates_dT(k_He4_Ca45_to_n_Ti48) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ca45_to_p_Sc48);
    rate_eval.screened_rates(k_He4_Ca45_to_p_Sc48) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca45_to_p_Sc48);
        rate_eval.dscreened_rates_dT(k_He4_Ca45_to_p_Sc48) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 20.0_rt, 46.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ca46_to_Sc47);
    rate_eval.screened_rates(k_p_Ca46_to_Sc47) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ca46_to_Sc47);
        rate_eval.dscreened_rates_dT(k_p_Ca46_to_Sc47) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ca46_to_n_Sc46);
    rate_eval.screened_rates(k_p_Ca46_to_n_Sc46) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ca46_to_n_Sc46);
        rate_eval.dscreened_rates_dT(k_p_Ca46_to_n_Sc46) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 20.0_rt, 46.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ca46_to_Ti50);
    rate_eval.screened_rates(k_He4_Ca46_to_Ti50) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca46_to_Ti50);
        rate_eval.dscreened_rates_dT(k_He4_Ca46_to_Ti50) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ca46_to_n_Ti49);
    rate_eval.screened_rates(k_He4_Ca46_to_n_Ti49) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca46_to_n_Ti49);
        rate_eval.dscreened_rates_dT(k_He4_Ca46_to_n_Ti49) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ca46_to_p_Sc49);
    rate_eval.screened_rates(k_He4_Ca46_to_p_Sc49) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca46_to_p_Sc49);
        rate_eval.dscreened_rates_dT(k_He4_Ca46_to_p_Sc49) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 20.0_rt, 47.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ca47_to_Sc48);
    rate_eval.screened_rates(k_p_Ca47_to_Sc48) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ca47_to_Sc48);
        rate_eval.dscreened_rates_dT(k_p_Ca47_to_Sc48) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ca47_to_n_Sc47);
    rate_eval.screened_rates(k_p_Ca47_to_n_Sc47) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ca47_to_n_Sc47);
        rate_eval.dscreened_rates_dT(k_p_Ca47_to_n_Sc47) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 20.0_rt, 47.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ca47_to_Ti51);
    rate_eval.screened_rates(k_He4_Ca47_to_Ti51) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca47_to_Ti51);
        rate_eval.dscreened_rates_dT(k_He4_Ca47_to_Ti51) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ca47_to_n_Ti50);
    rate_eval.screened_rates(k_He4_Ca47_to_n_Ti50) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca47_to_n_Ti50);
        rate_eval.dscreened_rates_dT(k_He4_Ca47_to_n_Ti50) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 20.0_rt, 48.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ca48_to_Sc49);
    rate_eval.screened_rates(k_p_Ca48_to_Sc49) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ca48_to_Sc49);
        rate_eval.dscreened_rates_dT(k_p_Ca48_to_Sc49) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ca48_to_n_Sc48);
    rate_eval.screened_rates(k_p_Ca48_to_n_Sc48) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ca48_to_n_Sc48);
        rate_eval.dscreened_rates_dT(k_p_Ca48_to_n_Sc48) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 21.0_rt, 43.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Sc43_to_Ti44);
    rate_eval.screened_rates(k_p_Sc43_to_Ti44) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Sc43_to_Ti44);
        rate_eval.dscreened_rates_dT(k_p_Sc43_to_Ti44) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Sc43_to_He4_Ca40);
    rate_eval.screened_rates(k_p_Sc43_to_He4_Ca40) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Sc43_to_He4_Ca40);
        rate_eval.dscreened_rates_dT(k_p_Sc43_to_He4_Ca40) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 21.0_rt, 43.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Sc43_to_V47);
    rate_eval.screened_rates(k_He4_Sc43_to_V47) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Sc43_to_V47);
        rate_eval.dscreened_rates_dT(k_He4_Sc43_to_V47) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Sc43_to_n_V46);
    rate_eval.screened_rates(k_He4_Sc43_to_n_V46) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Sc43_to_n_V46);
        rate_eval.dscreened_rates_dT(k_He4_Sc43_to_n_V46) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Sc43_to_p_Ti46);
    rate_eval.screened_rates(k_He4_Sc43_to_p_Ti46) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Sc43_to_p_Ti46);
        rate_eval.dscreened_rates_dT(k_He4_Sc43_to_p_Ti46) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 21.0_rt, 44.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Sc44_to_Ti45);
    rate_eval.screened_rates(k_p_Sc44_to_Ti45) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Sc44_to_Ti45);
        rate_eval.dscreened_rates_dT(k_p_Sc44_to_Ti45) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Sc44_to_n_Ti44);
    rate_eval.screened_rates(k_p_Sc44_to_n_Ti44) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Sc44_to_n_Ti44);
        rate_eval.dscreened_rates_dT(k_p_Sc44_to_n_Ti44) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Sc44_to_He4_Ca41);
    rate_eval.screened_rates(k_p_Sc44_to_He4_Ca41) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Sc44_to_He4_Ca41);
        rate_eval.dscreened_rates_dT(k_p_Sc44_to_He4_Ca41) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 21.0_rt, 44.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Sc44_to_V48);
    rate_eval.screened_rates(k_He4_Sc44_to_V48) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Sc44_to_V48);
        rate_eval.dscreened_rates_dT(k_He4_Sc44_to_V48) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Sc44_to_n_V47);
    rate_eval.screened_rates(k_He4_Sc44_to_n_V47) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Sc44_to_n_V47);
        rate_eval.dscreened_rates_dT(k_He4_Sc44_to_n_V47) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Sc44_to_p_Ti47);
    rate_eval.screened_rates(k_He4_Sc44_to_p_Ti47) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Sc44_to_p_Ti47);
        rate_eval.dscreened_rates_dT(k_He4_Sc44_to_p_Ti47) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 21.0_rt, 45.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Sc45_to_Ti46);
    rate_eval.screened_rates(k_p_Sc45_to_Ti46) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Sc45_to_Ti46);
        rate_eval.dscreened_rates_dT(k_p_Sc45_to_Ti46) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Sc45_to_n_Ti45);
    rate_eval.screened_rates(k_p_Sc45_to_n_Ti45) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Sc45_to_n_Ti45);
        rate_eval.dscreened_rates_dT(k_p_Sc45_to_n_Ti45) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Sc45_to_He4_Ca42);
    rate_eval.screened_rates(k_p_Sc45_to_He4_Ca42) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Sc45_to_He4_Ca42);
        rate_eval.dscreened_rates_dT(k_p_Sc45_to_He4_Ca42) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 21.0_rt, 45.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Sc45_to_V49);
    rate_eval.screened_rates(k_He4_Sc45_to_V49) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Sc45_to_V49);
        rate_eval.dscreened_rates_dT(k_He4_Sc45_to_V49) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Sc45_to_n_V48);
    rate_eval.screened_rates(k_He4_Sc45_to_n_V48) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Sc45_to_n_V48);
        rate_eval.dscreened_rates_dT(k_He4_Sc45_to_n_V48) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Sc45_to_p_Ti48);
    rate_eval.screened_rates(k_He4_Sc45_to_p_Ti48) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Sc45_to_p_Ti48);
        rate_eval.dscreened_rates_dT(k_He4_Sc45_to_p_Ti48) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 21.0_rt, 46.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Sc46_to_Ti47);
    rate_eval.screened_rates(k_p_Sc46_to_Ti47) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Sc46_to_Ti47);
        rate_eval.dscreened_rates_dT(k_p_Sc46_to_Ti47) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Sc46_to_n_Ti46);
    rate_eval.screened_rates(k_p_Sc46_to_n_Ti46) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Sc46_to_n_Ti46);
        rate_eval.dscreened_rates_dT(k_p_Sc46_to_n_Ti46) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Sc46_to_He4_Ca43);
    rate_eval.screened_rates(k_p_Sc46_to_He4_Ca43) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Sc46_to_He4_Ca43);
        rate_eval.dscreened_rates_dT(k_p_Sc46_to_He4_Ca43) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 21.0_rt, 46.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Sc46_to_V50);
    rate_eval.screened_rates(k_He4_Sc46_to_V50) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Sc46_to_V50);
        rate_eval.dscreened_rates_dT(k_He4_Sc46_to_V50) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Sc46_to_n_V49);
    rate_eval.screened_rates(k_He4_Sc46_to_n_V49) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Sc46_to_n_V49);
        rate_eval.dscreened_rates_dT(k_He4_Sc46_to_n_V49) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Sc46_to_p_Ti49);
    rate_eval.screened_rates(k_He4_Sc46_to_p_Ti49) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Sc46_to_p_Ti49);
        rate_eval.dscreened_rates_dT(k_He4_Sc46_to_p_Ti49) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 21.0_rt, 47.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Sc47_to_Ti48);
    rate_eval.screened_rates(k_p_Sc47_to_Ti48) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Sc47_to_Ti48);
        rate_eval.dscreened_rates_dT(k_p_Sc47_to_Ti48) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Sc47_to_n_Ti47);
    rate_eval.screened_rates(k_p_Sc47_to_n_Ti47) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Sc47_to_n_Ti47);
        rate_eval.dscreened_rates_dT(k_p_Sc47_to_n_Ti47) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Sc47_to_He4_Ca44);
    rate_eval.screened_rates(k_p_Sc47_to_He4_Ca44) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Sc47_to_He4_Ca44);
        rate_eval.dscreened_rates_dT(k_p_Sc47_to_He4_Ca44) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 21.0_rt, 47.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Sc47_to_V51);
    rate_eval.screened_rates(k_He4_Sc47_to_V51) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Sc47_to_V51);
        rate_eval.dscreened_rates_dT(k_He4_Sc47_to_V51) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Sc47_to_n_V50);
    rate_eval.screened_rates(k_He4_Sc47_to_n_V50) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Sc47_to_n_V50);
        rate_eval.dscreened_rates_dT(k_He4_Sc47_to_n_V50) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Sc47_to_p_Ti50);
    rate_eval.screened_rates(k_He4_Sc47_to_p_Ti50) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Sc47_to_p_Ti50);
        rate_eval.dscreened_rates_dT(k_He4_Sc47_to_p_Ti50) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 21.0_rt, 48.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Sc48_to_Ti49);
    rate_eval.screened_rates(k_p_Sc48_to_Ti49) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Sc48_to_Ti49);
        rate_eval.dscreened_rates_dT(k_p_Sc48_to_Ti49) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Sc48_to_n_Ti48);
    rate_eval.screened_rates(k_p_Sc48_to_n_Ti48) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Sc48_to_n_Ti48);
        rate_eval.dscreened_rates_dT(k_p_Sc48_to_n_Ti48) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Sc48_to_He4_Ca45);
    rate_eval.screened_rates(k_p_Sc48_to_He4_Ca45) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Sc48_to_He4_Ca45);
        rate_eval.dscreened_rates_dT(k_p_Sc48_to_He4_Ca45) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 21.0_rt, 48.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Sc48_to_V52);
    rate_eval.screened_rates(k_He4_Sc48_to_V52) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Sc48_to_V52);
        rate_eval.dscreened_rates_dT(k_He4_Sc48_to_V52) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Sc48_to_n_V51);
    rate_eval.screened_rates(k_He4_Sc48_to_n_V51) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Sc48_to_n_V51);
        rate_eval.dscreened_rates_dT(k_He4_Sc48_to_n_V51) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Sc48_to_p_Ti51);
    rate_eval.screened_rates(k_He4_Sc48_to_p_Ti51) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Sc48_to_p_Ti51);
        rate_eval.dscreened_rates_dT(k_He4_Sc48_to_p_Ti51) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 21.0_rt, 49.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Sc49_to_Ti50);
    rate_eval.screened_rates(k_p_Sc49_to_Ti50) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Sc49_to_Ti50);
        rate_eval.dscreened_rates_dT(k_p_Sc49_to_Ti50) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Sc49_to_n_Ti49);
    rate_eval.screened_rates(k_p_Sc49_to_n_Ti49) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Sc49_to_n_Ti49);
        rate_eval.dscreened_rates_dT(k_p_Sc49_to_n_Ti49) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Sc49_to_He4_Ca46);
    rate_eval.screened_rates(k_p_Sc49_to_He4_Ca46) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Sc49_to_He4_Ca46);
        rate_eval.dscreened_rates_dT(k_p_Sc49_to_He4_Ca46) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 22.0_rt, 44.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ti44_to_Cr48);
    rate_eval.screened_rates(k_He4_Ti44_to_Cr48) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ti44_to_Cr48);
        rate_eval.dscreened_rates_dT(k_He4_Ti44_to_Cr48) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ti44_to_p_V47);
    rate_eval.screened_rates(k_He4_Ti44_to_p_V47) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ti44_to_p_V47);
        rate_eval.dscreened_rates_dT(k_He4_Ti44_to_p_V47) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 22.0_rt, 45.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ti45_to_V46);
    rate_eval.screened_rates(k_p_Ti45_to_V46) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ti45_to_V46);
        rate_eval.dscreened_rates_dT(k_p_Ti45_to_V46) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 22.0_rt, 45.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ti45_to_Cr49);
    rate_eval.screened_rates(k_He4_Ti45_to_Cr49) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ti45_to_Cr49);
        rate_eval.dscreened_rates_dT(k_He4_Ti45_to_Cr49) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ti45_to_n_Cr48);
    rate_eval.screened_rates(k_He4_Ti45_to_n_Cr48) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ti45_to_n_Cr48);
        rate_eval.dscreened_rates_dT(k_He4_Ti45_to_n_Cr48) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ti45_to_p_V48);
    rate_eval.screened_rates(k_He4_Ti45_to_p_V48) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ti45_to_p_V48);
        rate_eval.dscreened_rates_dT(k_He4_Ti45_to_p_V48) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 22.0_rt, 46.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ti46_to_V47);
    rate_eval.screened_rates(k_p_Ti46_to_V47) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ti46_to_V47);
        rate_eval.dscreened_rates_dT(k_p_Ti46_to_V47) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ti46_to_n_V46);
    rate_eval.screened_rates(k_p_Ti46_to_n_V46) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ti46_to_n_V46);
        rate_eval.dscreened_rates_dT(k_p_Ti46_to_n_V46) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ti46_to_He4_Sc43);
    rate_eval.screened_rates(k_p_Ti46_to_He4_Sc43) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ti46_to_He4_Sc43);
        rate_eval.dscreened_rates_dT(k_p_Ti46_to_He4_Sc43) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 22.0_rt, 46.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ti46_to_Cr50);
    rate_eval.screened_rates(k_He4_Ti46_to_Cr50) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ti46_to_Cr50);
        rate_eval.dscreened_rates_dT(k_He4_Ti46_to_Cr50) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ti46_to_n_Cr49);
    rate_eval.screened_rates(k_He4_Ti46_to_n_Cr49) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ti46_to_n_Cr49);
        rate_eval.dscreened_rates_dT(k_He4_Ti46_to_n_Cr49) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ti46_to_p_V49);
    rate_eval.screened_rates(k_He4_Ti46_to_p_V49) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ti46_to_p_V49);
        rate_eval.dscreened_rates_dT(k_He4_Ti46_to_p_V49) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 22.0_rt, 47.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ti47_to_V48);
    rate_eval.screened_rates(k_p_Ti47_to_V48) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ti47_to_V48);
        rate_eval.dscreened_rates_dT(k_p_Ti47_to_V48) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ti47_to_n_V47);
    rate_eval.screened_rates(k_p_Ti47_to_n_V47) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ti47_to_n_V47);
        rate_eval.dscreened_rates_dT(k_p_Ti47_to_n_V47) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ti47_to_He4_Sc44);
    rate_eval.screened_rates(k_p_Ti47_to_He4_Sc44) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ti47_to_He4_Sc44);
        rate_eval.dscreened_rates_dT(k_p_Ti47_to_He4_Sc44) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 22.0_rt, 47.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ti47_to_Cr51);
    rate_eval.screened_rates(k_He4_Ti47_to_Cr51) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ti47_to_Cr51);
        rate_eval.dscreened_rates_dT(k_He4_Ti47_to_Cr51) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ti47_to_n_Cr50);
    rate_eval.screened_rates(k_He4_Ti47_to_n_Cr50) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ti47_to_n_Cr50);
        rate_eval.dscreened_rates_dT(k_He4_Ti47_to_n_Cr50) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ti47_to_p_V50);
    rate_eval.screened_rates(k_He4_Ti47_to_p_V50) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ti47_to_p_V50);
        rate_eval.dscreened_rates_dT(k_He4_Ti47_to_p_V50) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 22.0_rt, 48.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ti48_to_V49);
    rate_eval.screened_rates(k_p_Ti48_to_V49) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ti48_to_V49);
        rate_eval.dscreened_rates_dT(k_p_Ti48_to_V49) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ti48_to_n_V48);
    rate_eval.screened_rates(k_p_Ti48_to_n_V48) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ti48_to_n_V48);
        rate_eval.dscreened_rates_dT(k_p_Ti48_to_n_V48) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ti48_to_He4_Sc45);
    rate_eval.screened_rates(k_p_Ti48_to_He4_Sc45) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ti48_to_He4_Sc45);
        rate_eval.dscreened_rates_dT(k_p_Ti48_to_He4_Sc45) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 22.0_rt, 48.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ti48_to_Cr52);
    rate_eval.screened_rates(k_He4_Ti48_to_Cr52) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ti48_to_Cr52);
        rate_eval.dscreened_rates_dT(k_He4_Ti48_to_Cr52) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ti48_to_n_Cr51);
    rate_eval.screened_rates(k_He4_Ti48_to_n_Cr51) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ti48_to_n_Cr51);
        rate_eval.dscreened_rates_dT(k_He4_Ti48_to_n_Cr51) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ti48_to_p_V51);
    rate_eval.screened_rates(k_He4_Ti48_to_p_V51) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ti48_to_p_V51);
        rate_eval.dscreened_rates_dT(k_He4_Ti48_to_p_V51) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 22.0_rt, 49.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ti49_to_V50);
    rate_eval.screened_rates(k_p_Ti49_to_V50) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ti49_to_V50);
        rate_eval.dscreened_rates_dT(k_p_Ti49_to_V50) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ti49_to_n_V49);
    rate_eval.screened_rates(k_p_Ti49_to_n_V49) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ti49_to_n_V49);
        rate_eval.dscreened_rates_dT(k_p_Ti49_to_n_V49) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ti49_to_He4_Sc46);
    rate_eval.screened_rates(k_p_Ti49_to_He4_Sc46) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ti49_to_He4_Sc46);
        rate_eval.dscreened_rates_dT(k_p_Ti49_to_He4_Sc46) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 22.0_rt, 49.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ti49_to_Cr53);
    rate_eval.screened_rates(k_He4_Ti49_to_Cr53) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ti49_to_Cr53);
        rate_eval.dscreened_rates_dT(k_He4_Ti49_to_Cr53) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ti49_to_n_Cr52);
    rate_eval.screened_rates(k_He4_Ti49_to_n_Cr52) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ti49_to_n_Cr52);
        rate_eval.dscreened_rates_dT(k_He4_Ti49_to_n_Cr52) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ti49_to_p_V52);
    rate_eval.screened_rates(k_He4_Ti49_to_p_V52) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ti49_to_p_V52);
        rate_eval.dscreened_rates_dT(k_He4_Ti49_to_p_V52) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 22.0_rt, 50.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ti50_to_V51);
    rate_eval.screened_rates(k_p_Ti50_to_V51) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ti50_to_V51);
        rate_eval.dscreened_rates_dT(k_p_Ti50_to_V51) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ti50_to_n_V50);
    rate_eval.screened_rates(k_p_Ti50_to_n_V50) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ti50_to_n_V50);
        rate_eval.dscreened_rates_dT(k_p_Ti50_to_n_V50) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ti50_to_He4_Sc47);
    rate_eval.screened_rates(k_p_Ti50_to_He4_Sc47) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ti50_to_He4_Sc47);
        rate_eval.dscreened_rates_dT(k_p_Ti50_to_He4_Sc47) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 22.0_rt, 50.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ti50_to_Cr54);
    rate_eval.screened_rates(k_He4_Ti50_to_Cr54) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ti50_to_Cr54);
        rate_eval.dscreened_rates_dT(k_He4_Ti50_to_Cr54) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ti50_to_n_Cr53);
    rate_eval.screened_rates(k_He4_Ti50_to_n_Cr53) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ti50_to_n_Cr53);
        rate_eval.dscreened_rates_dT(k_He4_Ti50_to_n_Cr53) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 22.0_rt, 51.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ti51_to_V52);
    rate_eval.screened_rates(k_p_Ti51_to_V52) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ti51_to_V52);
        rate_eval.dscreened_rates_dT(k_p_Ti51_to_V52) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ti51_to_n_V51);
    rate_eval.screened_rates(k_p_Ti51_to_n_V51) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ti51_to_n_V51);
        rate_eval.dscreened_rates_dT(k_p_Ti51_to_n_V51) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ti51_to_He4_Sc48);
    rate_eval.screened_rates(k_p_Ti51_to_He4_Sc48) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ti51_to_He4_Sc48);
        rate_eval.dscreened_rates_dT(k_p_Ti51_to_He4_Sc48) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 23.0_rt, 46.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_V46_to_Mn50);
    rate_eval.screened_rates(k_He4_V46_to_Mn50) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_V46_to_Mn50);
        rate_eval.dscreened_rates_dT(k_He4_V46_to_Mn50) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_V46_to_p_Cr49);
    rate_eval.screened_rates(k_He4_V46_to_p_Cr49) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_V46_to_p_Cr49);
        rate_eval.dscreened_rates_dT(k_He4_V46_to_p_Cr49) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 23.0_rt, 47.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_V47_to_Cr48);
    rate_eval.screened_rates(k_p_V47_to_Cr48) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_V47_to_Cr48);
        rate_eval.dscreened_rates_dT(k_p_V47_to_Cr48) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_V47_to_He4_Ti44);
    rate_eval.screened_rates(k_p_V47_to_He4_Ti44) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_V47_to_He4_Ti44);
        rate_eval.dscreened_rates_dT(k_p_V47_to_He4_Ti44) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 23.0_rt, 47.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_V47_to_Mn51);
    rate_eval.screened_rates(k_He4_V47_to_Mn51) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_V47_to_Mn51);
        rate_eval.dscreened_rates_dT(k_He4_V47_to_Mn51) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_V47_to_n_Mn50);
    rate_eval.screened_rates(k_He4_V47_to_n_Mn50) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_V47_to_n_Mn50);
        rate_eval.dscreened_rates_dT(k_He4_V47_to_n_Mn50) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_V47_to_p_Cr50);
    rate_eval.screened_rates(k_He4_V47_to_p_Cr50) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_V47_to_p_Cr50);
        rate_eval.dscreened_rates_dT(k_He4_V47_to_p_Cr50) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 23.0_rt, 48.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_V48_to_Cr49);
    rate_eval.screened_rates(k_p_V48_to_Cr49) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_V48_to_Cr49);
        rate_eval.dscreened_rates_dT(k_p_V48_to_Cr49) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_V48_to_n_Cr48);
    rate_eval.screened_rates(k_p_V48_to_n_Cr48) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_V48_to_n_Cr48);
        rate_eval.dscreened_rates_dT(k_p_V48_to_n_Cr48) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_V48_to_He4_Ti45);
    rate_eval.screened_rates(k_p_V48_to_He4_Ti45) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_V48_to_He4_Ti45);
        rate_eval.dscreened_rates_dT(k_p_V48_to_He4_Ti45) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 23.0_rt, 48.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_V48_to_Mn52);
    rate_eval.screened_rates(k_He4_V48_to_Mn52) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_V48_to_Mn52);
        rate_eval.dscreened_rates_dT(k_He4_V48_to_Mn52) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_V48_to_n_Mn51);
    rate_eval.screened_rates(k_He4_V48_to_n_Mn51) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_V48_to_n_Mn51);
        rate_eval.dscreened_rates_dT(k_He4_V48_to_n_Mn51) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_V48_to_p_Cr51);
    rate_eval.screened_rates(k_He4_V48_to_p_Cr51) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_V48_to_p_Cr51);
        rate_eval.dscreened_rates_dT(k_He4_V48_to_p_Cr51) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 23.0_rt, 49.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_V49_to_Cr50);
    rate_eval.screened_rates(k_p_V49_to_Cr50) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_V49_to_Cr50);
        rate_eval.dscreened_rates_dT(k_p_V49_to_Cr50) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_V49_to_n_Cr49);
    rate_eval.screened_rates(k_p_V49_to_n_Cr49) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_V49_to_n_Cr49);
        rate_eval.dscreened_rates_dT(k_p_V49_to_n_Cr49) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_V49_to_He4_Ti46);
    rate_eval.screened_rates(k_p_V49_to_He4_Ti46) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_V49_to_He4_Ti46);
        rate_eval.dscreened_rates_dT(k_p_V49_to_He4_Ti46) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 23.0_rt, 49.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_V49_to_Mn53);
    rate_eval.screened_rates(k_He4_V49_to_Mn53) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_V49_to_Mn53);
        rate_eval.dscreened_rates_dT(k_He4_V49_to_Mn53) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_V49_to_n_Mn52);
    rate_eval.screened_rates(k_He4_V49_to_n_Mn52) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_V49_to_n_Mn52);
        rate_eval.dscreened_rates_dT(k_He4_V49_to_n_Mn52) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_V49_to_p_Cr52);
    rate_eval.screened_rates(k_He4_V49_to_p_Cr52) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_V49_to_p_Cr52);
        rate_eval.dscreened_rates_dT(k_He4_V49_to_p_Cr52) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 23.0_rt, 50.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_V50_to_Cr51);
    rate_eval.screened_rates(k_p_V50_to_Cr51) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_V50_to_Cr51);
        rate_eval.dscreened_rates_dT(k_p_V50_to_Cr51) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_V50_to_n_Cr50);
    rate_eval.screened_rates(k_p_V50_to_n_Cr50) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_V50_to_n_Cr50);
        rate_eval.dscreened_rates_dT(k_p_V50_to_n_Cr50) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_V50_to_He4_Ti47);
    rate_eval.screened_rates(k_p_V50_to_He4_Ti47) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_V50_to_He4_Ti47);
        rate_eval.dscreened_rates_dT(k_p_V50_to_He4_Ti47) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 23.0_rt, 50.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_V50_to_Mn54);
    rate_eval.screened_rates(k_He4_V50_to_Mn54) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_V50_to_Mn54);
        rate_eval.dscreened_rates_dT(k_He4_V50_to_Mn54) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_V50_to_n_Mn53);
    rate_eval.screened_rates(k_He4_V50_to_n_Mn53) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_V50_to_n_Mn53);
        rate_eval.dscreened_rates_dT(k_He4_V50_to_n_Mn53) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_V50_to_p_Cr53);
    rate_eval.screened_rates(k_He4_V50_to_p_Cr53) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_V50_to_p_Cr53);
        rate_eval.dscreened_rates_dT(k_He4_V50_to_p_Cr53) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 23.0_rt, 51.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_V51_to_Cr52);
    rate_eval.screened_rates(k_p_V51_to_Cr52) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_V51_to_Cr52);
        rate_eval.dscreened_rates_dT(k_p_V51_to_Cr52) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_V51_to_n_Cr51);
    rate_eval.screened_rates(k_p_V51_to_n_Cr51) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_V51_to_n_Cr51);
        rate_eval.dscreened_rates_dT(k_p_V51_to_n_Cr51) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_V51_to_He4_Ti48);
    rate_eval.screened_rates(k_p_V51_to_He4_Ti48) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_V51_to_He4_Ti48);
        rate_eval.dscreened_rates_dT(k_p_V51_to_He4_Ti48) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 23.0_rt, 51.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_V51_to_Mn55);
    rate_eval.screened_rates(k_He4_V51_to_Mn55) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_V51_to_Mn55);
        rate_eval.dscreened_rates_dT(k_He4_V51_to_Mn55) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_V51_to_n_Mn54);
    rate_eval.screened_rates(k_He4_V51_to_n_Mn54) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_V51_to_n_Mn54);
        rate_eval.dscreened_rates_dT(k_He4_V51_to_n_Mn54) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_V51_to_p_Cr54);
    rate_eval.screened_rates(k_He4_V51_to_p_Cr54) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_V51_to_p_Cr54);
        rate_eval.dscreened_rates_dT(k_He4_V51_to_p_Cr54) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 23.0_rt, 52.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_V52_to_Cr53);
    rate_eval.screened_rates(k_p_V52_to_Cr53) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_V52_to_Cr53);
        rate_eval.dscreened_rates_dT(k_p_V52_to_Cr53) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_V52_to_n_Cr52);
    rate_eval.screened_rates(k_p_V52_to_n_Cr52) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_V52_to_n_Cr52);
        rate_eval.dscreened_rates_dT(k_p_V52_to_n_Cr52) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_V52_to_He4_Ti49);
    rate_eval.screened_rates(k_p_V52_to_He4_Ti49) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_V52_to_He4_Ti49);
        rate_eval.dscreened_rates_dT(k_p_V52_to_He4_Ti49) = ratraw * dscor_dt + dratraw_dT * scor;
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

    ratraw = rate_eval.screened_rates(k_He4_Cr48_to_p_Mn51);
    rate_eval.screened_rates(k_He4_Cr48_to_p_Mn51) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cr48_to_p_Mn51);
        rate_eval.dscreened_rates_dT(k_He4_Cr48_to_p_Mn51) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 24.0_rt, 49.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Cr49_to_Mn50);
    rate_eval.screened_rates(k_p_Cr49_to_Mn50) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cr49_to_Mn50);
        rate_eval.dscreened_rates_dT(k_p_Cr49_to_Mn50) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cr49_to_He4_V46);
    rate_eval.screened_rates(k_p_Cr49_to_He4_V46) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cr49_to_He4_V46);
        rate_eval.dscreened_rates_dT(k_p_Cr49_to_He4_V46) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 24.0_rt, 49.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Cr49_to_Fe53);
    rate_eval.screened_rates(k_He4_Cr49_to_Fe53) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cr49_to_Fe53);
        rate_eval.dscreened_rates_dT(k_He4_Cr49_to_Fe53) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cr49_to_n_Fe52);
    rate_eval.screened_rates(k_He4_Cr49_to_n_Fe52) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cr49_to_n_Fe52);
        rate_eval.dscreened_rates_dT(k_He4_Cr49_to_n_Fe52) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cr49_to_p_Mn52);
    rate_eval.screened_rates(k_He4_Cr49_to_p_Mn52) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cr49_to_p_Mn52);
        rate_eval.dscreened_rates_dT(k_He4_Cr49_to_p_Mn52) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 24.0_rt, 50.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Cr50_to_Mn51);
    rate_eval.screened_rates(k_p_Cr50_to_Mn51) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cr50_to_Mn51);
        rate_eval.dscreened_rates_dT(k_p_Cr50_to_Mn51) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cr50_to_n_Mn50);
    rate_eval.screened_rates(k_p_Cr50_to_n_Mn50) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cr50_to_n_Mn50);
        rate_eval.dscreened_rates_dT(k_p_Cr50_to_n_Mn50) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cr50_to_He4_V47);
    rate_eval.screened_rates(k_p_Cr50_to_He4_V47) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cr50_to_He4_V47);
        rate_eval.dscreened_rates_dT(k_p_Cr50_to_He4_V47) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 24.0_rt, 50.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Cr50_to_Fe54);
    rate_eval.screened_rates(k_He4_Cr50_to_Fe54) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cr50_to_Fe54);
        rate_eval.dscreened_rates_dT(k_He4_Cr50_to_Fe54) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cr50_to_n_Fe53);
    rate_eval.screened_rates(k_He4_Cr50_to_n_Fe53) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cr50_to_n_Fe53);
        rate_eval.dscreened_rates_dT(k_He4_Cr50_to_n_Fe53) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cr50_to_p_Mn53);
    rate_eval.screened_rates(k_He4_Cr50_to_p_Mn53) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cr50_to_p_Mn53);
        rate_eval.dscreened_rates_dT(k_He4_Cr50_to_p_Mn53) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 24.0_rt, 51.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Cr51_to_Mn52);
    rate_eval.screened_rates(k_p_Cr51_to_Mn52) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cr51_to_Mn52);
        rate_eval.dscreened_rates_dT(k_p_Cr51_to_Mn52) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cr51_to_n_Mn51);
    rate_eval.screened_rates(k_p_Cr51_to_n_Mn51) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cr51_to_n_Mn51);
        rate_eval.dscreened_rates_dT(k_p_Cr51_to_n_Mn51) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cr51_to_He4_V48);
    rate_eval.screened_rates(k_p_Cr51_to_He4_V48) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cr51_to_He4_V48);
        rate_eval.dscreened_rates_dT(k_p_Cr51_to_He4_V48) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 24.0_rt, 51.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Cr51_to_Fe55);
    rate_eval.screened_rates(k_He4_Cr51_to_Fe55) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cr51_to_Fe55);
        rate_eval.dscreened_rates_dT(k_He4_Cr51_to_Fe55) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cr51_to_n_Fe54);
    rate_eval.screened_rates(k_He4_Cr51_to_n_Fe54) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cr51_to_n_Fe54);
        rate_eval.dscreened_rates_dT(k_He4_Cr51_to_n_Fe54) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cr51_to_p_Mn54);
    rate_eval.screened_rates(k_He4_Cr51_to_p_Mn54) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cr51_to_p_Mn54);
        rate_eval.dscreened_rates_dT(k_He4_Cr51_to_p_Mn54) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 24.0_rt, 52.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Cr52_to_Mn53);
    rate_eval.screened_rates(k_p_Cr52_to_Mn53) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cr52_to_Mn53);
        rate_eval.dscreened_rates_dT(k_p_Cr52_to_Mn53) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cr52_to_n_Mn52);
    rate_eval.screened_rates(k_p_Cr52_to_n_Mn52) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cr52_to_n_Mn52);
        rate_eval.dscreened_rates_dT(k_p_Cr52_to_n_Mn52) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cr52_to_He4_V49);
    rate_eval.screened_rates(k_p_Cr52_to_He4_V49) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cr52_to_He4_V49);
        rate_eval.dscreened_rates_dT(k_p_Cr52_to_He4_V49) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 24.0_rt, 52.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Cr52_to_Fe56);
    rate_eval.screened_rates(k_He4_Cr52_to_Fe56) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cr52_to_Fe56);
        rate_eval.dscreened_rates_dT(k_He4_Cr52_to_Fe56) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cr52_to_n_Fe55);
    rate_eval.screened_rates(k_He4_Cr52_to_n_Fe55) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cr52_to_n_Fe55);
        rate_eval.dscreened_rates_dT(k_He4_Cr52_to_n_Fe55) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cr52_to_p_Mn55);
    rate_eval.screened_rates(k_He4_Cr52_to_p_Mn55) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cr52_to_p_Mn55);
        rate_eval.dscreened_rates_dT(k_He4_Cr52_to_p_Mn55) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 24.0_rt, 53.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Cr53_to_Mn54);
    rate_eval.screened_rates(k_p_Cr53_to_Mn54) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cr53_to_Mn54);
        rate_eval.dscreened_rates_dT(k_p_Cr53_to_Mn54) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cr53_to_n_Mn53);
    rate_eval.screened_rates(k_p_Cr53_to_n_Mn53) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cr53_to_n_Mn53);
        rate_eval.dscreened_rates_dT(k_p_Cr53_to_n_Mn53) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cr53_to_He4_V50);
    rate_eval.screened_rates(k_p_Cr53_to_He4_V50) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cr53_to_He4_V50);
        rate_eval.dscreened_rates_dT(k_p_Cr53_to_He4_V50) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 24.0_rt, 53.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Cr53_to_Fe57);
    rate_eval.screened_rates(k_He4_Cr53_to_Fe57) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cr53_to_Fe57);
        rate_eval.dscreened_rates_dT(k_He4_Cr53_to_Fe57) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cr53_to_n_Fe56);
    rate_eval.screened_rates(k_He4_Cr53_to_n_Fe56) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cr53_to_n_Fe56);
        rate_eval.dscreened_rates_dT(k_He4_Cr53_to_n_Fe56) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 24.0_rt, 54.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Cr54_to_Mn55);
    rate_eval.screened_rates(k_p_Cr54_to_Mn55) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cr54_to_Mn55);
        rate_eval.dscreened_rates_dT(k_p_Cr54_to_Mn55) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cr54_to_n_Mn54);
    rate_eval.screened_rates(k_p_Cr54_to_n_Mn54) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cr54_to_n_Mn54);
        rate_eval.dscreened_rates_dT(k_p_Cr54_to_n_Mn54) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cr54_to_He4_V51);
    rate_eval.screened_rates(k_p_Cr54_to_He4_V51) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cr54_to_He4_V51);
        rate_eval.dscreened_rates_dT(k_p_Cr54_to_He4_V51) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 24.0_rt, 54.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Cr54_to_Fe58);
    rate_eval.screened_rates(k_He4_Cr54_to_Fe58) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cr54_to_Fe58);
        rate_eval.dscreened_rates_dT(k_He4_Cr54_to_Fe58) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cr54_to_n_Fe57);
    rate_eval.screened_rates(k_He4_Cr54_to_n_Fe57) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cr54_to_n_Fe57);
        rate_eval.dscreened_rates_dT(k_He4_Cr54_to_n_Fe57) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 25.0_rt, 50.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Mn50_to_Co54);
    rate_eval.screened_rates(k_He4_Mn50_to_Co54) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mn50_to_Co54);
        rate_eval.dscreened_rates_dT(k_He4_Mn50_to_Co54) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Mn50_to_n_Co53);
    rate_eval.screened_rates(k_He4_Mn50_to_n_Co53) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mn50_to_n_Co53);
        rate_eval.dscreened_rates_dT(k_He4_Mn50_to_n_Co53) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Mn50_to_p_Fe53);
    rate_eval.screened_rates(k_He4_Mn50_to_p_Fe53) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mn50_to_p_Fe53);
        rate_eval.dscreened_rates_dT(k_He4_Mn50_to_p_Fe53) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 25.0_rt, 51.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Mn51_to_Fe52);
    rate_eval.screened_rates(k_p_Mn51_to_Fe52) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mn51_to_Fe52);
        rate_eval.dscreened_rates_dT(k_p_Mn51_to_Fe52) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Mn51_to_He4_Cr48);
    rate_eval.screened_rates(k_p_Mn51_to_He4_Cr48) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mn51_to_He4_Cr48);
        rate_eval.dscreened_rates_dT(k_p_Mn51_to_He4_Cr48) = ratraw * dscor_dt + dratraw_dT * scor;
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

    ratraw = rate_eval.screened_rates(k_He4_Mn51_to_n_Co54);
    rate_eval.screened_rates(k_He4_Mn51_to_n_Co54) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mn51_to_n_Co54);
        rate_eval.dscreened_rates_dT(k_He4_Mn51_to_n_Co54) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Mn51_to_p_Fe54);
    rate_eval.screened_rates(k_He4_Mn51_to_p_Fe54) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mn51_to_p_Fe54);
        rate_eval.dscreened_rates_dT(k_He4_Mn51_to_p_Fe54) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 25.0_rt, 52.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Mn52_to_Fe53);
    rate_eval.screened_rates(k_p_Mn52_to_Fe53) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mn52_to_Fe53);
        rate_eval.dscreened_rates_dT(k_p_Mn52_to_Fe53) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Mn52_to_n_Fe52);
    rate_eval.screened_rates(k_p_Mn52_to_n_Fe52) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mn52_to_n_Fe52);
        rate_eval.dscreened_rates_dT(k_p_Mn52_to_n_Fe52) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Mn52_to_He4_Cr49);
    rate_eval.screened_rates(k_p_Mn52_to_He4_Cr49) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mn52_to_He4_Cr49);
        rate_eval.dscreened_rates_dT(k_p_Mn52_to_He4_Cr49) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 25.0_rt, 52.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Mn52_to_Co56);
    rate_eval.screened_rates(k_He4_Mn52_to_Co56) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mn52_to_Co56);
        rate_eval.dscreened_rates_dT(k_He4_Mn52_to_Co56) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Mn52_to_n_Co55);
    rate_eval.screened_rates(k_He4_Mn52_to_n_Co55) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mn52_to_n_Co55);
        rate_eval.dscreened_rates_dT(k_He4_Mn52_to_n_Co55) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Mn52_to_p_Fe55);
    rate_eval.screened_rates(k_He4_Mn52_to_p_Fe55) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mn52_to_p_Fe55);
        rate_eval.dscreened_rates_dT(k_He4_Mn52_to_p_Fe55) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 25.0_rt, 53.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Mn53_to_Fe54);
    rate_eval.screened_rates(k_p_Mn53_to_Fe54) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mn53_to_Fe54);
        rate_eval.dscreened_rates_dT(k_p_Mn53_to_Fe54) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Mn53_to_n_Fe53);
    rate_eval.screened_rates(k_p_Mn53_to_n_Fe53) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mn53_to_n_Fe53);
        rate_eval.dscreened_rates_dT(k_p_Mn53_to_n_Fe53) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Mn53_to_He4_Cr50);
    rate_eval.screened_rates(k_p_Mn53_to_He4_Cr50) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mn53_to_He4_Cr50);
        rate_eval.dscreened_rates_dT(k_p_Mn53_to_He4_Cr50) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 25.0_rt, 53.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Mn53_to_Co57);
    rate_eval.screened_rates(k_He4_Mn53_to_Co57) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mn53_to_Co57);
        rate_eval.dscreened_rates_dT(k_He4_Mn53_to_Co57) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Mn53_to_n_Co56);
    rate_eval.screened_rates(k_He4_Mn53_to_n_Co56) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mn53_to_n_Co56);
        rate_eval.dscreened_rates_dT(k_He4_Mn53_to_n_Co56) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Mn53_to_p_Fe56);
    rate_eval.screened_rates(k_He4_Mn53_to_p_Fe56) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mn53_to_p_Fe56);
        rate_eval.dscreened_rates_dT(k_He4_Mn53_to_p_Fe56) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 25.0_rt, 54.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Mn54_to_Fe55);
    rate_eval.screened_rates(k_p_Mn54_to_Fe55) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mn54_to_Fe55);
        rate_eval.dscreened_rates_dT(k_p_Mn54_to_Fe55) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Mn54_to_n_Fe54);
    rate_eval.screened_rates(k_p_Mn54_to_n_Fe54) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mn54_to_n_Fe54);
        rate_eval.dscreened_rates_dT(k_p_Mn54_to_n_Fe54) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Mn54_to_He4_Cr51);
    rate_eval.screened_rates(k_p_Mn54_to_He4_Cr51) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mn54_to_He4_Cr51);
        rate_eval.dscreened_rates_dT(k_p_Mn54_to_He4_Cr51) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 25.0_rt, 54.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Mn54_to_Co58);
    rate_eval.screened_rates(k_He4_Mn54_to_Co58) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mn54_to_Co58);
        rate_eval.dscreened_rates_dT(k_He4_Mn54_to_Co58) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Mn54_to_n_Co57);
    rate_eval.screened_rates(k_He4_Mn54_to_n_Co57) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mn54_to_n_Co57);
        rate_eval.dscreened_rates_dT(k_He4_Mn54_to_n_Co57) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Mn54_to_p_Fe57);
    rate_eval.screened_rates(k_He4_Mn54_to_p_Fe57) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mn54_to_p_Fe57);
        rate_eval.dscreened_rates_dT(k_He4_Mn54_to_p_Fe57) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 25.0_rt, 55.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Mn55_to_Fe56);
    rate_eval.screened_rates(k_p_Mn55_to_Fe56) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mn55_to_Fe56);
        rate_eval.dscreened_rates_dT(k_p_Mn55_to_Fe56) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Mn55_to_n_Fe55);
    rate_eval.screened_rates(k_p_Mn55_to_n_Fe55) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mn55_to_n_Fe55);
        rate_eval.dscreened_rates_dT(k_p_Mn55_to_n_Fe55) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Mn55_to_He4_Cr52);
    rate_eval.screened_rates(k_p_Mn55_to_He4_Cr52) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Mn55_to_He4_Cr52);
        rate_eval.dscreened_rates_dT(k_p_Mn55_to_He4_Cr52) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 25.0_rt, 55.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Mn55_to_Co59);
    rate_eval.screened_rates(k_He4_Mn55_to_Co59) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mn55_to_Co59);
        rate_eval.dscreened_rates_dT(k_He4_Mn55_to_Co59) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Mn55_to_n_Co58);
    rate_eval.screened_rates(k_He4_Mn55_to_n_Co58) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mn55_to_n_Co58);
        rate_eval.dscreened_rates_dT(k_He4_Mn55_to_n_Co58) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Mn55_to_p_Fe58);
    rate_eval.screened_rates(k_He4_Mn55_to_p_Fe58) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mn55_to_p_Fe58);
        rate_eval.dscreened_rates_dT(k_He4_Mn55_to_p_Fe58) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 26.0_rt, 52.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Fe52_to_Co53);
    rate_eval.screened_rates(k_p_Fe52_to_Co53) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Fe52_to_Co53);
        rate_eval.dscreened_rates_dT(k_p_Fe52_to_Co53) = ratraw * dscor_dt + dratraw_dT * scor;
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

    ratraw = rate_eval.screened_rates(k_He4_Fe52_to_p_Co55);
    rate_eval.screened_rates(k_He4_Fe52_to_p_Co55) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe52_to_p_Co55);
        rate_eval.dscreened_rates_dT(k_He4_Fe52_to_p_Co55) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 26.0_rt, 53.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Fe53_to_Co54);
    rate_eval.screened_rates(k_p_Fe53_to_Co54) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Fe53_to_Co54);
        rate_eval.dscreened_rates_dT(k_p_Fe53_to_Co54) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Fe53_to_n_Co53);
    rate_eval.screened_rates(k_p_Fe53_to_n_Co53) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Fe53_to_n_Co53);
        rate_eval.dscreened_rates_dT(k_p_Fe53_to_n_Co53) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Fe53_to_He4_Mn50);
    rate_eval.screened_rates(k_p_Fe53_to_He4_Mn50) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Fe53_to_He4_Mn50);
        rate_eval.dscreened_rates_dT(k_p_Fe53_to_He4_Mn50) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 26.0_rt, 53.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Fe53_to_Ni57);
    rate_eval.screened_rates(k_He4_Fe53_to_Ni57) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe53_to_Ni57);
        rate_eval.dscreened_rates_dT(k_He4_Fe53_to_Ni57) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Fe53_to_n_Ni56);
    rate_eval.screened_rates(k_He4_Fe53_to_n_Ni56) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe53_to_n_Ni56);
        rate_eval.dscreened_rates_dT(k_He4_Fe53_to_n_Ni56) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Fe53_to_p_Co56);
    rate_eval.screened_rates(k_He4_Fe53_to_p_Co56) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe53_to_p_Co56);
        rate_eval.dscreened_rates_dT(k_He4_Fe53_to_p_Co56) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 26.0_rt, 54.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Fe54_to_Co55);
    rate_eval.screened_rates(k_p_Fe54_to_Co55) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Fe54_to_Co55);
        rate_eval.dscreened_rates_dT(k_p_Fe54_to_Co55) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Fe54_to_n_Co54);
    rate_eval.screened_rates(k_p_Fe54_to_n_Co54) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Fe54_to_n_Co54);
        rate_eval.dscreened_rates_dT(k_p_Fe54_to_n_Co54) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Fe54_to_He4_Mn51);
    rate_eval.screened_rates(k_p_Fe54_to_He4_Mn51) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Fe54_to_He4_Mn51);
        rate_eval.dscreened_rates_dT(k_p_Fe54_to_He4_Mn51) = ratraw * dscor_dt + dratraw_dT * scor;
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

    ratraw = rate_eval.screened_rates(k_He4_Fe54_to_n_Ni57);
    rate_eval.screened_rates(k_He4_Fe54_to_n_Ni57) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe54_to_n_Ni57);
        rate_eval.dscreened_rates_dT(k_He4_Fe54_to_n_Ni57) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Fe54_to_p_Co57);
    rate_eval.screened_rates(k_He4_Fe54_to_p_Co57) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe54_to_p_Co57);
        rate_eval.dscreened_rates_dT(k_He4_Fe54_to_p_Co57) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 26.0_rt, 55.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Fe55_to_Co56);
    rate_eval.screened_rates(k_p_Fe55_to_Co56) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Fe55_to_Co56);
        rate_eval.dscreened_rates_dT(k_p_Fe55_to_Co56) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Fe55_to_n_Co55);
    rate_eval.screened_rates(k_p_Fe55_to_n_Co55) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Fe55_to_n_Co55);
        rate_eval.dscreened_rates_dT(k_p_Fe55_to_n_Co55) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Fe55_to_He4_Mn52);
    rate_eval.screened_rates(k_p_Fe55_to_He4_Mn52) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Fe55_to_He4_Mn52);
        rate_eval.dscreened_rates_dT(k_p_Fe55_to_He4_Mn52) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 26.0_rt, 55.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Fe55_to_Ni59);
    rate_eval.screened_rates(k_He4_Fe55_to_Ni59) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe55_to_Ni59);
        rate_eval.dscreened_rates_dT(k_He4_Fe55_to_Ni59) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Fe55_to_n_Ni58);
    rate_eval.screened_rates(k_He4_Fe55_to_n_Ni58) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe55_to_n_Ni58);
        rate_eval.dscreened_rates_dT(k_He4_Fe55_to_n_Ni58) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Fe55_to_p_Co58);
    rate_eval.screened_rates(k_He4_Fe55_to_p_Co58) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe55_to_p_Co58);
        rate_eval.dscreened_rates_dT(k_He4_Fe55_to_p_Co58) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 26.0_rt, 56.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Fe56_to_Co57);
    rate_eval.screened_rates(k_p_Fe56_to_Co57) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Fe56_to_Co57);
        rate_eval.dscreened_rates_dT(k_p_Fe56_to_Co57) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Fe56_to_n_Co56);
    rate_eval.screened_rates(k_p_Fe56_to_n_Co56) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Fe56_to_n_Co56);
        rate_eval.dscreened_rates_dT(k_p_Fe56_to_n_Co56) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Fe56_to_He4_Mn53);
    rate_eval.screened_rates(k_p_Fe56_to_He4_Mn53) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Fe56_to_He4_Mn53);
        rate_eval.dscreened_rates_dT(k_p_Fe56_to_He4_Mn53) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 26.0_rt, 56.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Fe56_to_Ni60);
    rate_eval.screened_rates(k_He4_Fe56_to_Ni60) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe56_to_Ni60);
        rate_eval.dscreened_rates_dT(k_He4_Fe56_to_Ni60) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Fe56_to_n_Ni59);
    rate_eval.screened_rates(k_He4_Fe56_to_n_Ni59) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe56_to_n_Ni59);
        rate_eval.dscreened_rates_dT(k_He4_Fe56_to_n_Ni59) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Fe56_to_p_Co59);
    rate_eval.screened_rates(k_He4_Fe56_to_p_Co59) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe56_to_p_Co59);
        rate_eval.dscreened_rates_dT(k_He4_Fe56_to_p_Co59) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 26.0_rt, 57.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Fe57_to_Co58);
    rate_eval.screened_rates(k_p_Fe57_to_Co58) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Fe57_to_Co58);
        rate_eval.dscreened_rates_dT(k_p_Fe57_to_Co58) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Fe57_to_n_Co57);
    rate_eval.screened_rates(k_p_Fe57_to_n_Co57) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Fe57_to_n_Co57);
        rate_eval.dscreened_rates_dT(k_p_Fe57_to_n_Co57) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Fe57_to_He4_Mn54);
    rate_eval.screened_rates(k_p_Fe57_to_He4_Mn54) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Fe57_to_He4_Mn54);
        rate_eval.dscreened_rates_dT(k_p_Fe57_to_He4_Mn54) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 26.0_rt, 57.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Fe57_to_Ni61);
    rate_eval.screened_rates(k_He4_Fe57_to_Ni61) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe57_to_Ni61);
        rate_eval.dscreened_rates_dT(k_He4_Fe57_to_Ni61) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Fe57_to_n_Ni60);
    rate_eval.screened_rates(k_He4_Fe57_to_n_Ni60) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe57_to_n_Ni60);
        rate_eval.dscreened_rates_dT(k_He4_Fe57_to_n_Ni60) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 26.0_rt, 58.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Fe58_to_Co59);
    rate_eval.screened_rates(k_p_Fe58_to_Co59) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Fe58_to_Co59);
        rate_eval.dscreened_rates_dT(k_p_Fe58_to_Co59) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Fe58_to_n_Co58);
    rate_eval.screened_rates(k_p_Fe58_to_n_Co58) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Fe58_to_n_Co58);
        rate_eval.dscreened_rates_dT(k_p_Fe58_to_n_Co58) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Fe58_to_He4_Mn55);
    rate_eval.screened_rates(k_p_Fe58_to_He4_Mn55) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Fe58_to_He4_Mn55);
        rate_eval.dscreened_rates_dT(k_p_Fe58_to_He4_Mn55) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 26.0_rt, 58.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Fe58_to_Ni62);
    rate_eval.screened_rates(k_He4_Fe58_to_Ni62) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe58_to_Ni62);
        rate_eval.dscreened_rates_dT(k_He4_Fe58_to_Ni62) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Fe58_to_n_Ni61);
    rate_eval.screened_rates(k_He4_Fe58_to_n_Ni61) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Fe58_to_n_Ni61);
        rate_eval.dscreened_rates_dT(k_He4_Fe58_to_n_Ni61) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 27.0_rt, 53.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Co53_to_Cu57);
    rate_eval.screened_rates(k_He4_Co53_to_Cu57) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Co53_to_Cu57);
        rate_eval.dscreened_rates_dT(k_He4_Co53_to_Cu57) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Co53_to_p_Ni56);
    rate_eval.screened_rates(k_He4_Co53_to_p_Ni56) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Co53_to_p_Ni56);
        rate_eval.dscreened_rates_dT(k_He4_Co53_to_p_Ni56) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 27.0_rt, 54.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Co54_to_Cu58);
    rate_eval.screened_rates(k_He4_Co54_to_Cu58) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Co54_to_Cu58);
        rate_eval.dscreened_rates_dT(k_He4_Co54_to_Cu58) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Co54_to_n_Cu57);
    rate_eval.screened_rates(k_He4_Co54_to_n_Cu57) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Co54_to_n_Cu57);
        rate_eval.dscreened_rates_dT(k_He4_Co54_to_n_Cu57) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Co54_to_p_Ni57);
    rate_eval.screened_rates(k_He4_Co54_to_p_Ni57) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Co54_to_p_Ni57);
        rate_eval.dscreened_rates_dT(k_He4_Co54_to_p_Ni57) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 27.0_rt, 55.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Co55_to_Ni56);
    rate_eval.screened_rates(k_p_Co55_to_Ni56) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Co55_to_Ni56);
        rate_eval.dscreened_rates_dT(k_p_Co55_to_Ni56) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Co55_to_He4_Fe52);
    rate_eval.screened_rates(k_p_Co55_to_He4_Fe52) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Co55_to_He4_Fe52);
        rate_eval.dscreened_rates_dT(k_p_Co55_to_He4_Fe52) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 27.0_rt, 55.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Co55_to_Cu59);
    rate_eval.screened_rates(k_He4_Co55_to_Cu59) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Co55_to_Cu59);
        rate_eval.dscreened_rates_dT(k_He4_Co55_to_Cu59) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Co55_to_n_Cu58);
    rate_eval.screened_rates(k_He4_Co55_to_n_Cu58) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Co55_to_n_Cu58);
        rate_eval.dscreened_rates_dT(k_He4_Co55_to_n_Cu58) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Co55_to_p_Ni58);
    rate_eval.screened_rates(k_He4_Co55_to_p_Ni58) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Co55_to_p_Ni58);
        rate_eval.dscreened_rates_dT(k_He4_Co55_to_p_Ni58) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 27.0_rt, 56.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Co56_to_Ni57);
    rate_eval.screened_rates(k_p_Co56_to_Ni57) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Co56_to_Ni57);
        rate_eval.dscreened_rates_dT(k_p_Co56_to_Ni57) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Co56_to_n_Ni56);
    rate_eval.screened_rates(k_p_Co56_to_n_Ni56) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Co56_to_n_Ni56);
        rate_eval.dscreened_rates_dT(k_p_Co56_to_n_Ni56) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Co56_to_He4_Fe53);
    rate_eval.screened_rates(k_p_Co56_to_He4_Fe53) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Co56_to_He4_Fe53);
        rate_eval.dscreened_rates_dT(k_p_Co56_to_He4_Fe53) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 27.0_rt, 56.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Co56_to_Cu60);
    rate_eval.screened_rates(k_He4_Co56_to_Cu60) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Co56_to_Cu60);
        rate_eval.dscreened_rates_dT(k_He4_Co56_to_Cu60) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Co56_to_n_Cu59);
    rate_eval.screened_rates(k_He4_Co56_to_n_Cu59) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Co56_to_n_Cu59);
        rate_eval.dscreened_rates_dT(k_He4_Co56_to_n_Cu59) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Co56_to_p_Ni59);
    rate_eval.screened_rates(k_He4_Co56_to_p_Ni59) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Co56_to_p_Ni59);
        rate_eval.dscreened_rates_dT(k_He4_Co56_to_p_Ni59) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 27.0_rt, 57.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Co57_to_Ni58);
    rate_eval.screened_rates(k_p_Co57_to_Ni58) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Co57_to_Ni58);
        rate_eval.dscreened_rates_dT(k_p_Co57_to_Ni58) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Co57_to_n_Ni57);
    rate_eval.screened_rates(k_p_Co57_to_n_Ni57) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Co57_to_n_Ni57);
        rate_eval.dscreened_rates_dT(k_p_Co57_to_n_Ni57) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Co57_to_He4_Fe54);
    rate_eval.screened_rates(k_p_Co57_to_He4_Fe54) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Co57_to_He4_Fe54);
        rate_eval.dscreened_rates_dT(k_p_Co57_to_He4_Fe54) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 27.0_rt, 57.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Co57_to_Cu61);
    rate_eval.screened_rates(k_He4_Co57_to_Cu61) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Co57_to_Cu61);
        rate_eval.dscreened_rates_dT(k_He4_Co57_to_Cu61) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Co57_to_n_Cu60);
    rate_eval.screened_rates(k_He4_Co57_to_n_Cu60) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Co57_to_n_Cu60);
        rate_eval.dscreened_rates_dT(k_He4_Co57_to_n_Cu60) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Co57_to_p_Ni60);
    rate_eval.screened_rates(k_He4_Co57_to_p_Ni60) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Co57_to_p_Ni60);
        rate_eval.dscreened_rates_dT(k_He4_Co57_to_p_Ni60) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 27.0_rt, 58.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Co58_to_Ni59);
    rate_eval.screened_rates(k_p_Co58_to_Ni59) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Co58_to_Ni59);
        rate_eval.dscreened_rates_dT(k_p_Co58_to_Ni59) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Co58_to_n_Ni58);
    rate_eval.screened_rates(k_p_Co58_to_n_Ni58) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Co58_to_n_Ni58);
        rate_eval.dscreened_rates_dT(k_p_Co58_to_n_Ni58) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Co58_to_He4_Fe55);
    rate_eval.screened_rates(k_p_Co58_to_He4_Fe55) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Co58_to_He4_Fe55);
        rate_eval.dscreened_rates_dT(k_p_Co58_to_He4_Fe55) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 27.0_rt, 58.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Co58_to_Cu62);
    rate_eval.screened_rates(k_He4_Co58_to_Cu62) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Co58_to_Cu62);
        rate_eval.dscreened_rates_dT(k_He4_Co58_to_Cu62) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Co58_to_n_Cu61);
    rate_eval.screened_rates(k_He4_Co58_to_n_Cu61) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Co58_to_n_Cu61);
        rate_eval.dscreened_rates_dT(k_He4_Co58_to_n_Cu61) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Co58_to_p_Ni61);
    rate_eval.screened_rates(k_He4_Co58_to_p_Ni61) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Co58_to_p_Ni61);
        rate_eval.dscreened_rates_dT(k_He4_Co58_to_p_Ni61) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 27.0_rt, 59.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Co59_to_Ni60);
    rate_eval.screened_rates(k_p_Co59_to_Ni60) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Co59_to_Ni60);
        rate_eval.dscreened_rates_dT(k_p_Co59_to_Ni60) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Co59_to_n_Ni59);
    rate_eval.screened_rates(k_p_Co59_to_n_Ni59) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Co59_to_n_Ni59);
        rate_eval.dscreened_rates_dT(k_p_Co59_to_n_Ni59) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Co59_to_He4_Fe56);
    rate_eval.screened_rates(k_p_Co59_to_He4_Fe56) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Co59_to_He4_Fe56);
        rate_eval.dscreened_rates_dT(k_p_Co59_to_He4_Fe56) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 27.0_rt, 59.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Co59_to_Cu63);
    rate_eval.screened_rates(k_He4_Co59_to_Cu63) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Co59_to_Cu63);
        rate_eval.dscreened_rates_dT(k_He4_Co59_to_Cu63) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Co59_to_n_Cu62);
    rate_eval.screened_rates(k_He4_Co59_to_n_Cu62) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Co59_to_n_Cu62);
        rate_eval.dscreened_rates_dT(k_He4_Co59_to_n_Cu62) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Co59_to_p_Ni62);
    rate_eval.screened_rates(k_He4_Co59_to_p_Ni62) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Co59_to_p_Ni62);
        rate_eval.dscreened_rates_dT(k_He4_Co59_to_p_Ni62) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 28.0_rt, 56.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ni56_to_Cu57);
    rate_eval.screened_rates(k_p_Ni56_to_Cu57) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ni56_to_Cu57);
        rate_eval.dscreened_rates_dT(k_p_Ni56_to_Cu57) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ni56_to_He4_Co53);
    rate_eval.screened_rates(k_p_Ni56_to_He4_Co53) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ni56_to_He4_Co53);
        rate_eval.dscreened_rates_dT(k_p_Ni56_to_He4_Co53) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 28.0_rt, 56.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ni56_to_Zn60);
    rate_eval.screened_rates(k_He4_Ni56_to_Zn60) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ni56_to_Zn60);
        rate_eval.dscreened_rates_dT(k_He4_Ni56_to_Zn60) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ni56_to_n_Zn59);
    rate_eval.screened_rates(k_He4_Ni56_to_n_Zn59) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ni56_to_n_Zn59);
        rate_eval.dscreened_rates_dT(k_He4_Ni56_to_n_Zn59) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ni56_to_p_Cu59);
    rate_eval.screened_rates(k_He4_Ni56_to_p_Cu59) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ni56_to_p_Cu59);
        rate_eval.dscreened_rates_dT(k_He4_Ni56_to_p_Cu59) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 28.0_rt, 57.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ni57_to_Cu58);
    rate_eval.screened_rates(k_p_Ni57_to_Cu58) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ni57_to_Cu58);
        rate_eval.dscreened_rates_dT(k_p_Ni57_to_Cu58) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ni57_to_n_Cu57);
    rate_eval.screened_rates(k_p_Ni57_to_n_Cu57) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ni57_to_n_Cu57);
        rate_eval.dscreened_rates_dT(k_p_Ni57_to_n_Cu57) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ni57_to_He4_Co54);
    rate_eval.screened_rates(k_p_Ni57_to_He4_Co54) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ni57_to_He4_Co54);
        rate_eval.dscreened_rates_dT(k_p_Ni57_to_He4_Co54) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 28.0_rt, 57.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ni57_to_Zn61);
    rate_eval.screened_rates(k_He4_Ni57_to_Zn61) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ni57_to_Zn61);
        rate_eval.dscreened_rates_dT(k_He4_Ni57_to_Zn61) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ni57_to_n_Zn60);
    rate_eval.screened_rates(k_He4_Ni57_to_n_Zn60) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ni57_to_n_Zn60);
        rate_eval.dscreened_rates_dT(k_He4_Ni57_to_n_Zn60) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ni57_to_p_Cu60);
    rate_eval.screened_rates(k_He4_Ni57_to_p_Cu60) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ni57_to_p_Cu60);
        rate_eval.dscreened_rates_dT(k_He4_Ni57_to_p_Cu60) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 28.0_rt, 58.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ni58_to_Cu59);
    rate_eval.screened_rates(k_p_Ni58_to_Cu59) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ni58_to_Cu59);
        rate_eval.dscreened_rates_dT(k_p_Ni58_to_Cu59) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ni58_to_n_Cu58);
    rate_eval.screened_rates(k_p_Ni58_to_n_Cu58) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ni58_to_n_Cu58);
        rate_eval.dscreened_rates_dT(k_p_Ni58_to_n_Cu58) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ni58_to_He4_Co55);
    rate_eval.screened_rates(k_p_Ni58_to_He4_Co55) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ni58_to_He4_Co55);
        rate_eval.dscreened_rates_dT(k_p_Ni58_to_He4_Co55) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 28.0_rt, 58.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ni58_to_Zn62);
    rate_eval.screened_rates(k_He4_Ni58_to_Zn62) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ni58_to_Zn62);
        rate_eval.dscreened_rates_dT(k_He4_Ni58_to_Zn62) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ni58_to_n_Zn61);
    rate_eval.screened_rates(k_He4_Ni58_to_n_Zn61) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ni58_to_n_Zn61);
        rate_eval.dscreened_rates_dT(k_He4_Ni58_to_n_Zn61) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ni58_to_p_Cu61);
    rate_eval.screened_rates(k_He4_Ni58_to_p_Cu61) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ni58_to_p_Cu61);
        rate_eval.dscreened_rates_dT(k_He4_Ni58_to_p_Cu61) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 28.0_rt, 59.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ni59_to_Cu60);
    rate_eval.screened_rates(k_p_Ni59_to_Cu60) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ni59_to_Cu60);
        rate_eval.dscreened_rates_dT(k_p_Ni59_to_Cu60) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ni59_to_n_Cu59);
    rate_eval.screened_rates(k_p_Ni59_to_n_Cu59) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ni59_to_n_Cu59);
        rate_eval.dscreened_rates_dT(k_p_Ni59_to_n_Cu59) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ni59_to_He4_Co56);
    rate_eval.screened_rates(k_p_Ni59_to_He4_Co56) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ni59_to_He4_Co56);
        rate_eval.dscreened_rates_dT(k_p_Ni59_to_He4_Co56) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 28.0_rt, 59.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ni59_to_Zn63);
    rate_eval.screened_rates(k_He4_Ni59_to_Zn63) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ni59_to_Zn63);
        rate_eval.dscreened_rates_dT(k_He4_Ni59_to_Zn63) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ni59_to_n_Zn62);
    rate_eval.screened_rates(k_He4_Ni59_to_n_Zn62) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ni59_to_n_Zn62);
        rate_eval.dscreened_rates_dT(k_He4_Ni59_to_n_Zn62) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ni59_to_p_Cu62);
    rate_eval.screened_rates(k_He4_Ni59_to_p_Cu62) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ni59_to_p_Cu62);
        rate_eval.dscreened_rates_dT(k_He4_Ni59_to_p_Cu62) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 28.0_rt, 60.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ni60_to_Cu61);
    rate_eval.screened_rates(k_p_Ni60_to_Cu61) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ni60_to_Cu61);
        rate_eval.dscreened_rates_dT(k_p_Ni60_to_Cu61) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ni60_to_n_Cu60);
    rate_eval.screened_rates(k_p_Ni60_to_n_Cu60) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ni60_to_n_Cu60);
        rate_eval.dscreened_rates_dT(k_p_Ni60_to_n_Cu60) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ni60_to_He4_Co57);
    rate_eval.screened_rates(k_p_Ni60_to_He4_Co57) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ni60_to_He4_Co57);
        rate_eval.dscreened_rates_dT(k_p_Ni60_to_He4_Co57) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 28.0_rt, 60.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ni60_to_Zn64);
    rate_eval.screened_rates(k_He4_Ni60_to_Zn64) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ni60_to_Zn64);
        rate_eval.dscreened_rates_dT(k_He4_Ni60_to_Zn64) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ni60_to_n_Zn63);
    rate_eval.screened_rates(k_He4_Ni60_to_n_Zn63) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ni60_to_n_Zn63);
        rate_eval.dscreened_rates_dT(k_He4_Ni60_to_n_Zn63) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ni60_to_p_Cu63);
    rate_eval.screened_rates(k_He4_Ni60_to_p_Cu63) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ni60_to_p_Cu63);
        rate_eval.dscreened_rates_dT(k_He4_Ni60_to_p_Cu63) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 28.0_rt, 61.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ni61_to_Cu62);
    rate_eval.screened_rates(k_p_Ni61_to_Cu62) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ni61_to_Cu62);
        rate_eval.dscreened_rates_dT(k_p_Ni61_to_Cu62) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ni61_to_n_Cu61);
    rate_eval.screened_rates(k_p_Ni61_to_n_Cu61) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ni61_to_n_Cu61);
        rate_eval.dscreened_rates_dT(k_p_Ni61_to_n_Cu61) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ni61_to_He4_Co58);
    rate_eval.screened_rates(k_p_Ni61_to_He4_Co58) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ni61_to_He4_Co58);
        rate_eval.dscreened_rates_dT(k_p_Ni61_to_He4_Co58) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 28.0_rt, 61.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ni61_to_Zn65);
    rate_eval.screened_rates(k_He4_Ni61_to_Zn65) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ni61_to_Zn65);
        rate_eval.dscreened_rates_dT(k_He4_Ni61_to_Zn65) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ni61_to_n_Zn64);
    rate_eval.screened_rates(k_He4_Ni61_to_n_Zn64) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ni61_to_n_Zn64);
        rate_eval.dscreened_rates_dT(k_He4_Ni61_to_n_Zn64) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ni61_to_p_Cu64);
    rate_eval.screened_rates(k_He4_Ni61_to_p_Cu64) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ni61_to_p_Cu64);
        rate_eval.dscreened_rates_dT(k_He4_Ni61_to_p_Cu64) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 28.0_rt, 62.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ni62_to_Cu63);
    rate_eval.screened_rates(k_p_Ni62_to_Cu63) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ni62_to_Cu63);
        rate_eval.dscreened_rates_dT(k_p_Ni62_to_Cu63) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ni62_to_n_Cu62);
    rate_eval.screened_rates(k_p_Ni62_to_n_Cu62) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ni62_to_n_Cu62);
        rate_eval.dscreened_rates_dT(k_p_Ni62_to_n_Cu62) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ni62_to_He4_Co59);
    rate_eval.screened_rates(k_p_Ni62_to_He4_Co59) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ni62_to_He4_Co59);
        rate_eval.dscreened_rates_dT(k_p_Ni62_to_He4_Co59) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 28.0_rt, 62.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ni62_to_Zn66);
    rate_eval.screened_rates(k_He4_Ni62_to_Zn66) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ni62_to_Zn66);
        rate_eval.dscreened_rates_dT(k_He4_Ni62_to_Zn66) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ni62_to_n_Zn65);
    rate_eval.screened_rates(k_He4_Ni62_to_n_Zn65) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ni62_to_n_Zn65);
        rate_eval.dscreened_rates_dT(k_He4_Ni62_to_n_Zn65) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Ni62_to_p_Cu65);
    rate_eval.screened_rates(k_He4_Ni62_to_p_Cu65) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ni62_to_p_Cu65);
        rate_eval.dscreened_rates_dT(k_He4_Ni62_to_p_Cu65) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 28.0_rt, 63.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ni63_to_Cu64);
    rate_eval.screened_rates(k_p_Ni63_to_Cu64) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ni63_to_Cu64);
        rate_eval.dscreened_rates_dT(k_p_Ni63_to_Cu64) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ni63_to_n_Cu63);
    rate_eval.screened_rates(k_p_Ni63_to_n_Cu63) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ni63_to_n_Cu63);
        rate_eval.dscreened_rates_dT(k_p_Ni63_to_n_Cu63) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 28.0_rt, 64.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ni64_to_Cu65);
    rate_eval.screened_rates(k_p_Ni64_to_Cu65) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ni64_to_Cu65);
        rate_eval.dscreened_rates_dT(k_p_Ni64_to_Cu65) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ni64_to_n_Cu64);
    rate_eval.screened_rates(k_p_Ni64_to_n_Cu64) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ni64_to_n_Cu64);
        rate_eval.dscreened_rates_dT(k_p_Ni64_to_n_Cu64) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 29.0_rt, 58.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Cu58_to_Zn59);
    rate_eval.screened_rates(k_p_Cu58_to_Zn59) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cu58_to_Zn59);
        rate_eval.dscreened_rates_dT(k_p_Cu58_to_Zn59) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 29.0_rt, 58.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Cu58_to_Ga62);
    rate_eval.screened_rates(k_He4_Cu58_to_Ga62) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cu58_to_Ga62);
        rate_eval.dscreened_rates_dT(k_He4_Cu58_to_Ga62) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cu58_to_p_Zn61);
    rate_eval.screened_rates(k_He4_Cu58_to_p_Zn61) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cu58_to_p_Zn61);
        rate_eval.dscreened_rates_dT(k_He4_Cu58_to_p_Zn61) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 29.0_rt, 59.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Cu59_to_Zn60);
    rate_eval.screened_rates(k_p_Cu59_to_Zn60) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cu59_to_Zn60);
        rate_eval.dscreened_rates_dT(k_p_Cu59_to_Zn60) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cu59_to_n_Zn59);
    rate_eval.screened_rates(k_p_Cu59_to_n_Zn59) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cu59_to_n_Zn59);
        rate_eval.dscreened_rates_dT(k_p_Cu59_to_n_Zn59) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cu59_to_He4_Ni56);
    rate_eval.screened_rates(k_p_Cu59_to_He4_Ni56) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cu59_to_He4_Ni56);
        rate_eval.dscreened_rates_dT(k_p_Cu59_to_He4_Ni56) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 29.0_rt, 59.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Cu59_to_Ga63);
    rate_eval.screened_rates(k_He4_Cu59_to_Ga63) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cu59_to_Ga63);
        rate_eval.dscreened_rates_dT(k_He4_Cu59_to_Ga63) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cu59_to_n_Ga62);
    rate_eval.screened_rates(k_He4_Cu59_to_n_Ga62) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cu59_to_n_Ga62);
        rate_eval.dscreened_rates_dT(k_He4_Cu59_to_n_Ga62) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cu59_to_p_Zn62);
    rate_eval.screened_rates(k_He4_Cu59_to_p_Zn62) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cu59_to_p_Zn62);
        rate_eval.dscreened_rates_dT(k_He4_Cu59_to_p_Zn62) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 29.0_rt, 60.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Cu60_to_Zn61);
    rate_eval.screened_rates(k_p_Cu60_to_Zn61) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cu60_to_Zn61);
        rate_eval.dscreened_rates_dT(k_p_Cu60_to_Zn61) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cu60_to_n_Zn60);
    rate_eval.screened_rates(k_p_Cu60_to_n_Zn60) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cu60_to_n_Zn60);
        rate_eval.dscreened_rates_dT(k_p_Cu60_to_n_Zn60) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cu60_to_He4_Ni57);
    rate_eval.screened_rates(k_p_Cu60_to_He4_Ni57) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cu60_to_He4_Ni57);
        rate_eval.dscreened_rates_dT(k_p_Cu60_to_He4_Ni57) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 29.0_rt, 60.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Cu60_to_Ga64);
    rate_eval.screened_rates(k_He4_Cu60_to_Ga64) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cu60_to_Ga64);
        rate_eval.dscreened_rates_dT(k_He4_Cu60_to_Ga64) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cu60_to_n_Ga63);
    rate_eval.screened_rates(k_He4_Cu60_to_n_Ga63) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cu60_to_n_Ga63);
        rate_eval.dscreened_rates_dT(k_He4_Cu60_to_n_Ga63) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cu60_to_p_Zn63);
    rate_eval.screened_rates(k_He4_Cu60_to_p_Zn63) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cu60_to_p_Zn63);
        rate_eval.dscreened_rates_dT(k_He4_Cu60_to_p_Zn63) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 29.0_rt, 61.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Cu61_to_Zn62);
    rate_eval.screened_rates(k_p_Cu61_to_Zn62) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cu61_to_Zn62);
        rate_eval.dscreened_rates_dT(k_p_Cu61_to_Zn62) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cu61_to_n_Zn61);
    rate_eval.screened_rates(k_p_Cu61_to_n_Zn61) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cu61_to_n_Zn61);
        rate_eval.dscreened_rates_dT(k_p_Cu61_to_n_Zn61) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cu61_to_He4_Ni58);
    rate_eval.screened_rates(k_p_Cu61_to_He4_Ni58) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cu61_to_He4_Ni58);
        rate_eval.dscreened_rates_dT(k_p_Cu61_to_He4_Ni58) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 29.0_rt, 62.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Cu62_to_Zn63);
    rate_eval.screened_rates(k_p_Cu62_to_Zn63) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cu62_to_Zn63);
        rate_eval.dscreened_rates_dT(k_p_Cu62_to_Zn63) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cu62_to_n_Zn62);
    rate_eval.screened_rates(k_p_Cu62_to_n_Zn62) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cu62_to_n_Zn62);
        rate_eval.dscreened_rates_dT(k_p_Cu62_to_n_Zn62) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cu62_to_He4_Ni59);
    rate_eval.screened_rates(k_p_Cu62_to_He4_Ni59) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cu62_to_He4_Ni59);
        rate_eval.dscreened_rates_dT(k_p_Cu62_to_He4_Ni59) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 29.0_rt, 63.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Cu63_to_Zn64);
    rate_eval.screened_rates(k_p_Cu63_to_Zn64) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cu63_to_Zn64);
        rate_eval.dscreened_rates_dT(k_p_Cu63_to_Zn64) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cu63_to_n_Zn63);
    rate_eval.screened_rates(k_p_Cu63_to_n_Zn63) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cu63_to_n_Zn63);
        rate_eval.dscreened_rates_dT(k_p_Cu63_to_n_Zn63) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cu63_to_He4_Ni60);
    rate_eval.screened_rates(k_p_Cu63_to_He4_Ni60) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cu63_to_He4_Ni60);
        rate_eval.dscreened_rates_dT(k_p_Cu63_to_He4_Ni60) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 29.0_rt, 64.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Cu64_to_Zn65);
    rate_eval.screened_rates(k_p_Cu64_to_Zn65) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cu64_to_Zn65);
        rate_eval.dscreened_rates_dT(k_p_Cu64_to_Zn65) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cu64_to_n_Zn64);
    rate_eval.screened_rates(k_p_Cu64_to_n_Zn64) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cu64_to_n_Zn64);
        rate_eval.dscreened_rates_dT(k_p_Cu64_to_n_Zn64) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cu64_to_He4_Ni61);
    rate_eval.screened_rates(k_p_Cu64_to_He4_Ni61) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cu64_to_He4_Ni61);
        rate_eval.dscreened_rates_dT(k_p_Cu64_to_He4_Ni61) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 29.0_rt, 65.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Cu65_to_Zn66);
    rate_eval.screened_rates(k_p_Cu65_to_Zn66) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cu65_to_Zn66);
        rate_eval.dscreened_rates_dT(k_p_Cu65_to_Zn66) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cu65_to_n_Zn65);
    rate_eval.screened_rates(k_p_Cu65_to_n_Zn65) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cu65_to_n_Zn65);
        rate_eval.dscreened_rates_dT(k_p_Cu65_to_n_Zn65) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Cu65_to_He4_Ni62);
    rate_eval.screened_rates(k_p_Cu65_to_He4_Ni62) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Cu65_to_He4_Ni62);
        rate_eval.dscreened_rates_dT(k_p_Cu65_to_He4_Ni62) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 30.0_rt, 59.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Zn59_to_Ge63);
    rate_eval.screened_rates(k_He4_Zn59_to_Ge63) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Zn59_to_Ge63);
        rate_eval.dscreened_rates_dT(k_He4_Zn59_to_Ge63) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Zn59_to_p_Ga62);
    rate_eval.screened_rates(k_He4_Zn59_to_p_Ga62) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Zn59_to_p_Ga62);
        rate_eval.dscreened_rates_dT(k_He4_Zn59_to_p_Ga62) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 30.0_rt, 60.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Zn60_to_Ge64);
    rate_eval.screened_rates(k_He4_Zn60_to_Ge64) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Zn60_to_Ge64);
        rate_eval.dscreened_rates_dT(k_He4_Zn60_to_Ge64) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Zn60_to_n_Ge63);
    rate_eval.screened_rates(k_He4_Zn60_to_n_Ge63) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Zn60_to_n_Ge63);
        rate_eval.dscreened_rates_dT(k_He4_Zn60_to_n_Ge63) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Zn60_to_p_Ga63);
    rate_eval.screened_rates(k_He4_Zn60_to_p_Ga63) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Zn60_to_p_Ga63);
        rate_eval.dscreened_rates_dT(k_He4_Zn60_to_p_Ga63) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 30.0_rt, 61.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Zn61_to_Ga62);
    rate_eval.screened_rates(k_p_Zn61_to_Ga62) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Zn61_to_Ga62);
        rate_eval.dscreened_rates_dT(k_p_Zn61_to_Ga62) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Zn61_to_He4_Cu58);
    rate_eval.screened_rates(k_p_Zn61_to_He4_Cu58) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Zn61_to_He4_Cu58);
        rate_eval.dscreened_rates_dT(k_p_Zn61_to_He4_Cu58) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 30.0_rt, 62.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Zn62_to_Ga63);
    rate_eval.screened_rates(k_p_Zn62_to_Ga63) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Zn62_to_Ga63);
        rate_eval.dscreened_rates_dT(k_p_Zn62_to_Ga63) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Zn62_to_n_Ga62);
    rate_eval.screened_rates(k_p_Zn62_to_n_Ga62) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Zn62_to_n_Ga62);
        rate_eval.dscreened_rates_dT(k_p_Zn62_to_n_Ga62) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Zn62_to_He4_Cu59);
    rate_eval.screened_rates(k_p_Zn62_to_He4_Cu59) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Zn62_to_He4_Cu59);
        rate_eval.dscreened_rates_dT(k_p_Zn62_to_He4_Cu59) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 30.0_rt, 63.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Zn63_to_Ga64);
    rate_eval.screened_rates(k_p_Zn63_to_Ga64) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Zn63_to_Ga64);
        rate_eval.dscreened_rates_dT(k_p_Zn63_to_Ga64) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Zn63_to_n_Ga63);
    rate_eval.screened_rates(k_p_Zn63_to_n_Ga63) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Zn63_to_n_Ga63);
        rate_eval.dscreened_rates_dT(k_p_Zn63_to_n_Ga63) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Zn63_to_He4_Cu60);
    rate_eval.screened_rates(k_p_Zn63_to_He4_Cu60) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Zn63_to_He4_Cu60);
        rate_eval.dscreened_rates_dT(k_p_Zn63_to_He4_Cu60) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 31.0_rt, 62.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ga62_to_Ge63);
    rate_eval.screened_rates(k_p_Ga62_to_Ge63) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ga62_to_Ge63);
        rate_eval.dscreened_rates_dT(k_p_Ga62_to_Ge63) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ga62_to_He4_Zn59);
    rate_eval.screened_rates(k_p_Ga62_to_He4_Zn59) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ga62_to_He4_Zn59);
        rate_eval.dscreened_rates_dT(k_p_Ga62_to_He4_Zn59) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 31.0_rt, 63.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ga63_to_Ge64);
    rate_eval.screened_rates(k_p_Ga63_to_Ge64) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ga63_to_Ge64);
        rate_eval.dscreened_rates_dT(k_p_Ga63_to_Ge64) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ga63_to_n_Ge63);
    rate_eval.screened_rates(k_p_Ga63_to_n_Ge63) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ga63_to_n_Ge63);
        rate_eval.dscreened_rates_dT(k_p_Ga63_to_n_Ge63) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ga63_to_He4_Zn60);
    rate_eval.screened_rates(k_p_Ga63_to_He4_Zn60) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ga63_to_He4_Zn60);
        rate_eval.dscreened_rates_dT(k_p_Ga63_to_He4_Zn60) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 2.0_rt, 2.0_rt, 3.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_d_He3_to_p_He4);
    rate_eval.screened_rates(k_d_He3_to_p_He4) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_d_He3_to_p_He4);
        rate_eval.dscreened_rates_dT(k_d_He3_to_p_He4) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 2.0_rt, 4.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_He4_to_d_He3);
    rate_eval.screened_rates(k_p_He4_to_d_He3) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_He4_to_d_He3);
        rate_eval.dscreened_rates_dT(k_p_He4_to_d_He3) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_n_p_He4_to_Li6);
    rate_eval.screened_rates(k_n_p_He4_to_Li6) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_n_p_He4_to_Li6);
        rate_eval.dscreened_rates_dT(k_n_p_He4_to_Li6) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 2.0_rt, 4.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_He4_to_n_Be7);
    rate_eval.screened_rates(k_He4_He4_to_n_Be7) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_He4_to_n_Be7);
        rate_eval.dscreened_rates_dT(k_He4_He4_to_n_Be7) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_He4_to_p_Li7);
    rate_eval.screened_rates(k_He4_He4_to_p_Li7) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_He4_to_p_Li7);
        rate_eval.dscreened_rates_dT(k_He4_He4_to_p_Li7) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_n_He4_He4_to_Be9);
    rate_eval.screened_rates(k_n_He4_He4_to_Be9) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_n_He4_He4_to_Be9);
        rate_eval.dscreened_rates_dT(k_n_He4_He4_to_Be9) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_n_He4_He4_to_d_Li7);
    rate_eval.screened_rates(k_n_He4_He4_to_d_Li7) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_n_He4_He4_to_d_Li7);
        rate_eval.dscreened_rates_dT(k_n_He4_He4_to_d_Li7) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 2.0_rt, 3.0_rt, 6.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_d_Li6_to_n_Be7);
    rate_eval.screened_rates(k_d_Li6_to_n_Be7) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_d_Li6_to_n_Be7);
        rate_eval.dscreened_rates_dT(k_d_Li6_to_n_Be7) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_d_Li6_to_p_Li7);
    rate_eval.screened_rates(k_d_Li6_to_p_Li7) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_d_Li6_to_p_Li7);
        rate_eval.dscreened_rates_dT(k_d_Li6_to_p_Li7) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 3.0_rt, 7.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Li7_to_n_Be7);
    rate_eval.screened_rates(k_p_Li7_to_n_Be7) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Li7_to_n_Be7);
        rate_eval.dscreened_rates_dT(k_p_Li7_to_n_Be7) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Li7_to_d_Li6);
    rate_eval.screened_rates(k_p_Li7_to_d_Li6) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Li7_to_d_Li6);
        rate_eval.dscreened_rates_dT(k_p_Li7_to_d_Li6) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Li7_to_He4_He4);
    rate_eval.screened_rates(k_p_Li7_to_He4_He4) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Li7_to_He4_He4);
        rate_eval.dscreened_rates_dT(k_p_Li7_to_He4_He4) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 4.0_rt, 7.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Be7_to_p_B10);
    rate_eval.screened_rates(k_He4_Be7_to_p_B10) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Be7_to_p_B10);
        rate_eval.dscreened_rates_dT(k_He4_Be7_to_p_B10) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 4.0_rt, 9.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Be9_to_n_C12);
    rate_eval.screened_rates(k_He4_Be9_to_n_C12) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Be9_to_n_C12);
        rate_eval.dscreened_rates_dT(k_He4_Be9_to_n_C12) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 5.0_rt, 10.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_B10_to_He4_Be7);
    rate_eval.screened_rates(k_p_B10_to_He4_Be7) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_B10_to_He4_Be7);
        rate_eval.dscreened_rates_dT(k_p_B10_to_He4_Be7) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 5.0_rt, 10.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_B10_to_n_N13);
    rate_eval.screened_rates(k_He4_B10_to_n_N13) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_B10_to_n_N13);
        rate_eval.dscreened_rates_dT(k_He4_B10_to_n_N13) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_B10_to_p_C13);
    rate_eval.screened_rates(k_He4_B10_to_p_C13) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_B10_to_p_C13);
        rate_eval.dscreened_rates_dT(k_He4_B10_to_p_C13) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 5.0_rt, 11.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_B11_to_n_N14);
    rate_eval.screened_rates(k_He4_B11_to_n_N14) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_B11_to_n_N14);
        rate_eval.dscreened_rates_dT(k_He4_B11_to_n_N14) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_B11_to_p_C14);
    rate_eval.screened_rates(k_He4_B11_to_p_C14) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_B11_to_p_C14);
        rate_eval.dscreened_rates_dT(k_He4_B11_to_p_C14) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(6.0_rt, 12.0_rt, 6.0_rt, 12.0_rt);


        static_assert(scn_fac.z1 == 6.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_C12_C12_to_n_Mg23);
    rate_eval.screened_rates(k_C12_C12_to_n_Mg23) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_C12_to_n_Mg23);
        rate_eval.dscreened_rates_dT(k_C12_C12_to_n_Mg23) = ratraw * dscor_dt + dratraw_dT * scor;
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


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 2.0_rt, 6.0_rt, 13.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_d_C13_to_n_N14);
    rate_eval.screened_rates(k_d_C13_to_n_N14) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_d_C13_to_n_N14);
        rate_eval.dscreened_rates_dT(k_d_C13_to_n_N14) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 6.0_rt, 13.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_C13_to_n_O16);
    rate_eval.screened_rates(k_He4_C13_to_n_O16) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_C13_to_n_O16);
        rate_eval.dscreened_rates_dT(k_He4_C13_to_n_O16) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 2.0_rt, 6.0_rt, 14.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_d_C14_to_n_N15);
    rate_eval.screened_rates(k_d_C14_to_n_N15) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_d_C14_to_n_N15);
        rate_eval.dscreened_rates_dT(k_d_C14_to_n_N15) = ratraw * dscor_dt + dratraw_dT * scor;
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


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 10.0_rt, 18.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ne18_to_p_Na21);
    rate_eval.screened_rates(k_He4_Ne18_to_p_Na21) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ne18_to_p_Na21);
        rate_eval.dscreened_rates_dT(k_He4_Ne18_to_p_Na21) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(6.0_rt, 12.0_rt, 10.0_rt, 20.0_rt);


        static_assert(scn_fac.z1 == 6.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_C12_Ne20_to_p_P31);
    rate_eval.screened_rates(k_C12_Ne20_to_p_P31) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_Ne20_to_p_P31);
        rate_eval.dscreened_rates_dT(k_C12_Ne20_to_p_P31) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_C12_Ne20_to_He4_Si28);
    rate_eval.screened_rates(k_C12_Ne20_to_He4_Si28) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_C12_Ne20_to_He4_Si28);
        rate_eval.dscreened_rates_dT(k_C12_Ne20_to_He4_Si28) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 11.0_rt, 21.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Na21_to_He4_Ne18);
    rate_eval.screened_rates(k_p_Na21_to_He4_Ne18) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Na21_to_He4_Ne18);
        rate_eval.dscreened_rates_dT(k_p_Na21_to_He4_Ne18) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 12.0_rt, 23.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Mg23_to_p_Al26);
    rate_eval.screened_rates(k_He4_Mg23_to_p_Al26) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Mg23_to_p_Al26);
        rate_eval.dscreened_rates_dT(k_He4_Mg23_to_p_Al26) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 13.0_rt, 26.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Al26_to_He4_Mg23);
    rate_eval.screened_rates(k_p_Al26_to_He4_Mg23) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Al26_to_He4_Mg23);
        rate_eval.dscreened_rates_dT(k_p_Al26_to_He4_Mg23) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 19.0_rt, 37.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_K37_to_p_Ca40);
    rate_eval.screened_rates(k_He4_K37_to_p_Ca40) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_K37_to_p_Ca40);
        rate_eval.dscreened_rates_dT(k_He4_K37_to_p_Ca40) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 19.0_rt, 38.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_K38_to_p_Ca41);
    rate_eval.screened_rates(k_He4_K38_to_p_Ca41) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_K38_to_p_Ca41);
        rate_eval.dscreened_rates_dT(k_He4_K38_to_p_Ca41) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 20.0_rt, 40.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ca40_to_He4_K37);
    rate_eval.screened_rates(k_p_Ca40_to_He4_K37) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ca40_to_He4_K37);
        rate_eval.dscreened_rates_dT(k_p_Ca40_to_He4_K37) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 20.0_rt, 41.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ca41_to_He4_K38);
    rate_eval.screened_rates(k_p_Ca41_to_He4_K38) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ca41_to_He4_K38);
        rate_eval.dscreened_rates_dT(k_p_Ca41_to_He4_K38) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 20.0_rt, 48.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ca48_to_n_Ti51);
    rate_eval.screened_rates(k_He4_Ca48_to_n_Ti51) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ca48_to_n_Ti51);
        rate_eval.dscreened_rates_dT(k_He4_Ca48_to_n_Ti51) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 21.0_rt, 49.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Sc49_to_n_V52);
    rate_eval.screened_rates(k_He4_Sc49_to_n_V52) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Sc49_to_n_V52);
        rate_eval.dscreened_rates_dT(k_He4_Sc49_to_n_V52) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 22.0_rt, 51.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ti51_to_n_Cr54);
    rate_eval.screened_rates(k_He4_Ti51_to_n_Cr54) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ti51_to_n_Cr54);
        rate_eval.dscreened_rates_dT(k_He4_Ti51_to_n_Cr54) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 23.0_rt, 52.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_V52_to_n_Mn55);
    rate_eval.screened_rates(k_He4_V52_to_n_Mn55) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_V52_to_n_Mn55);
        rate_eval.dscreened_rates_dT(k_He4_V52_to_n_Mn55) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 28.0_rt, 63.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Ni63_to_n_Zn66);
    rate_eval.screened_rates(k_He4_Ni63_to_n_Zn66) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Ni63_to_n_Zn66);
        rate_eval.dscreened_rates_dT(k_He4_Ni63_to_n_Zn66) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 29.0_rt, 57.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Cu57_to_p_Zn60);
    rate_eval.screened_rates(k_He4_Cu57_to_p_Zn60) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cu57_to_p_Zn60);
        rate_eval.dscreened_rates_dT(k_He4_Cu57_to_p_Zn60) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 29.0_rt, 61.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Cu61_to_n_Ga64);
    rate_eval.screened_rates(k_He4_Cu61_to_n_Ga64) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cu61_to_n_Ga64);
        rate_eval.dscreened_rates_dT(k_He4_Cu61_to_n_Ga64) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Cu61_to_p_Zn64);
    rate_eval.screened_rates(k_He4_Cu61_to_p_Zn64) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cu61_to_p_Zn64);
        rate_eval.dscreened_rates_dT(k_He4_Cu61_to_p_Zn64) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 29.0_rt, 62.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Cu62_to_p_Zn65);
    rate_eval.screened_rates(k_He4_Cu62_to_p_Zn65) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cu62_to_p_Zn65);
        rate_eval.dscreened_rates_dT(k_He4_Cu62_to_p_Zn65) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 29.0_rt, 63.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Cu63_to_p_Zn66);
    rate_eval.screened_rates(k_He4_Cu63_to_p_Zn66) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Cu63_to_p_Zn66);
        rate_eval.dscreened_rates_dT(k_He4_Cu63_to_p_Zn66) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 30.0_rt, 60.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Zn60_to_He4_Cu57);
    rate_eval.screened_rates(k_p_Zn60_to_He4_Cu57) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Zn60_to_He4_Cu57);
        rate_eval.dscreened_rates_dT(k_p_Zn60_to_He4_Cu57) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 30.0_rt, 61.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_Zn61_to_n_Ge64);
    rate_eval.screened_rates(k_He4_Zn61_to_n_Ge64) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Zn61_to_n_Ge64);
        rate_eval.dscreened_rates_dT(k_He4_Zn61_to_n_Ge64) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_He4_Zn61_to_p_Ga64);
    rate_eval.screened_rates(k_He4_Zn61_to_p_Ga64) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_Zn61_to_p_Ga64);
        rate_eval.dscreened_rates_dT(k_He4_Zn61_to_p_Ga64) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 30.0_rt, 64.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Zn64_to_n_Ga64);
    rate_eval.screened_rates(k_p_Zn64_to_n_Ga64) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Zn64_to_n_Ga64);
        rate_eval.dscreened_rates_dT(k_p_Zn64_to_n_Ga64) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Zn64_to_He4_Cu61);
    rate_eval.screened_rates(k_p_Zn64_to_He4_Cu61) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Zn64_to_He4_Cu61);
        rate_eval.dscreened_rates_dT(k_p_Zn64_to_He4_Cu61) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 30.0_rt, 65.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Zn65_to_He4_Cu62);
    rate_eval.screened_rates(k_p_Zn65_to_He4_Cu62) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Zn65_to_He4_Cu62);
        rate_eval.dscreened_rates_dT(k_p_Zn65_to_He4_Cu62) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 30.0_rt, 66.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Zn66_to_He4_Cu63);
    rate_eval.screened_rates(k_p_Zn66_to_He4_Cu63) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Zn66_to_He4_Cu63);
        rate_eval.dscreened_rates_dT(k_p_Zn66_to_He4_Cu63) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 31.0_rt, 64.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_Ga64_to_n_Ge64);
    rate_eval.screened_rates(k_p_Ga64_to_n_Ge64) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ga64_to_n_Ge64);
        rate_eval.dscreened_rates_dT(k_p_Ga64_to_n_Ge64) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_Ga64_to_He4_Zn61);
    rate_eval.screened_rates(k_p_Ga64_to_He4_Zn61) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_Ga64_to_He4_Zn61);
        rate_eval.dscreened_rates_dT(k_p_Ga64_to_He4_Zn61) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 3.0_rt, 2.0_rt, 3.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He3_He3_to_p_p_He4);
    rate_eval.screened_rates(k_He3_He3_to_p_p_He4) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He3_He3_to_p_p_He4);
        rate_eval.dscreened_rates_dT(k_He3_He3_to_p_p_He4) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 2.0_rt, 3.0_rt, 7.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_d_Li7_to_n_He4_He4);
    rate_eval.screened_rates(k_d_Li7_to_n_He4_He4) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_d_Li7_to_n_He4_He4);
        rate_eval.dscreened_rates_dT(k_d_Li7_to_n_He4_He4) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 2.0_rt, 4.0_rt, 7.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_d_Be7_to_p_He4_He4);
    rate_eval.screened_rates(k_d_Be7_to_p_He4_He4) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_d_Be7_to_p_He4_He4);
        rate_eval.dscreened_rates_dT(k_d_Be7_to_p_He4_He4) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 3.0_rt, 3.0_rt, 7.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He3_Li7_to_n_p_He4_He4);
    rate_eval.screened_rates(k_He3_Li7_to_n_p_He4_He4) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He3_Li7_to_n_p_He4_He4);
        rate_eval.dscreened_rates_dT(k_He3_Li7_to_n_p_He4_He4) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 3.0_rt, 4.0_rt, 7.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He3_Be7_to_p_p_He4_He4);
    rate_eval.screened_rates(k_He3_Be7_to_p_p_He4_He4) *= scor;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He3_Be7_to_p_p_He4_He4);
        rate_eval.dscreened_rates_dT(k_He3_Be7_to_p_p_He4_He4) = ratraw * dscor_dt + dratraw_dT * scor;
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

