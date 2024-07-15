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

    using dual_t = std::conditional_t<do_T_derivatives, autodiff::dual, amrex::Real>;
    dual_t temp = state.T;
    if constexpr (do_T_derivatives) {
        // seed the dual number for temperature before calculating anything with it
        autodiff::seed(temp);
    }
    plasma_state_t<dual_t> pstate{};
    fill_plasma_state(pstate, temp, state.rho, Y);

    tf_t tfactors = evaluate_tfactors(state.T);

    fill_reaclib_rates<do_T_derivatives, T>(tfactors, rate_eval);



    // Evaluate screening factors

    amrex::Real ratraw, dratraw_dT;
    amrex::Real scor, dscor_dt;
    [[maybe_unused]] amrex::Real scor2, dscor2_dt;


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


    // Fill approximate rates

    fill_approx_rates<do_T_derivatives, T>(tfactors, rate_eval);

    // Calculate tabular rates

    [[maybe_unused]] amrex::Real rate, drate_dt, edot_nu, edot_gamma;

    rate_eval.enuc_weak = 0.0_rt;

    tabular_evaluate(j_F20_O20_meta, j_F20_O20_rhoy, j_F20_O20_temp, j_F20_O20_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_F20_to_O20) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_F20_to_O20) = drate_dt;
    }
    rate_eval.enuc_weak += C::Legacy::n_A * Y(F20) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Ne20_F20_meta, j_Ne20_F20_rhoy, j_Ne20_F20_temp, j_Ne20_F20_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Ne20_to_F20) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ne20_to_F20) = drate_dt;
    }
    rate_eval.enuc_weak += C::Legacy::n_A * Y(Ne20) * (edot_nu + edot_gamma);

    tabular_evaluate(j_O20_F20_meta, j_O20_F20_rhoy, j_O20_F20_temp, j_O20_F20_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_O20_to_F20) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_O20_to_F20) = drate_dt;
    }
    rate_eval.enuc_weak += C::Legacy::n_A * Y(O20) * (edot_nu + edot_gamma);

    tabular_evaluate(j_F20_Ne20_meta, j_F20_Ne20_rhoy, j_F20_Ne20_temp, j_F20_Ne20_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_F20_to_Ne20) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_F20_to_Ne20) = drate_dt;
    }
    rate_eval.enuc_weak += C::Legacy::n_A * Y(F20) * (edot_nu + edot_gamma);


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

    tabular_evaluate(j_F20_O20_meta, j_F20_O20_rhoy, j_F20_O20_temp, j_F20_O20_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_F20_to_O20) = rate;
    rate_eval.enuc_weak += C::Legacy::n_A * Y(F20) * (edot_nu + edot_gamma);

    tabular_evaluate(j_Ne20_F20_meta, j_Ne20_F20_rhoy, j_Ne20_F20_temp, j_Ne20_F20_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_Ne20_to_F20) = rate;
    rate_eval.enuc_weak += C::Legacy::n_A * Y(Ne20) * (edot_nu + edot_gamma);

    tabular_evaluate(j_O20_F20_meta, j_O20_F20_rhoy, j_O20_F20_temp, j_O20_F20_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_O20_to_F20) = rate;
    rate_eval.enuc_weak += C::Legacy::n_A * Y(O20) * (edot_nu + edot_gamma);

    tabular_evaluate(j_F20_Ne20_meta, j_F20_Ne20_rhoy, j_F20_Ne20_temp, j_F20_Ne20_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_F20_to_Ne20) = rate;
    rate_eval.enuc_weak += C::Legacy::n_A * Y(F20) * (edot_nu + edot_gamma);

    auto screened_rates = rate_eval.screened_rates;

    ydot_nuc(H1) = 0.0_rt;

    ydot_nuc(He4) = 0.0_rt;

    ydot_nuc(O16) = 0.0_rt;

    ydot_nuc(O20) =
        (-screened_rates(k_O20_to_F20)*Y(O20) + screened_rates(k_F20_to_O20)*Y(F20));

    ydot_nuc(F20) =
        (screened_rates(k_O20_to_F20)*Y(O20) + -screened_rates(k_F20_to_O20)*Y(F20)) +
        (-screened_rates(k_F20_to_Ne20)*Y(F20) + screened_rates(k_Ne20_to_F20)*Y(Ne20));

    ydot_nuc(Ne20) =
        (screened_rates(k_F20_to_Ne20)*Y(F20) + -screened_rates(k_Ne20_to_F20)*Y(Ne20));

    ydot_nuc(Mg24) = 0.0_rt;

    ydot_nuc(Al27) = 0.0_rt;

    ydot_nuc(Si28) = 0.0_rt;

    ydot_nuc(P31) = 0.0_rt;

    ydot_nuc(S32) = 0.0_rt;

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
        -screened_rates(k_p_Al27_to_Si28)*Y(Al27)*Y(H1)*state.rho +
        -screened_rates(k_p_P31_to_S32)*Y(P31)*Y(H1)*state.rho +
        0.5*screened_rates(k_O16_O16_to_p_P31)*amrex::Math::powi<2>(Y(O16))*state.rho +
        (-screened_rates(k_p_Al27_to_He4_Mg24)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_He4_Mg24_to_p_Al27)*Y(He4)*Y(Mg24)*state.rho) +
        (-screened_rates(k_p_P31_to_He4_Si28)*Y(P31)*Y(H1)*state.rho + screened_rates(k_He4_Si28_to_p_P31)*Y(He4)*Y(Si28)*state.rho);

    ydot_nuc(He4) =
        (-screened_rates(k_He4_O16_to_Ne20)*Y(He4)*Y(O16)*state.rho + screened_rates(k_Ne20_to_He4_O16)*Y(Ne20)) +
        -screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*Y(Ne20)*state.rho +
        -screened_rates(k_He4_Mg24_to_Si28)*Y(He4)*Y(Mg24)*state.rho +
        -screened_rates(k_He4_Al27_to_P31)*Y(Al27)*Y(He4)*state.rho +
        -screened_rates(k_He4_Si28_to_S32)*Y(He4)*Y(Si28)*state.rho +
        0.5*screened_rates(k_O16_O16_to_He4_Si28)*amrex::Math::powi<2>(Y(O16))*state.rho +
        (screened_rates(k_p_Al27_to_He4_Mg24)*Y(Al27)*Y(H1)*state.rho + -screened_rates(k_He4_Mg24_to_p_Al27)*Y(He4)*Y(Mg24)*state.rho) +
        (screened_rates(k_p_P31_to_He4_Si28)*Y(P31)*Y(H1)*state.rho + -screened_rates(k_He4_Si28_to_p_P31)*Y(He4)*Y(Si28)*state.rho);

    ydot_nuc(O16) =
        (-screened_rates(k_He4_O16_to_Ne20)*Y(He4)*Y(O16)*state.rho + screened_rates(k_Ne20_to_He4_O16)*Y(Ne20)) +
        -screened_rates(k_O16_O16_to_p_P31)*amrex::Math::powi<2>(Y(O16))*state.rho +
        -screened_rates(k_O16_O16_to_He4_Si28)*amrex::Math::powi<2>(Y(O16))*state.rho;

    ydot_nuc(O20) =
        (-screened_rates(k_O20_to_F20)*Y(O20) + screened_rates(k_F20_to_O20)*Y(F20));

    ydot_nuc(F20) =
        (screened_rates(k_O20_to_F20)*Y(O20) + -screened_rates(k_F20_to_O20)*Y(F20)) +
        (-screened_rates(k_F20_to_Ne20)*Y(F20) + screened_rates(k_Ne20_to_F20)*Y(Ne20));

    ydot_nuc(Ne20) =
        (screened_rates(k_F20_to_Ne20)*Y(F20) + -screened_rates(k_Ne20_to_F20)*Y(Ne20)) +
        (screened_rates(k_He4_O16_to_Ne20)*Y(He4)*Y(O16)*state.rho + -screened_rates(k_Ne20_to_He4_O16)*Y(Ne20)) +
        -screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*Y(Ne20)*state.rho;

    ydot_nuc(Mg24) =
        screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*Y(Ne20)*state.rho +
        -screened_rates(k_He4_Mg24_to_Si28)*Y(He4)*Y(Mg24)*state.rho +
        (screened_rates(k_p_Al27_to_He4_Mg24)*Y(Al27)*Y(H1)*state.rho + -screened_rates(k_He4_Mg24_to_p_Al27)*Y(He4)*Y(Mg24)*state.rho);

    ydot_nuc(Al27) =
        -screened_rates(k_p_Al27_to_Si28)*Y(Al27)*Y(H1)*state.rho +
        -screened_rates(k_He4_Al27_to_P31)*Y(Al27)*Y(He4)*state.rho +
        (-screened_rates(k_p_Al27_to_He4_Mg24)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_He4_Mg24_to_p_Al27)*Y(He4)*Y(Mg24)*state.rho);

    ydot_nuc(Si28) =
        screened_rates(k_He4_Mg24_to_Si28)*Y(He4)*Y(Mg24)*state.rho +
        screened_rates(k_p_Al27_to_Si28)*Y(Al27)*Y(H1)*state.rho +
        -screened_rates(k_He4_Si28_to_S32)*Y(He4)*Y(Si28)*state.rho +
        0.5*screened_rates(k_O16_O16_to_He4_Si28)*amrex::Math::powi<2>(Y(O16))*state.rho +
        (screened_rates(k_p_P31_to_He4_Si28)*Y(P31)*Y(H1)*state.rho + -screened_rates(k_He4_Si28_to_p_P31)*Y(He4)*Y(Si28)*state.rho);

    ydot_nuc(P31) =
        screened_rates(k_He4_Al27_to_P31)*Y(Al27)*Y(He4)*state.rho +
        -screened_rates(k_p_P31_to_S32)*Y(P31)*Y(H1)*state.rho +
        0.5*screened_rates(k_O16_O16_to_p_P31)*amrex::Math::powi<2>(Y(O16))*state.rho +
        (-screened_rates(k_p_P31_to_He4_Si28)*Y(P31)*Y(H1)*state.rho + screened_rates(k_He4_Si28_to_p_P31)*Y(He4)*Y(Si28)*state.rho);

    ydot_nuc(S32) =
        screened_rates(k_He4_Si28_to_S32)*Y(He4)*Y(Si28)*state.rho +
        screened_rates(k_p_P31_to_S32)*Y(P31)*Y(H1)*state.rho;

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

    scratch = -screened_rates(k_p_Al27_to_He4_Mg24)*Y(Al27)*state.rho - screened_rates(k_p_Al27_to_Si28)*Y(Al27)*state.rho - screened_rates(k_p_P31_to_He4_Si28)*Y(P31)*state.rho - screened_rates(k_p_P31_to_S32)*Y(P31)*state.rho;
    jac.set(H1, H1, scratch);

    scratch = screened_rates(k_He4_Mg24_to_p_Al27)*Y(Mg24)*state.rho + screened_rates(k_He4_Si28_to_p_P31)*Y(Si28)*state.rho;
    jac.set(H1, He4, scratch);

    scratch = 1.0*screened_rates(k_O16_O16_to_p_P31)*Y(O16)*state.rho;
    jac.set(H1, O16, scratch);

    scratch = screened_rates(k_He4_Mg24_to_p_Al27)*Y(He4)*state.rho;
    jac.set(H1, Mg24, scratch);

    scratch = -screened_rates(k_p_Al27_to_He4_Mg24)*Y(H1)*state.rho - screened_rates(k_p_Al27_to_Si28)*Y(H1)*state.rho;
    jac.set(H1, Al27, scratch);

    scratch = screened_rates(k_He4_Si28_to_p_P31)*Y(He4)*state.rho;
    jac.set(H1, Si28, scratch);

    scratch = -screened_rates(k_p_P31_to_He4_Si28)*Y(H1)*state.rho - screened_rates(k_p_P31_to_S32)*Y(H1)*state.rho;
    jac.set(H1, P31, scratch);

    scratch = screened_rates(k_p_Al27_to_He4_Mg24)*Y(Al27)*state.rho + screened_rates(k_p_P31_to_He4_Si28)*Y(P31)*state.rho;
    jac.set(He4, H1, scratch);

    scratch = -screened_rates(k_He4_Al27_to_P31)*Y(Al27)*state.rho - screened_rates(k_He4_Mg24_to_Si28)*Y(Mg24)*state.rho - screened_rates(k_He4_Mg24_to_p_Al27)*Y(Mg24)*state.rho - screened_rates(k_He4_Ne20_to_Mg24)*Y(Ne20)*state.rho - screened_rates(k_He4_O16_to_Ne20)*Y(O16)*state.rho - screened_rates(k_He4_Si28_to_S32)*Y(Si28)*state.rho - screened_rates(k_He4_Si28_to_p_P31)*Y(Si28)*state.rho;
    jac.set(He4, He4, scratch);

    scratch = -screened_rates(k_He4_O16_to_Ne20)*Y(He4)*state.rho + 1.0*screened_rates(k_O16_O16_to_He4_Si28)*Y(O16)*state.rho;
    jac.set(He4, O16, scratch);

    scratch = -screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*state.rho + screened_rates(k_Ne20_to_He4_O16);
    jac.set(He4, Ne20, scratch);

    scratch = -screened_rates(k_He4_Mg24_to_Si28)*Y(He4)*state.rho - screened_rates(k_He4_Mg24_to_p_Al27)*Y(He4)*state.rho;
    jac.set(He4, Mg24, scratch);

    scratch = -screened_rates(k_He4_Al27_to_P31)*Y(He4)*state.rho + screened_rates(k_p_Al27_to_He4_Mg24)*Y(H1)*state.rho;
    jac.set(He4, Al27, scratch);

    scratch = -screened_rates(k_He4_Si28_to_S32)*Y(He4)*state.rho - screened_rates(k_He4_Si28_to_p_P31)*Y(He4)*state.rho;
    jac.set(He4, Si28, scratch);

    scratch = screened_rates(k_p_P31_to_He4_Si28)*Y(H1)*state.rho;
    jac.set(He4, P31, scratch);

    scratch = -screened_rates(k_He4_O16_to_Ne20)*Y(O16)*state.rho;
    jac.set(O16, He4, scratch);

    scratch = -screened_rates(k_He4_O16_to_Ne20)*Y(He4)*state.rho - 2.0*screened_rates(k_O16_O16_to_He4_Si28)*Y(O16)*state.rho - 2.0*screened_rates(k_O16_O16_to_p_P31)*Y(O16)*state.rho;
    jac.set(O16, O16, scratch);

    scratch = screened_rates(k_Ne20_to_He4_O16);
    jac.set(O16, Ne20, scratch);

    scratch = -screened_rates(k_O20_to_F20);
    jac.set(O20, O20, scratch);

    scratch = screened_rates(k_F20_to_O20);
    jac.set(O20, F20, scratch);

    scratch = screened_rates(k_O20_to_F20);
    jac.set(F20, O20, scratch);

    scratch = -screened_rates(k_F20_to_Ne20) - screened_rates(k_F20_to_O20);
    jac.set(F20, F20, scratch);

    scratch = screened_rates(k_Ne20_to_F20);
    jac.set(F20, Ne20, scratch);

    scratch = -screened_rates(k_He4_Ne20_to_Mg24)*Y(Ne20)*state.rho + screened_rates(k_He4_O16_to_Ne20)*Y(O16)*state.rho;
    jac.set(Ne20, He4, scratch);

    scratch = screened_rates(k_He4_O16_to_Ne20)*Y(He4)*state.rho;
    jac.set(Ne20, O16, scratch);

    scratch = screened_rates(k_F20_to_Ne20);
    jac.set(Ne20, F20, scratch);

    scratch = -screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*state.rho - screened_rates(k_Ne20_to_F20) - screened_rates(k_Ne20_to_He4_O16);
    jac.set(Ne20, Ne20, scratch);

    scratch = screened_rates(k_p_Al27_to_He4_Mg24)*Y(Al27)*state.rho;
    jac.set(Mg24, H1, scratch);

    scratch = -screened_rates(k_He4_Mg24_to_Si28)*Y(Mg24)*state.rho - screened_rates(k_He4_Mg24_to_p_Al27)*Y(Mg24)*state.rho + screened_rates(k_He4_Ne20_to_Mg24)*Y(Ne20)*state.rho;
    jac.set(Mg24, He4, scratch);

    scratch = screened_rates(k_He4_Ne20_to_Mg24)*Y(He4)*state.rho;
    jac.set(Mg24, Ne20, scratch);

    scratch = -screened_rates(k_He4_Mg24_to_Si28)*Y(He4)*state.rho - screened_rates(k_He4_Mg24_to_p_Al27)*Y(He4)*state.rho;
    jac.set(Mg24, Mg24, scratch);

    scratch = screened_rates(k_p_Al27_to_He4_Mg24)*Y(H1)*state.rho;
    jac.set(Mg24, Al27, scratch);

    scratch = -screened_rates(k_p_Al27_to_He4_Mg24)*Y(Al27)*state.rho - screened_rates(k_p_Al27_to_Si28)*Y(Al27)*state.rho;
    jac.set(Al27, H1, scratch);

    scratch = -screened_rates(k_He4_Al27_to_P31)*Y(Al27)*state.rho + screened_rates(k_He4_Mg24_to_p_Al27)*Y(Mg24)*state.rho;
    jac.set(Al27, He4, scratch);

    scratch = screened_rates(k_He4_Mg24_to_p_Al27)*Y(He4)*state.rho;
    jac.set(Al27, Mg24, scratch);

    scratch = -screened_rates(k_He4_Al27_to_P31)*Y(He4)*state.rho - screened_rates(k_p_Al27_to_He4_Mg24)*Y(H1)*state.rho - screened_rates(k_p_Al27_to_Si28)*Y(H1)*state.rho;
    jac.set(Al27, Al27, scratch);

    scratch = screened_rates(k_p_Al27_to_Si28)*Y(Al27)*state.rho + screened_rates(k_p_P31_to_He4_Si28)*Y(P31)*state.rho;
    jac.set(Si28, H1, scratch);

    scratch = screened_rates(k_He4_Mg24_to_Si28)*Y(Mg24)*state.rho - screened_rates(k_He4_Si28_to_S32)*Y(Si28)*state.rho - screened_rates(k_He4_Si28_to_p_P31)*Y(Si28)*state.rho;
    jac.set(Si28, He4, scratch);

    scratch = 1.0*screened_rates(k_O16_O16_to_He4_Si28)*Y(O16)*state.rho;
    jac.set(Si28, O16, scratch);

    scratch = screened_rates(k_He4_Mg24_to_Si28)*Y(He4)*state.rho;
    jac.set(Si28, Mg24, scratch);

    scratch = screened_rates(k_p_Al27_to_Si28)*Y(H1)*state.rho;
    jac.set(Si28, Al27, scratch);

    scratch = -screened_rates(k_He4_Si28_to_S32)*Y(He4)*state.rho - screened_rates(k_He4_Si28_to_p_P31)*Y(He4)*state.rho;
    jac.set(Si28, Si28, scratch);

    scratch = screened_rates(k_p_P31_to_He4_Si28)*Y(H1)*state.rho;
    jac.set(Si28, P31, scratch);

    scratch = -screened_rates(k_p_P31_to_He4_Si28)*Y(P31)*state.rho - screened_rates(k_p_P31_to_S32)*Y(P31)*state.rho;
    jac.set(P31, H1, scratch);

    scratch = screened_rates(k_He4_Al27_to_P31)*Y(Al27)*state.rho + screened_rates(k_He4_Si28_to_p_P31)*Y(Si28)*state.rho;
    jac.set(P31, He4, scratch);

    scratch = 1.0*screened_rates(k_O16_O16_to_p_P31)*Y(O16)*state.rho;
    jac.set(P31, O16, scratch);

    scratch = screened_rates(k_He4_Al27_to_P31)*Y(He4)*state.rho;
    jac.set(P31, Al27, scratch);

    scratch = screened_rates(k_He4_Si28_to_p_P31)*Y(He4)*state.rho;
    jac.set(P31, Si28, scratch);

    scratch = -screened_rates(k_p_P31_to_He4_Si28)*Y(H1)*state.rho - screened_rates(k_p_P31_to_S32)*Y(H1)*state.rho;
    jac.set(P31, P31, scratch);

    scratch = screened_rates(k_p_P31_to_S32)*Y(P31)*state.rho;
    jac.set(S32, H1, scratch);

    scratch = screened_rates(k_He4_Si28_to_S32)*Y(Si28)*state.rho;
    jac.set(S32, He4, scratch);

    scratch = screened_rates(k_He4_Si28_to_S32)*Y(He4)*state.rho;
    jac.set(S32, Si28, scratch);

    scratch = screened_rates(k_p_P31_to_S32)*Y(H1)*state.rho;
    jac.set(S32, P31, scratch);


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
