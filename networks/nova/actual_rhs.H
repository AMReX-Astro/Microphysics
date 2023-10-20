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


    ratraw = rate_eval.screened_rates(k_He4_C12_to_O16);
    rate_eval.screened_rates(k_He4_C12_to_O16) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_C12_to_O16);
        rate_eval.dscreened_rates_dT(k_He4_C12_to_O16) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 6.0_rt, 13.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_C13_to_N14);
    rate_eval.screened_rates(k_p_C13_to_N14) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
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
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
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
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
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
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_N14_to_F18);
        rate_eval.dscreened_rates_dT(k_He4_N14_to_F18) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 7.0_rt, 15.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_N15_to_O16);
    rate_eval.screened_rates(k_p_N15_to_O16) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_N15_to_O16);
        rate_eval.dscreened_rates_dT(k_p_N15_to_O16) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_N15_to_He4_C12);
    rate_eval.screened_rates(k_p_N15_to_He4_C12) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_N15_to_He4_C12);
        rate_eval.dscreened_rates_dT(k_p_N15_to_He4_C12) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 8.0_rt, 16.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_O16_to_F17);
    rate_eval.screened_rates(k_p_O16_to_F17) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_O16_to_F17);
        rate_eval.dscreened_rates_dT(k_p_O16_to_F17) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 8.0_rt, 17.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_O17_to_F18);
    rate_eval.screened_rates(k_p_O17_to_F18) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_O17_to_F18);
        rate_eval.dscreened_rates_dT(k_p_O17_to_F18) = ratraw * dscor_dt + dratraw_dT * scor;
    }

    ratraw = rate_eval.screened_rates(k_p_O17_to_He4_N14);
    rate_eval.screened_rates(k_p_O17_to_He4_N14) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_O17_to_He4_N14);
        rate_eval.dscreened_rates_dT(k_p_O17_to_He4_N14) = ratraw * dscor_dt + dratraw_dT * scor;
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


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 8.0_rt, 14.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_He4_O14_to_p_F17);
    rate_eval.screened_rates(k_He4_O14_to_p_F17) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_O14_to_p_F17);
        rate_eval.dscreened_rates_dT(k_He4_O14_to_p_F17) = ratraw * dscor_dt + dratraw_dT * scor;
    }


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 9.0_rt, 18.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_F18_to_He4_O15);
    rate_eval.screened_rates(k_p_F18_to_He4_O15) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_p_F18_to_He4_O15);
        rate_eval.dscreened_rates_dT(k_p_F18_to_He4_O15) = ratraw * dscor_dt + dratraw_dT * scor;
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
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_He4_He4_He4_to_C12);
        rate_eval.dscreened_rates_dT(k_He4_He4_He4_to_C12) = ratraw * (scor * dscor2_dt + dscor_dt * scor2) + dratraw_dT * scor * scor2;
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

    ydot_nuc(H1) =
        -screened_rates(k_p_C12_to_N13)*Y(C12)*Y(H1)*state.rho +
        -screened_rates(k_p_C13_to_N14)*Y(C13)*Y(H1)*state.rho +
        -screened_rates(k_p_N13_to_O14)*Y(N13)*Y(H1)*state.rho +
        -screened_rates(k_p_N14_to_O15)*Y(N14)*Y(H1)*state.rho +
        -screened_rates(k_p_N15_to_O16)*Y(N15)*Y(H1)*state.rho +
        -screened_rates(k_p_O16_to_F17)*Y(O16)*Y(H1)*state.rho +
        -screened_rates(k_p_O17_to_F18)*Y(O17)*Y(H1)*state.rho +
        screened_rates(k_He4_N13_to_p_O16)*Y(He4)*Y(N13)*state.rho +
        -screened_rates(k_p_N15_to_He4_C12)*Y(N15)*Y(H1)*state.rho +
        screened_rates(k_He4_O14_to_p_F17)*Y(He4)*Y(O14)*state.rho +
        -screened_rates(k_p_O17_to_He4_N14)*Y(O17)*Y(H1)*state.rho +
        -screened_rates(k_p_F18_to_He4_O15)*Y(F18)*Y(H1)*state.rho;

    ydot_nuc(He4) =
        -screened_rates(k_He4_C12_to_O16)*Y(C12)*Y(He4)*state.rho +
        -screened_rates(k_He4_N14_to_F18)*Y(He4)*Y(N14)*state.rho +
        -screened_rates(k_He4_N13_to_p_O16)*Y(He4)*Y(N13)*state.rho +
        screened_rates(k_p_N15_to_He4_C12)*Y(N15)*Y(H1)*state.rho +
        -screened_rates(k_He4_O14_to_p_F17)*Y(He4)*Y(O14)*state.rho +
        screened_rates(k_p_O17_to_He4_N14)*Y(O17)*Y(H1)*state.rho +
        screened_rates(k_p_F18_to_He4_O15)*Y(F18)*Y(H1)*state.rho +
        -0.5*screened_rates(k_He4_He4_He4_to_C12)*std::pow(Y(He4), 3)*std::pow(state.rho, 2);

    ydot_nuc(C12) =
        -screened_rates(k_p_C12_to_N13)*Y(C12)*Y(H1)*state.rho +
        -screened_rates(k_He4_C12_to_O16)*Y(C12)*Y(He4)*state.rho +
        screened_rates(k_p_N15_to_He4_C12)*Y(N15)*Y(H1)*state.rho +
        0.16666666666666667*screened_rates(k_He4_He4_He4_to_C12)*std::pow(Y(He4), 3)*std::pow(state.rho, 2);

    ydot_nuc(C13) =
        screened_rates(k_N13_to_C13_weak_wc12)*Y(N13) +
        -screened_rates(k_p_C13_to_N14)*Y(C13)*Y(H1)*state.rho;

    ydot_nuc(N13) =
        -screened_rates(k_N13_to_C13_weak_wc12)*Y(N13) +
        screened_rates(k_p_C12_to_N13)*Y(C12)*Y(H1)*state.rho +
        -screened_rates(k_p_N13_to_O14)*Y(N13)*Y(H1)*state.rho +
        -screened_rates(k_He4_N13_to_p_O16)*Y(He4)*Y(N13)*state.rho;

    ydot_nuc(N14) =
        screened_rates(k_O14_to_N14_weak_wc12)*Y(O14) +
        screened_rates(k_p_C13_to_N14)*Y(C13)*Y(H1)*state.rho +
        -screened_rates(k_p_N14_to_O15)*Y(N14)*Y(H1)*state.rho +
        -screened_rates(k_He4_N14_to_F18)*Y(He4)*Y(N14)*state.rho +
        screened_rates(k_p_O17_to_He4_N14)*Y(O17)*Y(H1)*state.rho;

    ydot_nuc(N15) =
        screened_rates(k_O15_to_N15_weak_wc12)*Y(O15) +
        -screened_rates(k_p_N15_to_O16)*Y(N15)*Y(H1)*state.rho +
        -screened_rates(k_p_N15_to_He4_C12)*Y(N15)*Y(H1)*state.rho;

    ydot_nuc(O14) =
        -screened_rates(k_O14_to_N14_weak_wc12)*Y(O14) +
        screened_rates(k_p_N13_to_O14)*Y(N13)*Y(H1)*state.rho +
        -screened_rates(k_He4_O14_to_p_F17)*Y(He4)*Y(O14)*state.rho;

    ydot_nuc(O15) =
        -screened_rates(k_O15_to_N15_weak_wc12)*Y(O15) +
        screened_rates(k_p_N14_to_O15)*Y(N14)*Y(H1)*state.rho +
        screened_rates(k_p_F18_to_He4_O15)*Y(F18)*Y(H1)*state.rho;

    ydot_nuc(O16) =
        screened_rates(k_He4_C12_to_O16)*Y(C12)*Y(He4)*state.rho +
        screened_rates(k_p_N15_to_O16)*Y(N15)*Y(H1)*state.rho +
        -screened_rates(k_p_O16_to_F17)*Y(O16)*Y(H1)*state.rho +
        screened_rates(k_He4_N13_to_p_O16)*Y(He4)*Y(N13)*state.rho;

    ydot_nuc(O17) =
        screened_rates(k_F17_to_O17_weak_wc12)*Y(F17) +
        -screened_rates(k_p_O17_to_F18)*Y(O17)*Y(H1)*state.rho +
        -screened_rates(k_p_O17_to_He4_N14)*Y(O17)*Y(H1)*state.rho;

    ydot_nuc(F17) =
        -screened_rates(k_F17_to_O17_weak_wc12)*Y(F17) +
        screened_rates(k_p_O16_to_F17)*Y(O16)*Y(H1)*state.rho +
        screened_rates(k_He4_O14_to_p_F17)*Y(He4)*Y(O14)*state.rho;

    ydot_nuc(F18) =
        screened_rates(k_He4_N14_to_F18)*Y(He4)*Y(N14)*state.rho +
        screened_rates(k_p_O17_to_F18)*Y(O17)*Y(H1)*state.rho +
        -screened_rates(k_p_F18_to_He4_O15)*Y(F18)*Y(H1)*state.rho;

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

    scratch = -screened_rates(k_p_C12_to_N13)*Y(C12)*state.rho - screened_rates(k_p_C13_to_N14)*Y(C13)*state.rho - screened_rates(k_p_F18_to_He4_O15)*Y(F18)*state.rho - screened_rates(k_p_N13_to_O14)*Y(N13)*state.rho - screened_rates(k_p_N14_to_O15)*Y(N14)*state.rho - screened_rates(k_p_N15_to_He4_C12)*Y(N15)*state.rho - screened_rates(k_p_N15_to_O16)*Y(N15)*state.rho - screened_rates(k_p_O16_to_F17)*Y(O16)*state.rho - screened_rates(k_p_O17_to_F18)*Y(O17)*state.rho - screened_rates(k_p_O17_to_He4_N14)*Y(O17)*state.rho;
    jac.set(H1, H1, scratch);

    scratch = screened_rates(k_He4_N13_to_p_O16)*Y(N13)*state.rho + screened_rates(k_He4_O14_to_p_F17)*Y(O14)*state.rho;
    jac.set(H1, He4, scratch);

    scratch = -screened_rates(k_p_C12_to_N13)*Y(H1)*state.rho;
    jac.set(H1, C12, scratch);

    scratch = -screened_rates(k_p_C13_to_N14)*Y(H1)*state.rho;
    jac.set(H1, C13, scratch);

    scratch = screened_rates(k_He4_N13_to_p_O16)*Y(He4)*state.rho - screened_rates(k_p_N13_to_O14)*Y(H1)*state.rho;
    jac.set(H1, N13, scratch);

    scratch = -screened_rates(k_p_N14_to_O15)*Y(H1)*state.rho;
    jac.set(H1, N14, scratch);

    scratch = -screened_rates(k_p_N15_to_He4_C12)*Y(H1)*state.rho - screened_rates(k_p_N15_to_O16)*Y(H1)*state.rho;
    jac.set(H1, N15, scratch);

    scratch = screened_rates(k_He4_O14_to_p_F17)*Y(He4)*state.rho;
    jac.set(H1, O14, scratch);

    scratch = -screened_rates(k_p_O16_to_F17)*Y(H1)*state.rho;
    jac.set(H1, O16, scratch);

    scratch = -screened_rates(k_p_O17_to_F18)*Y(H1)*state.rho - screened_rates(k_p_O17_to_He4_N14)*Y(H1)*state.rho;
    jac.set(H1, O17, scratch);

    scratch = -screened_rates(k_p_F18_to_He4_O15)*Y(H1)*state.rho;
    jac.set(H1, F18, scratch);

    scratch = screened_rates(k_p_F18_to_He4_O15)*Y(F18)*state.rho + screened_rates(k_p_N15_to_He4_C12)*Y(N15)*state.rho + screened_rates(k_p_O17_to_He4_N14)*Y(O17)*state.rho;
    jac.set(He4, H1, scratch);

    scratch = -screened_rates(k_He4_C12_to_O16)*Y(C12)*state.rho - 1.5*screened_rates(k_He4_He4_He4_to_C12)*std::pow(Y(He4), 2)*std::pow(state.rho, 2) - screened_rates(k_He4_N13_to_p_O16)*Y(N13)*state.rho - screened_rates(k_He4_N14_to_F18)*Y(N14)*state.rho - screened_rates(k_He4_O14_to_p_F17)*Y(O14)*state.rho;
    jac.set(He4, He4, scratch);

    scratch = -screened_rates(k_He4_C12_to_O16)*Y(He4)*state.rho;
    jac.set(He4, C12, scratch);

    scratch = -screened_rates(k_He4_N13_to_p_O16)*Y(He4)*state.rho;
    jac.set(He4, N13, scratch);

    scratch = -screened_rates(k_He4_N14_to_F18)*Y(He4)*state.rho;
    jac.set(He4, N14, scratch);

    scratch = screened_rates(k_p_N15_to_He4_C12)*Y(H1)*state.rho;
    jac.set(He4, N15, scratch);

    scratch = -screened_rates(k_He4_O14_to_p_F17)*Y(He4)*state.rho;
    jac.set(He4, O14, scratch);

    scratch = screened_rates(k_p_O17_to_He4_N14)*Y(H1)*state.rho;
    jac.set(He4, O17, scratch);

    scratch = screened_rates(k_p_F18_to_He4_O15)*Y(H1)*state.rho;
    jac.set(He4, F18, scratch);

    scratch = -screened_rates(k_p_C12_to_N13)*Y(C12)*state.rho + screened_rates(k_p_N15_to_He4_C12)*Y(N15)*state.rho;
    jac.set(C12, H1, scratch);

    scratch = -screened_rates(k_He4_C12_to_O16)*Y(C12)*state.rho + 0.5*screened_rates(k_He4_He4_He4_to_C12)*std::pow(Y(He4), 2)*std::pow(state.rho, 2);
    jac.set(C12, He4, scratch);

    scratch = -screened_rates(k_He4_C12_to_O16)*Y(He4)*state.rho - screened_rates(k_p_C12_to_N13)*Y(H1)*state.rho;
    jac.set(C12, C12, scratch);

    scratch = screened_rates(k_p_N15_to_He4_C12)*Y(H1)*state.rho;
    jac.set(C12, N15, scratch);

    scratch = -screened_rates(k_p_C13_to_N14)*Y(C13)*state.rho;
    jac.set(C13, H1, scratch);

    scratch = -screened_rates(k_p_C13_to_N14)*Y(H1)*state.rho;
    jac.set(C13, C13, scratch);

    scratch = screened_rates(k_N13_to_C13_weak_wc12);
    jac.set(C13, N13, scratch);

    scratch = screened_rates(k_p_C12_to_N13)*Y(C12)*state.rho - screened_rates(k_p_N13_to_O14)*Y(N13)*state.rho;
    jac.set(N13, H1, scratch);

    scratch = -screened_rates(k_He4_N13_to_p_O16)*Y(N13)*state.rho;
    jac.set(N13, He4, scratch);

    scratch = screened_rates(k_p_C12_to_N13)*Y(H1)*state.rho;
    jac.set(N13, C12, scratch);

    scratch = -screened_rates(k_He4_N13_to_p_O16)*Y(He4)*state.rho - screened_rates(k_N13_to_C13_weak_wc12) - screened_rates(k_p_N13_to_O14)*Y(H1)*state.rho;
    jac.set(N13, N13, scratch);

    scratch = screened_rates(k_p_C13_to_N14)*Y(C13)*state.rho - screened_rates(k_p_N14_to_O15)*Y(N14)*state.rho + screened_rates(k_p_O17_to_He4_N14)*Y(O17)*state.rho;
    jac.set(N14, H1, scratch);

    scratch = -screened_rates(k_He4_N14_to_F18)*Y(N14)*state.rho;
    jac.set(N14, He4, scratch);

    scratch = screened_rates(k_p_C13_to_N14)*Y(H1)*state.rho;
    jac.set(N14, C13, scratch);

    scratch = -screened_rates(k_He4_N14_to_F18)*Y(He4)*state.rho - screened_rates(k_p_N14_to_O15)*Y(H1)*state.rho;
    jac.set(N14, N14, scratch);

    scratch = screened_rates(k_O14_to_N14_weak_wc12);
    jac.set(N14, O14, scratch);

    scratch = screened_rates(k_p_O17_to_He4_N14)*Y(H1)*state.rho;
    jac.set(N14, O17, scratch);

    scratch = -screened_rates(k_p_N15_to_He4_C12)*Y(N15)*state.rho - screened_rates(k_p_N15_to_O16)*Y(N15)*state.rho;
    jac.set(N15, H1, scratch);

    scratch = -screened_rates(k_p_N15_to_He4_C12)*Y(H1)*state.rho - screened_rates(k_p_N15_to_O16)*Y(H1)*state.rho;
    jac.set(N15, N15, scratch);

    scratch = screened_rates(k_O15_to_N15_weak_wc12);
    jac.set(N15, O15, scratch);

    scratch = screened_rates(k_p_N13_to_O14)*Y(N13)*state.rho;
    jac.set(O14, H1, scratch);

    scratch = -screened_rates(k_He4_O14_to_p_F17)*Y(O14)*state.rho;
    jac.set(O14, He4, scratch);

    scratch = screened_rates(k_p_N13_to_O14)*Y(H1)*state.rho;
    jac.set(O14, N13, scratch);

    scratch = -screened_rates(k_He4_O14_to_p_F17)*Y(He4)*state.rho - screened_rates(k_O14_to_N14_weak_wc12);
    jac.set(O14, O14, scratch);

    scratch = screened_rates(k_p_F18_to_He4_O15)*Y(F18)*state.rho + screened_rates(k_p_N14_to_O15)*Y(N14)*state.rho;
    jac.set(O15, H1, scratch);

    scratch = screened_rates(k_p_N14_to_O15)*Y(H1)*state.rho;
    jac.set(O15, N14, scratch);

    scratch = -screened_rates(k_O15_to_N15_weak_wc12);
    jac.set(O15, O15, scratch);

    scratch = screened_rates(k_p_F18_to_He4_O15)*Y(H1)*state.rho;
    jac.set(O15, F18, scratch);

    scratch = screened_rates(k_p_N15_to_O16)*Y(N15)*state.rho - screened_rates(k_p_O16_to_F17)*Y(O16)*state.rho;
    jac.set(O16, H1, scratch);

    scratch = screened_rates(k_He4_C12_to_O16)*Y(C12)*state.rho + screened_rates(k_He4_N13_to_p_O16)*Y(N13)*state.rho;
    jac.set(O16, He4, scratch);

    scratch = screened_rates(k_He4_C12_to_O16)*Y(He4)*state.rho;
    jac.set(O16, C12, scratch);

    scratch = screened_rates(k_He4_N13_to_p_O16)*Y(He4)*state.rho;
    jac.set(O16, N13, scratch);

    scratch = screened_rates(k_p_N15_to_O16)*Y(H1)*state.rho;
    jac.set(O16, N15, scratch);

    scratch = -screened_rates(k_p_O16_to_F17)*Y(H1)*state.rho;
    jac.set(O16, O16, scratch);

    scratch = -screened_rates(k_p_O17_to_F18)*Y(O17)*state.rho - screened_rates(k_p_O17_to_He4_N14)*Y(O17)*state.rho;
    jac.set(O17, H1, scratch);

    scratch = -screened_rates(k_p_O17_to_F18)*Y(H1)*state.rho - screened_rates(k_p_O17_to_He4_N14)*Y(H1)*state.rho;
    jac.set(O17, O17, scratch);

    scratch = screened_rates(k_F17_to_O17_weak_wc12);
    jac.set(O17, F17, scratch);

    scratch = screened_rates(k_p_O16_to_F17)*Y(O16)*state.rho;
    jac.set(F17, H1, scratch);

    scratch = screened_rates(k_He4_O14_to_p_F17)*Y(O14)*state.rho;
    jac.set(F17, He4, scratch);

    scratch = screened_rates(k_He4_O14_to_p_F17)*Y(He4)*state.rho;
    jac.set(F17, O14, scratch);

    scratch = screened_rates(k_p_O16_to_F17)*Y(H1)*state.rho;
    jac.set(F17, O16, scratch);

    scratch = -screened_rates(k_F17_to_O17_weak_wc12);
    jac.set(F17, F17, scratch);

    scratch = -screened_rates(k_p_F18_to_He4_O15)*Y(F18)*state.rho + screened_rates(k_p_O17_to_F18)*Y(O17)*state.rho;
    jac.set(F18, H1, scratch);

    scratch = screened_rates(k_He4_N14_to_F18)*Y(N14)*state.rho;
    jac.set(F18, He4, scratch);

    scratch = screened_rates(k_He4_N14_to_F18)*Y(He4)*state.rho;
    jac.set(F18, N14, scratch);

    scratch = screened_rates(k_p_O17_to_F18)*Y(H1)*state.rho;
    jac.set(F18, O17, scratch);

    scratch = -screened_rates(k_p_F18_to_He4_O15)*Y(H1)*state.rho;
    jac.set(F18, F18, scratch);


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
