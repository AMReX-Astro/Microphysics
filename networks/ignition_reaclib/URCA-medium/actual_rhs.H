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

    ratraw = rate_eval.screened_rates(k_he4_ne20_to_p_na23);
    rate_eval.screened_rates(k_he4_ne20_to_p_na23) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_ne20_to_p_na23);
        rate_eval.dscreened_rates_dT(k_he4_ne20_to_p_na23) = ratraw * dscor_dt + dratraw_dT * scor;
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
        constexpr auto scn_fac = scrn::calculate_screen_factor(6.0_rt, 12.0_rt, 6.0_rt, 12.0_rt);


        static_assert(scn_fac.z1 == 6.0_rt);


        actual_screen<do_T_derivatives>(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_c12_c12_to_n_mg23);
    rate_eval.screened_rates(k_c12_c12_to_n_mg23) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_c12_c12_to_n_mg23);
        rate_eval.dscreened_rates_dT(k_c12_c12_to_n_mg23) = ratraw * dscor_dt + dratraw_dT * scor;
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


    ratraw = rate_eval.screened_rates(k_c12_o16_to_he4_mg24);
    rate_eval.screened_rates(k_c12_o16_to_he4_mg24) *= scor;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        dratraw_dT = rate_eval.dscreened_rates_dT(k_c12_o16_to_he4_mg24);
        rate_eval.dscreened_rates_dT(k_c12_o16_to_he4_mg24) = ratraw * dscor_dt + dratraw_dT * scor;
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


    // Fill approximate rates

    fill_approx_rates<do_T_derivatives, T>(tfactors, rate_eval);

    // Calculate tabular rates

    [[maybe_unused]] Real rate, drate_dt, edot_nu, edot_gamma;

    rate_eval.enuc_weak = 0.0;

    tabular_evaluate(j_na23_ne23_meta, j_na23_ne23_rhoy, j_na23_ne23_temp, j_na23_ne23_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_na23_to_ne23) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_na23_to_ne23) = drate_dt;
    }
    rate_eval.enuc_weak += C::Legacy::n_A * Y(Na23) * (edot_nu + edot_gamma);

    tabular_evaluate(j_ne23_na23_meta, j_ne23_na23_rhoy, j_ne23_na23_temp, j_ne23_na23_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_ne23_to_na23) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_ne23_to_na23) = drate_dt;
    }
    rate_eval.enuc_weak += C::Legacy::n_A * Y(Ne23) * (edot_nu + edot_gamma);

    tabular_evaluate(j_mg23_na23_meta, j_mg23_na23_rhoy, j_mg23_na23_temp, j_mg23_na23_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_mg23_to_na23) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_mg23_to_na23) = drate_dt;
    }
    rate_eval.enuc_weak += C::Legacy::n_A * Y(Mg23) * (edot_nu + edot_gamma);

    tabular_evaluate(j_n_p_meta, j_n_p_rhoy, j_n_p_temp, j_n_p_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_n_to_p) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_n_to_p) = drate_dt;
    }
    rate_eval.enuc_weak += C::Legacy::n_A * Y(N) * (edot_nu + edot_gamma);

    tabular_evaluate(j_p_n_meta, j_p_n_rhoy, j_p_n_temp, j_p_n_data,
                     rhoy, state.T, rate, drate_dt, edot_nu, edot_gamma);
    rate_eval.screened_rates(k_p_to_n) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_to_n) = drate_dt;
    }
    rate_eval.enuc_weak += C::Legacy::n_A * Y(H1) * (edot_nu + edot_gamma);


}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rhs_nuc(const burn_t& state,
             Array1D<Real, 1, neqs>& ydot_nuc,
             const Array1D<Real, 1, NumSpec>& Y,
             const Array1D<Real, 1, NumRates>& screened_rates) {

    using namespace Rates;

    ydot_nuc(N) =
        -screened_rates(k_n_mg23_to_mg24)*Y(Mg23)*Y(N)*state.rho +
        -screened_rates(k_n_mg23_to_p_na23)*Y(Mg23)*Y(N)*state.rho +
        -screened_rates(k_n_mg23_to_he4_ne20)*Y(Mg23)*Y(N)*state.rho +
        (-screened_rates(k_n_mg23_to_c12_c12)*Y(Mg23)*Y(N)*state.rho + 0.5*screened_rates(k_c12_c12_to_n_mg23)*std::pow(Y(C12), 2)*state.rho) +
        -screened_rates(k_n_to_p)*Y(N) +
        screened_rates(k_p_to_n)*Y(H1);

    ydot_nuc(H1) =
        -screened_rates(k_p_c12_to_n13)*Y(C12)*Y(H1)*state.rho +
        -screened_rates(k_p_na23_to_mg24)*Y(Na23)*Y(H1)*state.rho +
        0.5*screened_rates(k_c12_c12_to_p_na23)*std::pow(Y(C12), 2)*state.rho +
        screened_rates(k_he4_n13_to_p_o16)*Y(He4)*Y(N13)*state.rho +
        (-screened_rates(k_p_na23_to_he4_ne20)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_he4_ne20_to_p_na23)*Y(He4)*Y(Ne20)*state.rho) +
        screened_rates(k_n_mg23_to_p_na23)*Y(Mg23)*Y(N)*state.rho +
        screened_rates(k_n_to_p)*Y(N) +
        -screened_rates(k_p_to_n)*Y(H1);

    ydot_nuc(He4) =
        -screened_rates(k_he4_c12_to_o16)*Y(C12)*Y(He4)*state.rho +
        -screened_rates(k_he4_o16_to_ne20)*Y(He4)*Y(O16)*state.rho +
        -screened_rates(k_he4_ne20_to_mg24)*Y(He4)*Y(Ne20)*state.rho +
        0.5*screened_rates(k_c12_c12_to_he4_ne20)*std::pow(Y(C12), 2)*state.rho +
        -screened_rates(k_he4_n13_to_p_o16)*Y(He4)*Y(N13)*state.rho +
        screened_rates(k_c12_o16_to_he4_mg24)*Y(C12)*Y(O16)*state.rho +
        (screened_rates(k_p_na23_to_he4_ne20)*Y(Na23)*Y(H1)*state.rho + -screened_rates(k_he4_ne20_to_p_na23)*Y(He4)*Y(Ne20)*state.rho) +
        screened_rates(k_n_mg23_to_he4_ne20)*Y(Mg23)*Y(N)*state.rho +
        -0.5*screened_rates(k_he4_he4_he4_to_c12)*std::pow(Y(He4), 3)*std::pow(state.rho, 2);

    ydot_nuc(C12) =
        -screened_rates(k_p_c12_to_n13)*Y(C12)*Y(H1)*state.rho +
        -screened_rates(k_he4_c12_to_o16)*Y(C12)*Y(He4)*state.rho +
        -screened_rates(k_c12_c12_to_p_na23)*std::pow(Y(C12), 2)*state.rho +
        -screened_rates(k_c12_c12_to_he4_ne20)*std::pow(Y(C12), 2)*state.rho +
        -screened_rates(k_c12_o16_to_he4_mg24)*Y(C12)*Y(O16)*state.rho +
        (2.0*screened_rates(k_n_mg23_to_c12_c12)*Y(Mg23)*Y(N)*state.rho + -screened_rates(k_c12_c12_to_n_mg23)*std::pow(Y(C12), 2)*state.rho) +
        0.16666666666666667*screened_rates(k_he4_he4_he4_to_c12)*std::pow(Y(He4), 3)*std::pow(state.rho, 2);

    ydot_nuc(N13) =
        screened_rates(k_p_c12_to_n13)*Y(C12)*Y(H1)*state.rho +
        -screened_rates(k_he4_n13_to_p_o16)*Y(He4)*Y(N13)*state.rho;

    ydot_nuc(O16) =
        screened_rates(k_he4_c12_to_o16)*Y(C12)*Y(He4)*state.rho +
        -screened_rates(k_he4_o16_to_ne20)*Y(He4)*Y(O16)*state.rho +
        screened_rates(k_he4_n13_to_p_o16)*Y(He4)*Y(N13)*state.rho +
        -screened_rates(k_c12_o16_to_he4_mg24)*Y(C12)*Y(O16)*state.rho;

    ydot_nuc(Ne20) =
        screened_rates(k_he4_o16_to_ne20)*Y(He4)*Y(O16)*state.rho +
        -screened_rates(k_he4_ne20_to_mg24)*Y(He4)*Y(Ne20)*state.rho +
        0.5*screened_rates(k_c12_c12_to_he4_ne20)*std::pow(Y(C12), 2)*state.rho +
        (screened_rates(k_p_na23_to_he4_ne20)*Y(Na23)*Y(H1)*state.rho + -screened_rates(k_he4_ne20_to_p_na23)*Y(He4)*Y(Ne20)*state.rho) +
        screened_rates(k_n_mg23_to_he4_ne20)*Y(Mg23)*Y(N)*state.rho;

    ydot_nuc(Ne23) =
        (-screened_rates(k_ne23_to_na23)*Y(Ne23) + screened_rates(k_na23_to_ne23)*Y(Na23));

    ydot_nuc(Na23) =
        -screened_rates(k_p_na23_to_mg24)*Y(Na23)*Y(H1)*state.rho +
        0.5*screened_rates(k_c12_c12_to_p_na23)*std::pow(Y(C12), 2)*state.rho +
        (-screened_rates(k_p_na23_to_he4_ne20)*Y(Na23)*Y(H1)*state.rho + screened_rates(k_he4_ne20_to_p_na23)*Y(He4)*Y(Ne20)*state.rho) +
        screened_rates(k_n_mg23_to_p_na23)*Y(Mg23)*Y(N)*state.rho +
        (screened_rates(k_ne23_to_na23)*Y(Ne23) + -screened_rates(k_na23_to_ne23)*Y(Na23)) +
        screened_rates(k_mg23_to_na23)*Y(Mg23);

    ydot_nuc(Mg23) =
        -screened_rates(k_n_mg23_to_mg24)*Y(Mg23)*Y(N)*state.rho +
        -screened_rates(k_n_mg23_to_p_na23)*Y(Mg23)*Y(N)*state.rho +
        -screened_rates(k_n_mg23_to_he4_ne20)*Y(Mg23)*Y(N)*state.rho +
        (-screened_rates(k_n_mg23_to_c12_c12)*Y(Mg23)*Y(N)*state.rho + 0.5*screened_rates(k_c12_c12_to_n_mg23)*std::pow(Y(C12), 2)*state.rho) +
        -screened_rates(k_mg23_to_na23)*Y(Mg23);

    ydot_nuc(Mg24) =
        screened_rates(k_he4_ne20_to_mg24)*Y(He4)*Y(Ne20)*state.rho +
        screened_rates(k_p_na23_to_mg24)*Y(Na23)*Y(H1)*state.rho +
        screened_rates(k_n_mg23_to_mg24)*Y(Mg23)*Y(N)*state.rho +
        screened_rates(k_c12_o16_to_he4_mg24)*Y(C12)*Y(O16)*state.rho;

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

    scratch = -screened_rates(k_n_mg23_to_c12_c12)*Y(Mg23)*state.rho - screened_rates(k_n_mg23_to_he4_ne20)*Y(Mg23)*state.rho - screened_rates(k_n_mg23_to_mg24)*Y(Mg23)*state.rho - screened_rates(k_n_mg23_to_p_na23)*Y(Mg23)*state.rho - screened_rates(k_n_to_p);
    jac.set(N, N, scratch);

    scratch = screened_rates(k_p_to_n);
    jac.set(N, H1, scratch);

    scratch = 1.0*screened_rates(k_c12_c12_to_n_mg23)*Y(C12)*state.rho;
    jac.set(N, C12, scratch);

    scratch = -screened_rates(k_n_mg23_to_c12_c12)*Y(N)*state.rho - screened_rates(k_n_mg23_to_he4_ne20)*Y(N)*state.rho - screened_rates(k_n_mg23_to_mg24)*Y(N)*state.rho - screened_rates(k_n_mg23_to_p_na23)*Y(N)*state.rho;
    jac.set(N, Mg23, scratch);

    scratch = screened_rates(k_n_mg23_to_p_na23)*Y(Mg23)*state.rho + screened_rates(k_n_to_p);
    jac.set(H1, N, scratch);

    scratch = -screened_rates(k_p_c12_to_n13)*Y(C12)*state.rho - screened_rates(k_p_na23_to_he4_ne20)*Y(Na23)*state.rho - screened_rates(k_p_na23_to_mg24)*Y(Na23)*state.rho - screened_rates(k_p_to_n);
    jac.set(H1, H1, scratch);

    scratch = screened_rates(k_he4_n13_to_p_o16)*Y(N13)*state.rho + screened_rates(k_he4_ne20_to_p_na23)*Y(Ne20)*state.rho;
    jac.set(H1, He4, scratch);

    scratch = 1.0*screened_rates(k_c12_c12_to_p_na23)*Y(C12)*state.rho - screened_rates(k_p_c12_to_n13)*Y(H1)*state.rho;
    jac.set(H1, C12, scratch);

    scratch = screened_rates(k_he4_n13_to_p_o16)*Y(He4)*state.rho;
    jac.set(H1, N13, scratch);

    scratch = screened_rates(k_he4_ne20_to_p_na23)*Y(He4)*state.rho;
    jac.set(H1, Ne20, scratch);

    scratch = -screened_rates(k_p_na23_to_he4_ne20)*Y(H1)*state.rho - screened_rates(k_p_na23_to_mg24)*Y(H1)*state.rho;
    jac.set(H1, Na23, scratch);

    scratch = screened_rates(k_n_mg23_to_p_na23)*Y(N)*state.rho;
    jac.set(H1, Mg23, scratch);

    scratch = screened_rates(k_n_mg23_to_he4_ne20)*Y(Mg23)*state.rho;
    jac.set(He4, N, scratch);

    scratch = screened_rates(k_p_na23_to_he4_ne20)*Y(Na23)*state.rho;
    jac.set(He4, H1, scratch);

    scratch = -screened_rates(k_he4_c12_to_o16)*Y(C12)*state.rho - 1.5*screened_rates(k_he4_he4_he4_to_c12)*std::pow(Y(He4), 2)*std::pow(state.rho, 2) - screened_rates(k_he4_n13_to_p_o16)*Y(N13)*state.rho - screened_rates(k_he4_ne20_to_mg24)*Y(Ne20)*state.rho - screened_rates(k_he4_ne20_to_p_na23)*Y(Ne20)*state.rho - screened_rates(k_he4_o16_to_ne20)*Y(O16)*state.rho;
    jac.set(He4, He4, scratch);

    scratch = 1.0*screened_rates(k_c12_c12_to_he4_ne20)*Y(C12)*state.rho + screened_rates(k_c12_o16_to_he4_mg24)*Y(O16)*state.rho - screened_rates(k_he4_c12_to_o16)*Y(He4)*state.rho;
    jac.set(He4, C12, scratch);

    scratch = -screened_rates(k_he4_n13_to_p_o16)*Y(He4)*state.rho;
    jac.set(He4, N13, scratch);

    scratch = screened_rates(k_c12_o16_to_he4_mg24)*Y(C12)*state.rho - screened_rates(k_he4_o16_to_ne20)*Y(He4)*state.rho;
    jac.set(He4, O16, scratch);

    scratch = -screened_rates(k_he4_ne20_to_mg24)*Y(He4)*state.rho - screened_rates(k_he4_ne20_to_p_na23)*Y(He4)*state.rho;
    jac.set(He4, Ne20, scratch);

    scratch = screened_rates(k_p_na23_to_he4_ne20)*Y(H1)*state.rho;
    jac.set(He4, Na23, scratch);

    scratch = screened_rates(k_n_mg23_to_he4_ne20)*Y(N)*state.rho;
    jac.set(He4, Mg23, scratch);

    scratch = 2.0*screened_rates(k_n_mg23_to_c12_c12)*Y(Mg23)*state.rho;
    jac.set(C12, N, scratch);

    scratch = -screened_rates(k_p_c12_to_n13)*Y(C12)*state.rho;
    jac.set(C12, H1, scratch);

    scratch = -screened_rates(k_he4_c12_to_o16)*Y(C12)*state.rho + 0.5*screened_rates(k_he4_he4_he4_to_c12)*std::pow(Y(He4), 2)*std::pow(state.rho, 2);
    jac.set(C12, He4, scratch);

    scratch = -2.0*screened_rates(k_c12_c12_to_he4_ne20)*Y(C12)*state.rho - 2.0*screened_rates(k_c12_c12_to_n_mg23)*Y(C12)*state.rho - 2.0*screened_rates(k_c12_c12_to_p_na23)*Y(C12)*state.rho - screened_rates(k_c12_o16_to_he4_mg24)*Y(O16)*state.rho - screened_rates(k_he4_c12_to_o16)*Y(He4)*state.rho - screened_rates(k_p_c12_to_n13)*Y(H1)*state.rho;
    jac.set(C12, C12, scratch);

    scratch = -screened_rates(k_c12_o16_to_he4_mg24)*Y(C12)*state.rho;
    jac.set(C12, O16, scratch);

    scratch = 2.0*screened_rates(k_n_mg23_to_c12_c12)*Y(N)*state.rho;
    jac.set(C12, Mg23, scratch);

    scratch = screened_rates(k_p_c12_to_n13)*Y(C12)*state.rho;
    jac.set(N13, H1, scratch);

    scratch = -screened_rates(k_he4_n13_to_p_o16)*Y(N13)*state.rho;
    jac.set(N13, He4, scratch);

    scratch = screened_rates(k_p_c12_to_n13)*Y(H1)*state.rho;
    jac.set(N13, C12, scratch);

    scratch = -screened_rates(k_he4_n13_to_p_o16)*Y(He4)*state.rho;
    jac.set(N13, N13, scratch);

    scratch = screened_rates(k_he4_c12_to_o16)*Y(C12)*state.rho + screened_rates(k_he4_n13_to_p_o16)*Y(N13)*state.rho - screened_rates(k_he4_o16_to_ne20)*Y(O16)*state.rho;
    jac.set(O16, He4, scratch);

    scratch = -screened_rates(k_c12_o16_to_he4_mg24)*Y(O16)*state.rho + screened_rates(k_he4_c12_to_o16)*Y(He4)*state.rho;
    jac.set(O16, C12, scratch);

    scratch = screened_rates(k_he4_n13_to_p_o16)*Y(He4)*state.rho;
    jac.set(O16, N13, scratch);

    scratch = -screened_rates(k_c12_o16_to_he4_mg24)*Y(C12)*state.rho - screened_rates(k_he4_o16_to_ne20)*Y(He4)*state.rho;
    jac.set(O16, O16, scratch);

    scratch = screened_rates(k_n_mg23_to_he4_ne20)*Y(Mg23)*state.rho;
    jac.set(Ne20, N, scratch);

    scratch = screened_rates(k_p_na23_to_he4_ne20)*Y(Na23)*state.rho;
    jac.set(Ne20, H1, scratch);

    scratch = -screened_rates(k_he4_ne20_to_mg24)*Y(Ne20)*state.rho - screened_rates(k_he4_ne20_to_p_na23)*Y(Ne20)*state.rho + screened_rates(k_he4_o16_to_ne20)*Y(O16)*state.rho;
    jac.set(Ne20, He4, scratch);

    scratch = 1.0*screened_rates(k_c12_c12_to_he4_ne20)*Y(C12)*state.rho;
    jac.set(Ne20, C12, scratch);

    scratch = screened_rates(k_he4_o16_to_ne20)*Y(He4)*state.rho;
    jac.set(Ne20, O16, scratch);

    scratch = -screened_rates(k_he4_ne20_to_mg24)*Y(He4)*state.rho - screened_rates(k_he4_ne20_to_p_na23)*Y(He4)*state.rho;
    jac.set(Ne20, Ne20, scratch);

    scratch = screened_rates(k_p_na23_to_he4_ne20)*Y(H1)*state.rho;
    jac.set(Ne20, Na23, scratch);

    scratch = screened_rates(k_n_mg23_to_he4_ne20)*Y(N)*state.rho;
    jac.set(Ne20, Mg23, scratch);

    scratch = -screened_rates(k_ne23_to_na23);
    jac.set(Ne23, Ne23, scratch);

    scratch = screened_rates(k_na23_to_ne23);
    jac.set(Ne23, Na23, scratch);

    scratch = screened_rates(k_n_mg23_to_p_na23)*Y(Mg23)*state.rho;
    jac.set(Na23, N, scratch);

    scratch = -screened_rates(k_p_na23_to_he4_ne20)*Y(Na23)*state.rho - screened_rates(k_p_na23_to_mg24)*Y(Na23)*state.rho;
    jac.set(Na23, H1, scratch);

    scratch = screened_rates(k_he4_ne20_to_p_na23)*Y(Ne20)*state.rho;
    jac.set(Na23, He4, scratch);

    scratch = 1.0*screened_rates(k_c12_c12_to_p_na23)*Y(C12)*state.rho;
    jac.set(Na23, C12, scratch);

    scratch = screened_rates(k_he4_ne20_to_p_na23)*Y(He4)*state.rho;
    jac.set(Na23, Ne20, scratch);

    scratch = screened_rates(k_ne23_to_na23);
    jac.set(Na23, Ne23, scratch);

    scratch = -screened_rates(k_na23_to_ne23) - screened_rates(k_p_na23_to_he4_ne20)*Y(H1)*state.rho - screened_rates(k_p_na23_to_mg24)*Y(H1)*state.rho;
    jac.set(Na23, Na23, scratch);

    scratch = screened_rates(k_mg23_to_na23) + screened_rates(k_n_mg23_to_p_na23)*Y(N)*state.rho;
    jac.set(Na23, Mg23, scratch);

    scratch = -screened_rates(k_n_mg23_to_c12_c12)*Y(Mg23)*state.rho - screened_rates(k_n_mg23_to_he4_ne20)*Y(Mg23)*state.rho - screened_rates(k_n_mg23_to_mg24)*Y(Mg23)*state.rho - screened_rates(k_n_mg23_to_p_na23)*Y(Mg23)*state.rho;
    jac.set(Mg23, N, scratch);

    scratch = 1.0*screened_rates(k_c12_c12_to_n_mg23)*Y(C12)*state.rho;
    jac.set(Mg23, C12, scratch);

    scratch = -screened_rates(k_mg23_to_na23) - screened_rates(k_n_mg23_to_c12_c12)*Y(N)*state.rho - screened_rates(k_n_mg23_to_he4_ne20)*Y(N)*state.rho - screened_rates(k_n_mg23_to_mg24)*Y(N)*state.rho - screened_rates(k_n_mg23_to_p_na23)*Y(N)*state.rho;
    jac.set(Mg23, Mg23, scratch);

    scratch = screened_rates(k_n_mg23_to_mg24)*Y(Mg23)*state.rho;
    jac.set(Mg24, N, scratch);

    scratch = screened_rates(k_p_na23_to_mg24)*Y(Na23)*state.rho;
    jac.set(Mg24, H1, scratch);

    scratch = screened_rates(k_he4_ne20_to_mg24)*Y(Ne20)*state.rho;
    jac.set(Mg24, He4, scratch);

    scratch = screened_rates(k_c12_o16_to_he4_mg24)*Y(O16)*state.rho;
    jac.set(Mg24, C12, scratch);

    scratch = screened_rates(k_c12_o16_to_he4_mg24)*Y(C12)*state.rho;
    jac.set(Mg24, O16, scratch);

    scratch = screened_rates(k_he4_ne20_to_mg24)*Y(He4)*state.rho;
    jac.set(Mg24, Ne20, scratch);

    scratch = screened_rates(k_p_na23_to_mg24)*Y(H1)*state.rho;
    jac.set(Mg24, Na23, scratch);

    scratch = screened_rates(k_n_mg23_to_mg24)*Y(N)*state.rho;
    jac.set(Mg24, Mg23, scratch);


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
