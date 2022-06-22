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
#include <extern_parameters.H>

using namespace amrex;
using namespace ArrayUtil;

using namespace Species;
using namespace Rates;

using namespace rate_tables;

struct rate_eval_t {
    Array1D<Real, 1, NumRates>  screened_rates;
    Array1D<Real, 1, NumRates>  dscreened_rates_dT;
    Array1D<Real, NrateReaclib+1, NrateReaclib+NrateTabular> add_energy_rate;
};



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


AMREX_GPU_HOST_DEVICE AMREX_INLINE
void evaluate_rates(const burn_t& state, rate_eval_t& rate_eval) {


    // create molar fractions

    Array1D<Real, 1, NumSpec> Y;
    for (int n = 1; n <= NumSpec; ++n) {
        Y(n) = state.xn[n-1] * aion_inv[n-1];
    }

    Real rhoy = state.rho * state.y_e;

    // Calculate Reaclib rates

    plasma_state_t pstate;
    fill_plasma_state(pstate, state.T, state.rho, Y);

    tf_t tfactors = evaluate_tfactors(state.T);

    for (int i = 1; i <= NrateReaclib; ++i) {


        if (jacobian == 1) {
            Real rate;
            Real drate_dt;
            constexpr int do_T_derivatives = 1;
            reaclib_evaluate<do_T_derivatives>(pstate, tfactors, i, rate, drate_dt);
            rate_eval.screened_rates(i) = rate;
            rate_eval.dscreened_rates_dT(i) = drate_dt;
        } else {
            Real rate;
            Real drate_dt;
            // drate_dt is not used in this case
            constexpr int do_T_derivatives = 0;
            reaclib_evaluate<do_T_derivatives>(pstate, tfactors, i, rate, drate_dt);
            rate_eval.screened_rates(i) = rate;
            rate_eval.dscreened_rates_dT(i) = drate_dt;
        }
    }

    // Evaluate screening factors

    Real ratraw, dratraw_dT;
    Real scor, dscor_dt, dscor_dd;
    Real scor2, dscor2_dt, dscor2_dd;


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


    ratraw = rate_eval.screened_rates(k_he4_he4_he4__c12);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_he4_he4__c12);
    rate_eval.screened_rates(k_he4_he4_he4__c12) *= scor * scor2;
    rate_eval.dscreened_rates_dT(k_he4_he4_he4__c12) = ratraw * (scor * dscor2_dt + dscor_dt * scor2) + dratraw_dT * scor * scor2;


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 6.0_rt, 12.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_c12__o16);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_c12__o16);
    rate_eval.screened_rates(k_he4_c12__o16) *= scor;
    rate_eval.dscreened_rates_dT(k_he4_c12__o16) = ratraw * dscor_dt + dratraw_dT * scor;


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 7.0_rt, 14.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_n14__f18);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_n14__f18);
    rate_eval.screened_rates(k_he4_n14__f18) *= scor;
    rate_eval.dscreened_rates_dT(k_he4_n14__f18) = ratraw * dscor_dt + dratraw_dT * scor;


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 9.0_rt, 18.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_f18__p_ne21);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_f18__p_ne21);
    rate_eval.screened_rates(k_he4_f18__p_ne21) *= scor;
    rate_eval.dscreened_rates_dT(k_he4_f18__p_ne21) = ratraw * dscor_dt + dratraw_dT * scor;


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 6.0_rt, 12.0_rt);


        static_assert(scn_fac.z1 == 1.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_p_c12__n13);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_p_c12__n13);
    rate_eval.screened_rates(k_p_c12__n13) *= scor;
    rate_eval.dscreened_rates_dT(k_p_c12__n13) = ratraw * dscor_dt + dratraw_dT * scor;


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 7.0_rt, 13.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_n13__p_o16);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_n13__p_o16);
    rate_eval.screened_rates(k_he4_n13__p_o16) *= scor;
    rate_eval.dscreened_rates_dT(k_he4_n13__p_o16) = ratraw * dscor_dt + dratraw_dT * scor;


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 8.0_rt, 16.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_o16__ne20);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_o16__ne20);
    rate_eval.screened_rates(k_he4_o16__ne20) *= scor;
    rate_eval.dscreened_rates_dT(k_he4_o16__ne20) = ratraw * dscor_dt + dratraw_dT * scor;


    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 6.0_rt, 14.0_rt);


        static_assert(scn_fac.z1 == 2.0_rt);


        actual_screen(pstate, scn_fac, scor, dscor_dt);
    }


    ratraw = rate_eval.screened_rates(k_he4_c14__o18);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_c14__o18);
    rate_eval.screened_rates(k_he4_c14__o18) *= scor;
    rate_eval.dscreened_rates_dT(k_he4_c14__o18) = ratraw * dscor_dt + dratraw_dT * scor;


    // Calculate tabular rates

    Real rate, drate_dt, edot_nu;


}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rhs_nuc(const burn_t& state,
             Array1D<Real, 1, neqs>& ydot_nuc,
             const Array1D<Real, 1, NumSpec>& Y,
             const Array1D<Real, 1, NumRates>& screened_rates) {

    using namespace Rates;

    ydot_nuc(H1) =
        screened_rates(k_he4_f18__p_ne21)*Y(F18)*Y(He4)*state.rho +
        -screened_rates(k_p_c12__n13)*Y(C12)*Y(H1)*state.rho +
        screened_rates(k_he4_n13__p_o16)*Y(He4)*Y(N13)*state.rho;

    ydot_nuc(He4) =
        -0.5*screened_rates(k_he4_he4_he4__c12)*std::pow(Y(He4), 3)*std::pow(state.rho, 2) +
        -screened_rates(k_he4_c12__o16)*Y(C12)*Y(He4)*state.rho +
        -screened_rates(k_he4_n14__f18)*Y(He4)*Y(N14)*state.rho +
        -screened_rates(k_he4_f18__p_ne21)*Y(F18)*Y(He4)*state.rho +
        -screened_rates(k_he4_n13__p_o16)*Y(He4)*Y(N13)*state.rho +
        -screened_rates(k_he4_o16__ne20)*Y(He4)*Y(O16)*state.rho +
        -screened_rates(k_he4_c14__o18)*Y(C14)*Y(He4)*state.rho;

    ydot_nuc(C12) =
        0.16666666666666667*screened_rates(k_he4_he4_he4__c12)*std::pow(Y(He4), 3)*std::pow(state.rho, 2) +
        -screened_rates(k_he4_c12__o16)*Y(C12)*Y(He4)*state.rho +
        -screened_rates(k_p_c12__n13)*Y(C12)*Y(H1)*state.rho;

    ydot_nuc(C14) =
        -screened_rates(k_he4_c14__o18)*Y(C14)*Y(He4)*state.rho;

    ydot_nuc(N13) =
        screened_rates(k_p_c12__n13)*Y(C12)*Y(H1)*state.rho +
        -screened_rates(k_he4_n13__p_o16)*Y(He4)*Y(N13)*state.rho;

    ydot_nuc(N14) =
        -screened_rates(k_he4_n14__f18)*Y(He4)*Y(N14)*state.rho;

    ydot_nuc(O16) =
        screened_rates(k_he4_c12__o16)*Y(C12)*Y(He4)*state.rho +
        screened_rates(k_he4_n13__p_o16)*Y(He4)*Y(N13)*state.rho +
        -screened_rates(k_he4_o16__ne20)*Y(He4)*Y(O16)*state.rho;

    ydot_nuc(O18) =
        screened_rates(k_he4_c14__o18)*Y(C14)*Y(He4)*state.rho;

    ydot_nuc(F18) =
        screened_rates(k_he4_n14__f18)*Y(He4)*Y(N14)*state.rho +
        -screened_rates(k_he4_f18__p_ne21)*Y(F18)*Y(He4)*state.rho;

    ydot_nuc(Ne20) =
        screened_rates(k_he4_o16__ne20)*Y(He4)*Y(O16)*state.rho;

    ydot_nuc(Ne21) =
        screened_rates(k_he4_f18__p_ne21)*Y(F18)*Y(He4)*state.rho;

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

    rate_eval_t rate_eval;

    evaluate_rates(state, rate_eval);

    rhs_nuc(state, ydot, Y, rate_eval.screened_rates);

    // ion binding energy contributions

    Real enuc;
    ener_gener_rate(ydot, enuc);

    // include reaction neutrino losses (non-thermal)

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

    scratch = -screened_rates(k_p_c12__n13)*Y(C12)*state.rho;
    jac.set(H1, H1, scratch);

    scratch = screened_rates(k_he4_f18__p_ne21)*Y(F18)*state.rho + screened_rates(k_he4_n13__p_o16)*Y(N13)*state.rho;
    jac.set(H1, He4, scratch);

    scratch = -screened_rates(k_p_c12__n13)*Y(H1)*state.rho;
    jac.set(H1, C12, scratch);

    scratch = screened_rates(k_he4_n13__p_o16)*Y(He4)*state.rho;
    jac.set(H1, N13, scratch);

    scratch = screened_rates(k_he4_f18__p_ne21)*Y(He4)*state.rho;
    jac.set(H1, F18, scratch);

    scratch = -screened_rates(k_he4_c12__o16)*Y(C12)*state.rho - screened_rates(k_he4_c14__o18)*Y(C14)*state.rho - screened_rates(k_he4_f18__p_ne21)*Y(F18)*state.rho - 1.5*screened_rates(k_he4_he4_he4__c12)*std::pow(Y(He4), 2)*std::pow(state.rho, 2) - screened_rates(k_he4_n13__p_o16)*Y(N13)*state.rho - screened_rates(k_he4_n14__f18)*Y(N14)*state.rho - screened_rates(k_he4_o16__ne20)*Y(O16)*state.rho;
    jac.set(He4, He4, scratch);

    scratch = -screened_rates(k_he4_c12__o16)*Y(He4)*state.rho;
    jac.set(He4, C12, scratch);

    scratch = -screened_rates(k_he4_c14__o18)*Y(He4)*state.rho;
    jac.set(He4, C14, scratch);

    scratch = -screened_rates(k_he4_n13__p_o16)*Y(He4)*state.rho;
    jac.set(He4, N13, scratch);

    scratch = -screened_rates(k_he4_n14__f18)*Y(He4)*state.rho;
    jac.set(He4, N14, scratch);

    scratch = -screened_rates(k_he4_o16__ne20)*Y(He4)*state.rho;
    jac.set(He4, O16, scratch);

    scratch = -screened_rates(k_he4_f18__p_ne21)*Y(He4)*state.rho;
    jac.set(He4, F18, scratch);

    scratch = -screened_rates(k_p_c12__n13)*Y(C12)*state.rho;
    jac.set(C12, H1, scratch);

    scratch = -screened_rates(k_he4_c12__o16)*Y(C12)*state.rho + 0.5*screened_rates(k_he4_he4_he4__c12)*std::pow(Y(He4), 2)*std::pow(state.rho, 2);
    jac.set(C12, He4, scratch);

    scratch = -screened_rates(k_he4_c12__o16)*Y(He4)*state.rho - screened_rates(k_p_c12__n13)*Y(H1)*state.rho;
    jac.set(C12, C12, scratch);

    scratch = -screened_rates(k_he4_c14__o18)*Y(C14)*state.rho;
    jac.set(C14, He4, scratch);

    scratch = -screened_rates(k_he4_c14__o18)*Y(He4)*state.rho;
    jac.set(C14, C14, scratch);

    scratch = screened_rates(k_p_c12__n13)*Y(C12)*state.rho;
    jac.set(N13, H1, scratch);

    scratch = -screened_rates(k_he4_n13__p_o16)*Y(N13)*state.rho;
    jac.set(N13, He4, scratch);

    scratch = screened_rates(k_p_c12__n13)*Y(H1)*state.rho;
    jac.set(N13, C12, scratch);

    scratch = -screened_rates(k_he4_n13__p_o16)*Y(He4)*state.rho;
    jac.set(N13, N13, scratch);

    scratch = -screened_rates(k_he4_n14__f18)*Y(N14)*state.rho;
    jac.set(N14, He4, scratch);

    scratch = -screened_rates(k_he4_n14__f18)*Y(He4)*state.rho;
    jac.set(N14, N14, scratch);

    scratch = screened_rates(k_he4_c12__o16)*Y(C12)*state.rho + screened_rates(k_he4_n13__p_o16)*Y(N13)*state.rho - screened_rates(k_he4_o16__ne20)*Y(O16)*state.rho;
    jac.set(O16, He4, scratch);

    scratch = screened_rates(k_he4_c12__o16)*Y(He4)*state.rho;
    jac.set(O16, C12, scratch);

    scratch = screened_rates(k_he4_n13__p_o16)*Y(He4)*state.rho;
    jac.set(O16, N13, scratch);

    scratch = -screened_rates(k_he4_o16__ne20)*Y(He4)*state.rho;
    jac.set(O16, O16, scratch);

    scratch = screened_rates(k_he4_c14__o18)*Y(C14)*state.rho;
    jac.set(O18, He4, scratch);

    scratch = screened_rates(k_he4_c14__o18)*Y(He4)*state.rho;
    jac.set(O18, C14, scratch);

    scratch = -screened_rates(k_he4_f18__p_ne21)*Y(F18)*state.rho + screened_rates(k_he4_n14__f18)*Y(N14)*state.rho;
    jac.set(F18, He4, scratch);

    scratch = screened_rates(k_he4_n14__f18)*Y(He4)*state.rho;
    jac.set(F18, N14, scratch);

    scratch = -screened_rates(k_he4_f18__p_ne21)*Y(He4)*state.rho;
    jac.set(F18, F18, scratch);

    scratch = screened_rates(k_he4_o16__ne20)*Y(O16)*state.rho;
    jac.set(Ne20, He4, scratch);

    scratch = screened_rates(k_he4_o16__ne20)*Y(He4)*state.rho;
    jac.set(Ne20, O16, scratch);

    scratch = screened_rates(k_he4_f18__p_ne21)*Y(F18)*state.rho;
    jac.set(Ne21, He4, scratch);

    scratch = screened_rates(k_he4_f18__p_ne21)*Y(He4)*state.rho;
    jac.set(Ne21, F18, scratch);


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

    rate_eval_t rate_eval;

    evaluate_rates(state, rate_eval);

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

    init_reaclib();

    init_tabular();

}


#endif
