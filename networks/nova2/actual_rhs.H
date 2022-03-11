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

    if (jacobian == 1) {
        for (int i = 1; i <= NrateReaclib; ++i) {
            Real rate;
            Real drate_dt;
            constexpr int do_T_derivatives = 1;
            reaclib_evaluate<do_T_derivatives>(pstate, tfactors, i, rate, drate_dt);
            rate_eval.screened_rates(i) = rate;
            rate_eval.dscreened_rates_dT(i) = drate_dt;
        }
    } else {
        for (int i = 1; i <= NrateReaclib; ++i) {
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


    screen5(pstate, 0, zion[H1-1], aion[H1-1], zion[H1-1], aion[H1-1], scor, dscor_dt, dscor_dd);


    ratraw = rate_eval.screened_rates(k_p_p__d__weak__bet_pos_);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_p_p__d__weak__bet_pos_);
    rate_eval.screened_rates(k_p_p__d__weak__bet_pos_) *= scor;
    rate_eval.dscreened_rates_dT(k_p_p__d__weak__bet_pos_) = ratraw * dscor_dt + dratraw_dT * scor;

    ratraw = rate_eval.screened_rates(k_p_p__d__weak__electron_capture);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_p_p__d__weak__electron_capture);
    rate_eval.screened_rates(k_p_p__d__weak__electron_capture) *= scor;
    rate_eval.dscreened_rates_dT(k_p_p__d__weak__electron_capture) = ratraw * dscor_dt + dratraw_dT * scor;


    screen5(pstate, 1, zion[H1-1], aion[H1-1], zion[H2-1], aion[H2-1], scor, dscor_dt, dscor_dd);


    ratraw = rate_eval.screened_rates(k_p_d__he3);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_p_d__he3);
    rate_eval.screened_rates(k_p_d__he3) *= scor;
    rate_eval.dscreened_rates_dT(k_p_d__he3) = ratraw * dscor_dt + dratraw_dT * scor;


    screen5(pstate, 2, zion[H2-1], aion[H2-1], zion[H2-1], aion[H2-1], scor, dscor_dt, dscor_dd);


    ratraw = rate_eval.screened_rates(k_d_d__he4);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_d_d__he4);
    rate_eval.screened_rates(k_d_d__he4) *= scor;
    rate_eval.dscreened_rates_dT(k_d_d__he4) = ratraw * dscor_dt + dratraw_dT * scor;


    screen5(pstate, 3, zion[H1-1], aion[H1-1], zion[He3-1], aion[He3-1], scor, dscor_dt, dscor_dd);


    ratraw = rate_eval.screened_rates(k_p_he3__he4__weak__bet_pos_);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_p_he3__he4__weak__bet_pos_);
    rate_eval.screened_rates(k_p_he3__he4__weak__bet_pos_) *= scor;
    rate_eval.dscreened_rates_dT(k_p_he3__he4__weak__bet_pos_) = ratraw * dscor_dt + dratraw_dT * scor;


    screen5(pstate, 4, zion[He4-1], aion[He4-1], zion[He3-1], aion[He3-1], scor, dscor_dt, dscor_dd);


    ratraw = rate_eval.screened_rates(k_he4_he3__be7);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_he3__be7);
    rate_eval.screened_rates(k_he4_he3__be7) *= scor;
    rate_eval.dscreened_rates_dT(k_he4_he3__be7) = ratraw * dscor_dt + dratraw_dT * scor;


    screen5(pstate, 5, zion[H1-1], aion[H1-1], zion[Be7-1], aion[Be7-1], scor, dscor_dt, dscor_dd);


    ratraw = rate_eval.screened_rates(k_p_be7__b8);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_p_be7__b8);
    rate_eval.screened_rates(k_p_be7__b8) *= scor;
    rate_eval.dscreened_rates_dT(k_p_be7__b8) = ratraw * dscor_dt + dratraw_dT * scor;


    screen5(pstate, 6, zion[H1-1], aion[H1-1], zion[C12-1], aion[C12-1], scor, dscor_dt, dscor_dd);


    ratraw = rate_eval.screened_rates(k_p_c12__n13);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_p_c12__n13);
    rate_eval.screened_rates(k_p_c12__n13) *= scor;
    rate_eval.dscreened_rates_dT(k_p_c12__n13) = ratraw * dscor_dt + dratraw_dT * scor;


    screen5(pstate, 7, zion[He4-1], aion[He4-1], zion[C12-1], aion[C12-1], scor, dscor_dt, dscor_dd);


    ratraw = rate_eval.screened_rates(k_he4_c12__o16);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_c12__o16);
    rate_eval.screened_rates(k_he4_c12__o16) *= scor;
    rate_eval.dscreened_rates_dT(k_he4_c12__o16) = ratraw * dscor_dt + dratraw_dT * scor;


    screen5(pstate, 8, zion[H1-1], aion[H1-1], zion[C13-1], aion[C13-1], scor, dscor_dt, dscor_dd);


    ratraw = rate_eval.screened_rates(k_p_c13__n14);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_p_c13__n14);
    rate_eval.screened_rates(k_p_c13__n14) *= scor;
    rate_eval.dscreened_rates_dT(k_p_c13__n14) = ratraw * dscor_dt + dratraw_dT * scor;


    screen5(pstate, 9, zion[H1-1], aion[H1-1], zion[N13-1], aion[N13-1], scor, dscor_dt, dscor_dd);


    ratraw = rate_eval.screened_rates(k_p_n13__o14);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_p_n13__o14);
    rate_eval.screened_rates(k_p_n13__o14) *= scor;
    rate_eval.dscreened_rates_dT(k_p_n13__o14) = ratraw * dscor_dt + dratraw_dT * scor;


    screen5(pstate, 10, zion[H1-1], aion[H1-1], zion[N14-1], aion[N14-1], scor, dscor_dt, dscor_dd);


    ratraw = rate_eval.screened_rates(k_p_n14__o15);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_p_n14__o15);
    rate_eval.screened_rates(k_p_n14__o15) *= scor;
    rate_eval.dscreened_rates_dT(k_p_n14__o15) = ratraw * dscor_dt + dratraw_dT * scor;


    screen5(pstate, 11, zion[He4-1], aion[He4-1], zion[N14-1], aion[N14-1], scor, dscor_dt, dscor_dd);


    ratraw = rate_eval.screened_rates(k_he4_n14__f18);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_n14__f18);
    rate_eval.screened_rates(k_he4_n14__f18) *= scor;
    rate_eval.dscreened_rates_dT(k_he4_n14__f18) = ratraw * dscor_dt + dratraw_dT * scor;


    screen5(pstate, 12, zion[H1-1], aion[H1-1], zion[N15-1], aion[N15-1], scor, dscor_dt, dscor_dd);


    ratraw = rate_eval.screened_rates(k_p_n15__o16);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_p_n15__o16);
    rate_eval.screened_rates(k_p_n15__o16) *= scor;
    rate_eval.dscreened_rates_dT(k_p_n15__o16) = ratraw * dscor_dt + dratraw_dT * scor;

    ratraw = rate_eval.screened_rates(k_p_n15__he4_c12);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_p_n15__he4_c12);
    rate_eval.screened_rates(k_p_n15__he4_c12) *= scor;
    rate_eval.dscreened_rates_dT(k_p_n15__he4_c12) = ratraw * dscor_dt + dratraw_dT * scor;


    screen5(pstate, 13, zion[H1-1], aion[H1-1], zion[O16-1], aion[O16-1], scor, dscor_dt, dscor_dd);


    ratraw = rate_eval.screened_rates(k_p_o16__f17);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_p_o16__f17);
    rate_eval.screened_rates(k_p_o16__f17) *= scor;
    rate_eval.dscreened_rates_dT(k_p_o16__f17) = ratraw * dscor_dt + dratraw_dT * scor;


    screen5(pstate, 14, zion[H1-1], aion[H1-1], zion[O17-1], aion[O17-1], scor, dscor_dt, dscor_dd);


    ratraw = rate_eval.screened_rates(k_p_o17__f18);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_p_o17__f18);
    rate_eval.screened_rates(k_p_o17__f18) *= scor;
    rate_eval.dscreened_rates_dT(k_p_o17__f18) = ratraw * dscor_dt + dratraw_dT * scor;

    ratraw = rate_eval.screened_rates(k_p_o17__he4_n14);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_p_o17__he4_n14);
    rate_eval.screened_rates(k_p_o17__he4_n14) *= scor;
    rate_eval.dscreened_rates_dT(k_p_o17__he4_n14) = ratraw * dscor_dt + dratraw_dT * scor;


    screen5(pstate, 15, zion[H2-1], aion[H2-1], zion[He3-1], aion[He3-1], scor, dscor_dt, dscor_dd);


    ratraw = rate_eval.screened_rates(k_d_he3__p_he4);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_d_he3__p_he4);
    rate_eval.screened_rates(k_d_he3__p_he4) *= scor;
    rate_eval.dscreened_rates_dT(k_d_he3__p_he4) = ratraw * dscor_dt + dratraw_dT * scor;


    screen5(pstate, 16, zion[He4-1], aion[He4-1], zion[N13-1], aion[N13-1], scor, dscor_dt, dscor_dd);


    ratraw = rate_eval.screened_rates(k_he4_n13__p_o16);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_n13__p_o16);
    rate_eval.screened_rates(k_he4_n13__p_o16) *= scor;
    rate_eval.dscreened_rates_dT(k_he4_n13__p_o16) = ratraw * dscor_dt + dratraw_dT * scor;


    screen5(pstate, 17, zion[He4-1], aion[He4-1], zion[O14-1], aion[O14-1], scor, dscor_dt, dscor_dd);


    ratraw = rate_eval.screened_rates(k_he4_o14__p_f17);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_o14__p_f17);
    rate_eval.screened_rates(k_he4_o14__p_f17) *= scor;
    rate_eval.dscreened_rates_dT(k_he4_o14__p_f17) = ratraw * dscor_dt + dratraw_dT * scor;


    screen5(pstate, 18, zion[H1-1], aion[H1-1], zion[F18-1], aion[F18-1], scor, dscor_dt, dscor_dd);


    ratraw = rate_eval.screened_rates(k_p_f18__he4_o15);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_p_f18__he4_o15);
    rate_eval.screened_rates(k_p_f18__he4_o15) *= scor;
    rate_eval.dscreened_rates_dT(k_p_f18__he4_o15) = ratraw * dscor_dt + dratraw_dT * scor;


    screen5(pstate, 19, zion[He3-1], aion[He3-1], zion[He3-1], aion[He3-1], scor, dscor_dt, dscor_dd);


    ratraw = rate_eval.screened_rates(k_he3_he3__p_p_he4);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_he3_he3__p_p_he4);
    rate_eval.screened_rates(k_he3_he3__p_p_he4) *= scor;
    rate_eval.dscreened_rates_dT(k_he3_he3__p_p_he4) = ratraw * dscor_dt + dratraw_dT * scor;


    screen5(pstate, 20, zion[H2-1], aion[H2-1], zion[Be7-1], aion[Be7-1], scor, dscor_dt, dscor_dd);


    ratraw = rate_eval.screened_rates(k_d_be7__p_he4_he4);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_d_be7__p_he4_he4);
    rate_eval.screened_rates(k_d_be7__p_he4_he4) *= scor;
    rate_eval.dscreened_rates_dT(k_d_be7__p_he4_he4) = ratraw * dscor_dt + dratraw_dT * scor;


    screen5(pstate, 21, zion[He3-1], aion[He3-1], zion[Be7-1], aion[Be7-1], scor, dscor_dt, dscor_dd);


    ratraw = rate_eval.screened_rates(k_he3_be7__p_p_he4_he4);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_he3_be7__p_p_he4_he4);
    rate_eval.screened_rates(k_he3_be7__p_p_he4_he4) *= scor;
    rate_eval.dscreened_rates_dT(k_he3_be7__p_p_he4_he4) = ratraw * dscor_dt + dratraw_dT * scor;

    screen5(pstate, 22, zion[He4-1], aion[He4-1], zion[He4-1], aion[He4-1], scor, dscor_dt, dscor_dd);

    screen5(pstate, 23, zion[He4-1], aion[He4-1], 4.0_rt, 8.0_rt, scor2, dscor2_dt, dscor2_dd);

    ratraw = rate_eval.screened_rates(k_he4_he4_he4__c12);
    dratraw_dT = rate_eval.dscreened_rates_dT(k_he4_he4_he4__c12);
    rate_eval.screened_rates(k_he4_he4_he4__c12) *= scor * scor2;
    rate_eval.dscreened_rates_dT(k_he4_he4_he4__c12) = ratraw * (scor * dscor2_dt + dscor_dt * scor2) + dratraw_dT * scor;


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
        -screened_rates(k_p_p__d__weak__bet_pos_)*std::pow(Y(H1), 2)*state.rho +
        -screened_rates(k_p_p__d__weak__electron_capture)*std::pow(Y(H1), 2)*std::pow(state.rho, 2)*state.y_e +
        -screened_rates(k_p_d__he3)*Y(H2)*Y(H1)*state.rho +
        -screened_rates(k_p_he3__he4__weak__bet_pos_)*Y(He3)*Y(H1)*state.rho +
        -screened_rates(k_p_be7__b8)*Y(Be7)*Y(H1)*state.rho +
        -screened_rates(k_p_c12__n13)*Y(C12)*Y(H1)*state.rho +
        -screened_rates(k_p_c13__n14)*Y(C13)*Y(H1)*state.rho +
        -screened_rates(k_p_n13__o14)*Y(N13)*Y(H1)*state.rho +
        -screened_rates(k_p_n14__o15)*Y(N14)*Y(H1)*state.rho +
        -screened_rates(k_p_n15__o16)*Y(N15)*Y(H1)*state.rho +
        -screened_rates(k_p_o16__f17)*Y(O16)*Y(H1)*state.rho +
        -screened_rates(k_p_o17__f18)*Y(O17)*Y(H1)*state.rho +
        -screened_rates(k_p_n15__he4_c12)*Y(N15)*Y(H1)*state.rho +
        -screened_rates(k_p_o17__he4_n14)*Y(O17)*Y(H1)*state.rho +
        -screened_rates(k_p_f18__he4_o15)*Y(F18)*Y(H1)*state.rho +
        screened_rates(k_d_he3__p_he4)*Y(H2)*Y(He3)*state.rho +
        screened_rates(k_he4_n13__p_o16)*Y(He4)*Y(N13)*state.rho +
        screened_rates(k_he4_o14__p_f17)*Y(He4)*Y(O14)*state.rho +
        screened_rates(k_he3_he3__p_p_he4)*std::pow(Y(He3), 2)*state.rho +
        screened_rates(k_d_be7__p_he4_he4)*Y(Be7)*Y(H2)*state.rho +
        2.0*screened_rates(k_he3_be7__p_p_he4_he4)*Y(Be7)*Y(He3)*state.rho;

    ydot_nuc(H2) =
        -screened_rates(k_p_d__he3)*Y(H2)*Y(H1)*state.rho +
        -screened_rates(k_d_d__he4)*std::pow(Y(H2), 2)*state.rho +
        -screened_rates(k_d_he3__p_he4)*Y(H2)*Y(He3)*state.rho +
        -screened_rates(k_d_be7__p_he4_he4)*Y(Be7)*Y(H2)*state.rho +
        0.5*screened_rates(k_p_p__d__weak__bet_pos_)*std::pow(Y(H1), 2)*state.rho +
        0.5*screened_rates(k_p_p__d__weak__electron_capture)*std::pow(Y(H1), 2)*std::pow(state.rho, 2)*state.y_e;

    ydot_nuc(He3) =
        -screened_rates(k_p_he3__he4__weak__bet_pos_)*Y(He3)*Y(H1)*state.rho +
        -screened_rates(k_he4_he3__be7)*Y(He3)*Y(He4)*state.rho +
        -screened_rates(k_d_he3__p_he4)*Y(H2)*Y(He3)*state.rho +
        -screened_rates(k_he3_he3__p_p_he4)*std::pow(Y(He3), 2)*state.rho +
        -screened_rates(k_he3_be7__p_p_he4_he4)*Y(Be7)*Y(He3)*state.rho +
        screened_rates(k_p_d__he3)*Y(H2)*Y(H1)*state.rho;

    ydot_nuc(He4) =
        -screened_rates(k_he4_he3__be7)*Y(He3)*Y(He4)*state.rho +
        -screened_rates(k_he4_c12__o16)*Y(C12)*Y(He4)*state.rho +
        -screened_rates(k_he4_n14__f18)*Y(He4)*Y(N14)*state.rho +
        -screened_rates(k_he4_n13__p_o16)*Y(He4)*Y(N13)*state.rho +
        -screened_rates(k_he4_o14__p_f17)*Y(He4)*Y(O14)*state.rho +
        -0.5*screened_rates(k_he4_he4_he4__c12)*std::pow(Y(He4), 3)*std::pow(state.rho, 2) +
        2.0*screened_rates(k_b8__he4_he4__weak__wc12)*Y(B8) +
        0.5*screened_rates(k_d_d__he4)*std::pow(Y(H2), 2)*state.rho +
        screened_rates(k_p_he3__he4__weak__bet_pos_)*Y(He3)*Y(H1)*state.rho +
        screened_rates(k_d_he3__p_he4)*Y(H2)*Y(He3)*state.rho +
        screened_rates(k_p_n15__he4_c12)*Y(N15)*Y(H1)*state.rho +
        screened_rates(k_p_o17__he4_n14)*Y(O17)*Y(H1)*state.rho +
        screened_rates(k_p_f18__he4_o15)*Y(F18)*Y(H1)*state.rho +
        0.5*screened_rates(k_he3_he3__p_p_he4)*std::pow(Y(He3), 2)*state.rho +
        2.0*screened_rates(k_d_be7__p_he4_he4)*Y(Be7)*Y(H2)*state.rho +
        2.0*screened_rates(k_he3_be7__p_p_he4_he4)*Y(Be7)*Y(He3)*state.rho;

    ydot_nuc(Be7) =
        -screened_rates(k_p_be7__b8)*Y(Be7)*Y(H1)*state.rho +
        -screened_rates(k_d_be7__p_he4_he4)*Y(Be7)*Y(H2)*state.rho +
        -screened_rates(k_he3_be7__p_p_he4_he4)*Y(Be7)*Y(He3)*state.rho +
        screened_rates(k_he4_he3__be7)*Y(He3)*Y(He4)*state.rho;

    ydot_nuc(B8) =
        -screened_rates(k_b8__he4_he4__weak__wc12)*Y(B8) +
        screened_rates(k_p_be7__b8)*Y(Be7)*Y(H1)*state.rho;

    ydot_nuc(C12) =
        -screened_rates(k_p_c12__n13)*Y(C12)*Y(H1)*state.rho +
        -screened_rates(k_he4_c12__o16)*Y(C12)*Y(He4)*state.rho +
        screened_rates(k_p_n15__he4_c12)*Y(N15)*Y(H1)*state.rho +
        0.16666666666666667*screened_rates(k_he4_he4_he4__c12)*std::pow(Y(He4), 3)*std::pow(state.rho, 2);

    ydot_nuc(C13) =
        -screened_rates(k_p_c13__n14)*Y(C13)*Y(H1)*state.rho +
        screened_rates(k_n13__c13__weak__wc12)*Y(N13);

    ydot_nuc(N13) =
        -screened_rates(k_n13__c13__weak__wc12)*Y(N13) +
        -screened_rates(k_p_n13__o14)*Y(N13)*Y(H1)*state.rho +
        -screened_rates(k_he4_n13__p_o16)*Y(He4)*Y(N13)*state.rho +
        screened_rates(k_p_c12__n13)*Y(C12)*Y(H1)*state.rho;

    ydot_nuc(N14) =
        -screened_rates(k_p_n14__o15)*Y(N14)*Y(H1)*state.rho +
        -screened_rates(k_he4_n14__f18)*Y(He4)*Y(N14)*state.rho +
        screened_rates(k_o14__n14__weak__wc12)*Y(O14) +
        screened_rates(k_p_c13__n14)*Y(C13)*Y(H1)*state.rho +
        screened_rates(k_p_o17__he4_n14)*Y(O17)*Y(H1)*state.rho;

    ydot_nuc(N15) =
        -screened_rates(k_p_n15__o16)*Y(N15)*Y(H1)*state.rho +
        -screened_rates(k_p_n15__he4_c12)*Y(N15)*Y(H1)*state.rho +
        screened_rates(k_o15__n15__weak__wc12)*Y(O15);

    ydot_nuc(O14) =
        -screened_rates(k_o14__n14__weak__wc12)*Y(O14) +
        -screened_rates(k_he4_o14__p_f17)*Y(He4)*Y(O14)*state.rho +
        screened_rates(k_p_n13__o14)*Y(N13)*Y(H1)*state.rho;

    ydot_nuc(O15) =
        -screened_rates(k_o15__n15__weak__wc12)*Y(O15) +
        screened_rates(k_p_n14__o15)*Y(N14)*Y(H1)*state.rho +
        screened_rates(k_p_f18__he4_o15)*Y(F18)*Y(H1)*state.rho;

    ydot_nuc(O16) =
        -screened_rates(k_p_o16__f17)*Y(O16)*Y(H1)*state.rho +
        screened_rates(k_he4_c12__o16)*Y(C12)*Y(He4)*state.rho +
        screened_rates(k_p_n15__o16)*Y(N15)*Y(H1)*state.rho +
        screened_rates(k_he4_n13__p_o16)*Y(He4)*Y(N13)*state.rho;

    ydot_nuc(O17) =
        -screened_rates(k_p_o17__f18)*Y(O17)*Y(H1)*state.rho +
        -screened_rates(k_p_o17__he4_n14)*Y(O17)*Y(H1)*state.rho +
        screened_rates(k_f17__o17__weak__wc12)*Y(F17);

    ydot_nuc(F17) =
        -screened_rates(k_f17__o17__weak__wc12)*Y(F17) +
        screened_rates(k_p_o16__f17)*Y(O16)*Y(H1)*state.rho +
        screened_rates(k_he4_o14__p_f17)*Y(He4)*Y(O14)*state.rho;

    ydot_nuc(F18) =
        -screened_rates(k_p_f18__he4_o15)*Y(F18)*Y(H1)*state.rho +
        screened_rates(k_he4_n14__f18)*Y(He4)*Y(N14)*state.rho +
        screened_rates(k_p_o17__f18)*Y(O17)*Y(H1)*state.rho;

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

    scratch = -screened_rates(k_p_be7__b8)*Y(Be7)*state.rho - screened_rates(k_p_c12__n13)*Y(C12)*state.rho - screened_rates(k_p_c13__n14)*Y(C13)*state.rho - screened_rates(k_p_d__he3)*Y(H2)*state.rho - screened_rates(k_p_f18__he4_o15)*Y(F18)*state.rho - screened_rates(k_p_he3__he4__weak__bet_pos_)*Y(He3)*state.rho - screened_rates(k_p_n13__o14)*Y(N13)*state.rho - screened_rates(k_p_n14__o15)*Y(N14)*state.rho - screened_rates(k_p_n15__he4_c12)*Y(N15)*state.rho - screened_rates(k_p_n15__o16)*Y(N15)*state.rho - screened_rates(k_p_o16__f17)*Y(O16)*state.rho - screened_rates(k_p_o17__f18)*Y(O17)*state.rho - screened_rates(k_p_o17__he4_n14)*Y(O17)*state.rho - 2.0*screened_rates(k_p_p__d__weak__bet_pos_)*Y(H1)*state.rho - 2.0*screened_rates(k_p_p__d__weak__electron_capture)*Y(H1)*std::pow(state.rho, 2)*state.y_e;
    jac.set(H1, H1, scratch);

    scratch = screened_rates(k_d_be7__p_he4_he4)*Y(Be7)*state.rho + screened_rates(k_d_he3__p_he4)*Y(He3)*state.rho - screened_rates(k_p_d__he3)*Y(H1)*state.rho;
    jac.set(H1, H2, scratch);

    scratch = screened_rates(k_d_he3__p_he4)*Y(H2)*state.rho + 2.0*screened_rates(k_he3_be7__p_p_he4_he4)*Y(Be7)*state.rho + 2.0*screened_rates(k_he3_he3__p_p_he4)*Y(He3)*state.rho - screened_rates(k_p_he3__he4__weak__bet_pos_)*Y(H1)*state.rho;
    jac.set(H1, He3, scratch);

    scratch = screened_rates(k_he4_n13__p_o16)*Y(N13)*state.rho + screened_rates(k_he4_o14__p_f17)*Y(O14)*state.rho;
    jac.set(H1, He4, scratch);

    scratch = screened_rates(k_d_be7__p_he4_he4)*Y(H2)*state.rho + 2.0*screened_rates(k_he3_be7__p_p_he4_he4)*Y(He3)*state.rho - screened_rates(k_p_be7__b8)*Y(H1)*state.rho;
    jac.set(H1, Be7, scratch);

    scratch = -screened_rates(k_p_c12__n13)*Y(H1)*state.rho;
    jac.set(H1, C12, scratch);

    scratch = -screened_rates(k_p_c13__n14)*Y(H1)*state.rho;
    jac.set(H1, C13, scratch);

    scratch = screened_rates(k_he4_n13__p_o16)*Y(He4)*state.rho - screened_rates(k_p_n13__o14)*Y(H1)*state.rho;
    jac.set(H1, N13, scratch);

    scratch = -screened_rates(k_p_n14__o15)*Y(H1)*state.rho;
    jac.set(H1, N14, scratch);

    scratch = -screened_rates(k_p_n15__he4_c12)*Y(H1)*state.rho - screened_rates(k_p_n15__o16)*Y(H1)*state.rho;
    jac.set(H1, N15, scratch);

    scratch = screened_rates(k_he4_o14__p_f17)*Y(He4)*state.rho;
    jac.set(H1, O14, scratch);

    scratch = -screened_rates(k_p_o16__f17)*Y(H1)*state.rho;
    jac.set(H1, O16, scratch);

    scratch = -screened_rates(k_p_o17__f18)*Y(H1)*state.rho - screened_rates(k_p_o17__he4_n14)*Y(H1)*state.rho;
    jac.set(H1, O17, scratch);

    scratch = -screened_rates(k_p_f18__he4_o15)*Y(H1)*state.rho;
    jac.set(H1, F18, scratch);

    scratch = -screened_rates(k_p_d__he3)*Y(H2)*state.rho + 1.0*screened_rates(k_p_p__d__weak__bet_pos_)*Y(H1)*state.rho + 1.0*screened_rates(k_p_p__d__weak__electron_capture)*Y(H1)*std::pow(state.rho, 2)*state.y_e;
    jac.set(H2, H1, scratch);

    scratch = -screened_rates(k_d_be7__p_he4_he4)*Y(Be7)*state.rho - 2.0*screened_rates(k_d_d__he4)*Y(H2)*state.rho - screened_rates(k_d_he3__p_he4)*Y(He3)*state.rho - screened_rates(k_p_d__he3)*Y(H1)*state.rho;
    jac.set(H2, H2, scratch);

    scratch = -screened_rates(k_d_he3__p_he4)*Y(H2)*state.rho;
    jac.set(H2, He3, scratch);

    scratch = -screened_rates(k_d_be7__p_he4_he4)*Y(H2)*state.rho;
    jac.set(H2, Be7, scratch);

    scratch = screened_rates(k_p_d__he3)*Y(H2)*state.rho - screened_rates(k_p_he3__he4__weak__bet_pos_)*Y(He3)*state.rho;
    jac.set(He3, H1, scratch);

    scratch = -screened_rates(k_d_he3__p_he4)*Y(He3)*state.rho + screened_rates(k_p_d__he3)*Y(H1)*state.rho;
    jac.set(He3, H2, scratch);

    scratch = -screened_rates(k_d_he3__p_he4)*Y(H2)*state.rho - screened_rates(k_he3_be7__p_p_he4_he4)*Y(Be7)*state.rho - 2.0*screened_rates(k_he3_he3__p_p_he4)*Y(He3)*state.rho - screened_rates(k_he4_he3__be7)*Y(He4)*state.rho - screened_rates(k_p_he3__he4__weak__bet_pos_)*Y(H1)*state.rho;
    jac.set(He3, He3, scratch);

    scratch = -screened_rates(k_he4_he3__be7)*Y(He3)*state.rho;
    jac.set(He3, He4, scratch);

    scratch = -screened_rates(k_he3_be7__p_p_he4_he4)*Y(He3)*state.rho;
    jac.set(He3, Be7, scratch);

    scratch = screened_rates(k_p_f18__he4_o15)*Y(F18)*state.rho + screened_rates(k_p_he3__he4__weak__bet_pos_)*Y(He3)*state.rho + screened_rates(k_p_n15__he4_c12)*Y(N15)*state.rho + screened_rates(k_p_o17__he4_n14)*Y(O17)*state.rho;
    jac.set(He4, H1, scratch);

    scratch = 2.0*screened_rates(k_d_be7__p_he4_he4)*Y(Be7)*state.rho + 1.0*screened_rates(k_d_d__he4)*Y(H2)*state.rho + screened_rates(k_d_he3__p_he4)*Y(He3)*state.rho;
    jac.set(He4, H2, scratch);

    scratch = screened_rates(k_d_he3__p_he4)*Y(H2)*state.rho + 2.0*screened_rates(k_he3_be7__p_p_he4_he4)*Y(Be7)*state.rho + 1.0*screened_rates(k_he3_he3__p_p_he4)*Y(He3)*state.rho - screened_rates(k_he4_he3__be7)*Y(He4)*state.rho + screened_rates(k_p_he3__he4__weak__bet_pos_)*Y(H1)*state.rho;
    jac.set(He4, He3, scratch);

    scratch = -screened_rates(k_he4_c12__o16)*Y(C12)*state.rho - screened_rates(k_he4_he3__be7)*Y(He3)*state.rho - 1.5*screened_rates(k_he4_he4_he4__c12)*std::pow(Y(He4), 2)*std::pow(state.rho, 2) - screened_rates(k_he4_n13__p_o16)*Y(N13)*state.rho - screened_rates(k_he4_n14__f18)*Y(N14)*state.rho - screened_rates(k_he4_o14__p_f17)*Y(O14)*state.rho;
    jac.set(He4, He4, scratch);

    scratch = 2.0*screened_rates(k_d_be7__p_he4_he4)*Y(H2)*state.rho + 2.0*screened_rates(k_he3_be7__p_p_he4_he4)*Y(He3)*state.rho;
    jac.set(He4, Be7, scratch);

    scratch = 2.0*screened_rates(k_b8__he4_he4__weak__wc12);
    jac.set(He4, B8, scratch);

    scratch = -screened_rates(k_he4_c12__o16)*Y(He4)*state.rho;
    jac.set(He4, C12, scratch);

    scratch = -screened_rates(k_he4_n13__p_o16)*Y(He4)*state.rho;
    jac.set(He4, N13, scratch);

    scratch = -screened_rates(k_he4_n14__f18)*Y(He4)*state.rho;
    jac.set(He4, N14, scratch);

    scratch = screened_rates(k_p_n15__he4_c12)*Y(H1)*state.rho;
    jac.set(He4, N15, scratch);

    scratch = -screened_rates(k_he4_o14__p_f17)*Y(He4)*state.rho;
    jac.set(He4, O14, scratch);

    scratch = screened_rates(k_p_o17__he4_n14)*Y(H1)*state.rho;
    jac.set(He4, O17, scratch);

    scratch = screened_rates(k_p_f18__he4_o15)*Y(H1)*state.rho;
    jac.set(He4, F18, scratch);

    scratch = -screened_rates(k_p_be7__b8)*Y(Be7)*state.rho;
    jac.set(Be7, H1, scratch);

    scratch = -screened_rates(k_d_be7__p_he4_he4)*Y(Be7)*state.rho;
    jac.set(Be7, H2, scratch);

    scratch = -screened_rates(k_he3_be7__p_p_he4_he4)*Y(Be7)*state.rho + screened_rates(k_he4_he3__be7)*Y(He4)*state.rho;
    jac.set(Be7, He3, scratch);

    scratch = screened_rates(k_he4_he3__be7)*Y(He3)*state.rho;
    jac.set(Be7, He4, scratch);

    scratch = -screened_rates(k_d_be7__p_he4_he4)*Y(H2)*state.rho - screened_rates(k_he3_be7__p_p_he4_he4)*Y(He3)*state.rho - screened_rates(k_p_be7__b8)*Y(H1)*state.rho;
    jac.set(Be7, Be7, scratch);

    scratch = screened_rates(k_p_be7__b8)*Y(Be7)*state.rho;
    jac.set(B8, H1, scratch);

    scratch = screened_rates(k_p_be7__b8)*Y(H1)*state.rho;
    jac.set(B8, Be7, scratch);

    scratch = -screened_rates(k_b8__he4_he4__weak__wc12);
    jac.set(B8, B8, scratch);

    scratch = -screened_rates(k_p_c12__n13)*Y(C12)*state.rho + screened_rates(k_p_n15__he4_c12)*Y(N15)*state.rho;
    jac.set(C12, H1, scratch);

    scratch = -screened_rates(k_he4_c12__o16)*Y(C12)*state.rho + 0.5*screened_rates(k_he4_he4_he4__c12)*std::pow(Y(He4), 2)*std::pow(state.rho, 2);
    jac.set(C12, He4, scratch);

    scratch = -screened_rates(k_he4_c12__o16)*Y(He4)*state.rho - screened_rates(k_p_c12__n13)*Y(H1)*state.rho;
    jac.set(C12, C12, scratch);

    scratch = screened_rates(k_p_n15__he4_c12)*Y(H1)*state.rho;
    jac.set(C12, N15, scratch);

    scratch = -screened_rates(k_p_c13__n14)*Y(C13)*state.rho;
    jac.set(C13, H1, scratch);

    scratch = -screened_rates(k_p_c13__n14)*Y(H1)*state.rho;
    jac.set(C13, C13, scratch);

    scratch = screened_rates(k_n13__c13__weak__wc12);
    jac.set(C13, N13, scratch);

    scratch = screened_rates(k_p_c12__n13)*Y(C12)*state.rho - screened_rates(k_p_n13__o14)*Y(N13)*state.rho;
    jac.set(N13, H1, scratch);

    scratch = -screened_rates(k_he4_n13__p_o16)*Y(N13)*state.rho;
    jac.set(N13, He4, scratch);

    scratch = screened_rates(k_p_c12__n13)*Y(H1)*state.rho;
    jac.set(N13, C12, scratch);

    scratch = -screened_rates(k_he4_n13__p_o16)*Y(He4)*state.rho - screened_rates(k_n13__c13__weak__wc12) - screened_rates(k_p_n13__o14)*Y(H1)*state.rho;
    jac.set(N13, N13, scratch);

    scratch = screened_rates(k_p_c13__n14)*Y(C13)*state.rho - screened_rates(k_p_n14__o15)*Y(N14)*state.rho + screened_rates(k_p_o17__he4_n14)*Y(O17)*state.rho;
    jac.set(N14, H1, scratch);

    scratch = -screened_rates(k_he4_n14__f18)*Y(N14)*state.rho;
    jac.set(N14, He4, scratch);

    scratch = screened_rates(k_p_c13__n14)*Y(H1)*state.rho;
    jac.set(N14, C13, scratch);

    scratch = -screened_rates(k_he4_n14__f18)*Y(He4)*state.rho - screened_rates(k_p_n14__o15)*Y(H1)*state.rho;
    jac.set(N14, N14, scratch);

    scratch = screened_rates(k_o14__n14__weak__wc12);
    jac.set(N14, O14, scratch);

    scratch = screened_rates(k_p_o17__he4_n14)*Y(H1)*state.rho;
    jac.set(N14, O17, scratch);

    scratch = -screened_rates(k_p_n15__he4_c12)*Y(N15)*state.rho - screened_rates(k_p_n15__o16)*Y(N15)*state.rho;
    jac.set(N15, H1, scratch);

    scratch = -screened_rates(k_p_n15__he4_c12)*Y(H1)*state.rho - screened_rates(k_p_n15__o16)*Y(H1)*state.rho;
    jac.set(N15, N15, scratch);

    scratch = screened_rates(k_o15__n15__weak__wc12);
    jac.set(N15, O15, scratch);

    scratch = screened_rates(k_p_n13__o14)*Y(N13)*state.rho;
    jac.set(O14, H1, scratch);

    scratch = -screened_rates(k_he4_o14__p_f17)*Y(O14)*state.rho;
    jac.set(O14, He4, scratch);

    scratch = screened_rates(k_p_n13__o14)*Y(H1)*state.rho;
    jac.set(O14, N13, scratch);

    scratch = -screened_rates(k_he4_o14__p_f17)*Y(He4)*state.rho - screened_rates(k_o14__n14__weak__wc12);
    jac.set(O14, O14, scratch);

    scratch = screened_rates(k_p_f18__he4_o15)*Y(F18)*state.rho + screened_rates(k_p_n14__o15)*Y(N14)*state.rho;
    jac.set(O15, H1, scratch);

    scratch = screened_rates(k_p_n14__o15)*Y(H1)*state.rho;
    jac.set(O15, N14, scratch);

    scratch = -screened_rates(k_o15__n15__weak__wc12);
    jac.set(O15, O15, scratch);

    scratch = screened_rates(k_p_f18__he4_o15)*Y(H1)*state.rho;
    jac.set(O15, F18, scratch);

    scratch = screened_rates(k_p_n15__o16)*Y(N15)*state.rho - screened_rates(k_p_o16__f17)*Y(O16)*state.rho;
    jac.set(O16, H1, scratch);

    scratch = screened_rates(k_he4_c12__o16)*Y(C12)*state.rho + screened_rates(k_he4_n13__p_o16)*Y(N13)*state.rho;
    jac.set(O16, He4, scratch);

    scratch = screened_rates(k_he4_c12__o16)*Y(He4)*state.rho;
    jac.set(O16, C12, scratch);

    scratch = screened_rates(k_he4_n13__p_o16)*Y(He4)*state.rho;
    jac.set(O16, N13, scratch);

    scratch = screened_rates(k_p_n15__o16)*Y(H1)*state.rho;
    jac.set(O16, N15, scratch);

    scratch = -screened_rates(k_p_o16__f17)*Y(H1)*state.rho;
    jac.set(O16, O16, scratch);

    scratch = -screened_rates(k_p_o17__f18)*Y(O17)*state.rho - screened_rates(k_p_o17__he4_n14)*Y(O17)*state.rho;
    jac.set(O17, H1, scratch);

    scratch = -screened_rates(k_p_o17__f18)*Y(H1)*state.rho - screened_rates(k_p_o17__he4_n14)*Y(H1)*state.rho;
    jac.set(O17, O17, scratch);

    scratch = screened_rates(k_f17__o17__weak__wc12);
    jac.set(O17, F17, scratch);

    scratch = screened_rates(k_p_o16__f17)*Y(O16)*state.rho;
    jac.set(F17, H1, scratch);

    scratch = screened_rates(k_he4_o14__p_f17)*Y(O14)*state.rho;
    jac.set(F17, He4, scratch);

    scratch = screened_rates(k_he4_o14__p_f17)*Y(He4)*state.rho;
    jac.set(F17, O14, scratch);

    scratch = screened_rates(k_p_o16__f17)*Y(H1)*state.rho;
    jac.set(F17, O16, scratch);

    scratch = -screened_rates(k_f17__o17__weak__wc12);
    jac.set(F17, F17, scratch);

    scratch = -screened_rates(k_p_f18__he4_o15)*Y(F18)*state.rho + screened_rates(k_p_o17__f18)*Y(O17)*state.rho;
    jac.set(F18, H1, scratch);

    scratch = screened_rates(k_he4_n14__f18)*Y(N14)*state.rho;
    jac.set(F18, He4, scratch);

    scratch = screened_rates(k_he4_n14__f18)*Y(He4)*state.rho;
    jac.set(F18, N14, scratch);

    scratch = screened_rates(k_p_o17__f18)*Y(H1)*state.rho;
    jac.set(F18, O17, scratch);

    scratch = -screened_rates(k_p_f18__he4_o15)*Y(H1)*state.rho;
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

    screening_init();

    net_screening_init();

}


#endif
