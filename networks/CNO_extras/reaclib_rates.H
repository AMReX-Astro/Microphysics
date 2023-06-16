#ifndef REACLIB_RATES_H
#define REACLIB_RATES_H

#include <AMReX.H>
#include <AMReX_Print.H>

#include <tfactors.H>
#include <actual_network.H>
#include <partition_functions.H>

using namespace Rates;
using namespace Species;

struct rate_t {
    Array1D<Real, 1, NumRates>  screened_rates;
};

struct rate_derivs_t {
    Array1D<Real, 1, NumRates>  screened_rates;
    Array1D<Real, 1, NumRates>  dscreened_rates_dT;
};


template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_n13__c13__weak__wc12(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // n13 --> c13

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // wc12w
    ln_set_rate =  -6.7601;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 = 0.0;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_o14__n14__weak__wc12(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // o14 --> n14

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // wc12w
    ln_set_rate =  -4.62354;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 = 0.0;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_o15__n15__weak__wc12(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // o15 --> n15

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // wc12w
    ln_set_rate =  -5.17053;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 = 0.0;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_f17__o17__weak__wc12(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // f17 --> o17

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // wc12w
    ln_set_rate =  -4.53318;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 = 0.0;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_f18__o18__weak__wc12(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // f18 --> o18

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // wc12w
    ln_set_rate =  -9.15982;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 = 0.0;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_ne18__f18__weak__wc12(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ne18 --> f18

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // wc12w
    ln_set_rate =  -0.879336;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 = 0.0;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_ne19__f19__weak__wc12(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ne19 --> f19

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // wc12w
    ln_set_rate =  -3.21142;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 = 0.0;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_c12__n13(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // c12 + p --> n13

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ls09n
    ln_set_rate =  17.1482 + -13.692 * tfactors.T913i + -0.230881 * tfactors.T913
                         + 4.44362 * tfactors.T9 + -3.15898 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -13.692 * tfactors.T943i + (1.0/3.0) * -0.230881 * tfactors.T923i
                                  + 4.44362 + (5.0/3.0) * -3.15898 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // ls09r
    ln_set_rate =  17.5428 + -3.77849 * tfactors.T9i + -5.10735 * tfactors.T913i + -2.24111 * tfactors.T913
                         + 0.148883 * tfactors.T9 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  3.77849 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -5.10735 * tfactors.T943i + (1.0/3.0) * -2.24111 * tfactors.T923i
                                  + 0.148883 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_c12__o16(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // c12 + he4 --> o16

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // nac2 
    ln_set_rate =  254.634 + -1.84097 * tfactors.T9i + 103.411 * tfactors.T913i + -420.567 * tfactors.T913
                         + 64.0874 * tfactors.T9 + -12.4624 * tfactors.T953 + 137.303 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.84097 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 103.411 * tfactors.T943i + (1.0/3.0) * -420.567 * tfactors.T923i
                                  + 64.0874 + (5.0/3.0) * -12.4624 * tfactors.T923 + 137.303 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // nac2 
    ln_set_rate =  69.6526 + -1.39254 * tfactors.T9i + 58.9128 * tfactors.T913i + -148.273 * tfactors.T913
                         + 9.08324 * tfactors.T9 + -0.541041 * tfactors.T953 + 70.3554 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.39254 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 58.9128 * tfactors.T943i + (1.0/3.0) * -148.273 * tfactors.T923i
                                  + 9.08324 + (5.0/3.0) * -0.541041 * tfactors.T923 + 70.3554 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_c13__n14(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // c13 + p --> n14

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // nacrr
    ln_set_rate =  15.1825 + -13.5543 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  13.5543 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // nacrn
    ln_set_rate =  18.5155 + -13.72 * tfactors.T913i + -0.450018 * tfactors.T913
                         + 3.70823 * tfactors.T9 + -1.70545 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -13.72 * tfactors.T943i + (1.0/3.0) * -0.450018 * tfactors.T923i
                                  + 3.70823 + (5.0/3.0) * -1.70545 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // nacrr
    ln_set_rate =  13.9637 + -5.78147 * tfactors.T9i + -0.196703 * tfactors.T913
                         + 0.142126 * tfactors.T9 + -0.0238912 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  5.78147 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.196703 * tfactors.T923i
                                  + 0.142126 + (5.0/3.0) * -0.0238912 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_n13__o14(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // n13 + p --> o14

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // lg06r
    ln_set_rate =  10.9971 + -6.12602 * tfactors.T9i + 1.57122 * tfactors.T913i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  6.12602 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 1.57122 * tfactors.T943i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // lg06n
    ln_set_rate =  18.1356 + -15.1676 * tfactors.T913i + 0.0955166 * tfactors.T913
                         + 3.0659 * tfactors.T9 + -0.507339 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -15.1676 * tfactors.T943i + (1.0/3.0) * 0.0955166 * tfactors.T923i
                                  + 3.0659 + (5.0/3.0) * -0.507339 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_n14__o15(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // n14 + p --> o15

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // im05n
    ln_set_rate =  17.01 + -15.193 * tfactors.T913i + -0.161954 * tfactors.T913
                         + -7.52123 * tfactors.T9 + -0.987565 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -15.193 * tfactors.T943i + (1.0/3.0) * -0.161954 * tfactors.T923i
                                  + -7.52123 + (5.0/3.0) * -0.987565 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // im05r
    ln_set_rate =  6.73578 + -4.891 * tfactors.T9i
                         + 0.0682 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  4.891 * tfactors.T9i * tfactors.T9i
                                  + 0.0682 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // im05r
    ln_set_rate =  7.65444 + -2.998 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  2.998 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // im05n
    ln_set_rate =  20.1169 + -15.193 * tfactors.T913i + -4.63975 * tfactors.T913
                         + 9.73458 * tfactors.T9 + -9.55051 * tfactors.T953 + 0.333333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -15.193 * tfactors.T943i + (1.0/3.0) * -4.63975 * tfactors.T923i
                                  + 9.73458 + (5.0/3.0) * -9.55051 * tfactors.T923 + 0.333333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_n14__f18(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // n14 + he4 --> f18

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10n
    ln_set_rate =  21.5339 + -36.2504 * tfactors.T913i
                         + -5.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -36.2504 * tfactors.T943i
                                  + (5.0/3.0) * -5.0 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  13.8995 + -10.9656 * tfactors.T9i + -5.6227 * tfactors.T913i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  10.9656 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -5.6227 * tfactors.T943i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  0.196838 + -5.16034 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  5.16034 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_n15__o16(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // n15 + p --> o16

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // li10n
    ln_set_rate =  20.0176 + -15.24 * tfactors.T913i + 0.334926 * tfactors.T913
                         + 4.59088 * tfactors.T9 + -4.78468 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -15.24 * tfactors.T943i + (1.0/3.0) * 0.334926 * tfactors.T923i
                                  + 4.59088 + (5.0/3.0) * -4.78468 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // li10r
    ln_set_rate =  14.5444 + -10.2295 * tfactors.T9i
                         + 0.0459037 * tfactors.T9 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  10.2295 * tfactors.T9i * tfactors.T9i
                                  + 0.0459037 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // li10r
    ln_set_rate =  6.59056 + -2.92315 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  2.92315 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_n15__f19(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // n15 + he4 --> f19

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10r
    ln_set_rate =  -9.41892 + -4.17795 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  4.17795 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10n
    ln_set_rate =  25.3916 + -36.2324 * tfactors.T913i
                         + -2.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -36.2324 * tfactors.T943i
                                  + (5.0/3.0) * -2.0 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  -28.7989 + -4.19986 * tfactors.T9i + 35.4292 * tfactors.T913
                         + -5.5767 * tfactors.T9 + 0.441293 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  4.19986 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 35.4292 * tfactors.T923i
                                  + -5.5767 + (5.0/3.0) * 0.441293 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  3.5342 + -6.98462 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  6.98462 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_o14__ne18(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // o14 + he4 --> ne18

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // wh87r
    ln_set_rate =  -4.69948 + -12.159 * tfactors.T9i
                         + 5.0 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  12.159 * tfactors.T9i * tfactors.T9i
                                  + 5.0 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // wh87r
    ln_set_rate =  3.52636 + -22.61 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  22.61 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // wh87r
    ln_set_rate =  -2.15417 + -11.73 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  11.73 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // wh87n
    ln_set_rate =  26.4429 + -39.38 * tfactors.T913i + -0.0772187 * tfactors.T913
                         + -0.635361 * tfactors.T9 + 0.106236 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -39.38 * tfactors.T943i + (1.0/3.0) * -0.0772187 * tfactors.T923i
                                  + -0.635361 + (5.0/3.0) * 0.106236 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_o15__ne19(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // o15 + he4 --> ne19

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // dc11r
    ln_set_rate =  -32.2496 + -4.20439 * tfactors.T9i + -3.24609 * tfactors.T913i + 44.4647 * tfactors.T913
                         + -9.79962 * tfactors.T9 + 0.841782 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  4.20439 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -3.24609 * tfactors.T943i + (1.0/3.0) * 44.4647 * tfactors.T923i
                                  + -9.79962 + (5.0/3.0) * 0.841782 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // dc11r
    ln_set_rate =  -0.0452465 + -5.88439 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  5.88439 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // dc11n
    ln_set_rate =  26.2914 + -39.578 * tfactors.T913i
                         + -3.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -39.578 * tfactors.T943i
                                  + (5.0/3.0) * -3.0 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_o16__f17(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // o16 + p --> f17

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ia08n
    ln_set_rate =  19.0904 + -16.696 * tfactors.T913i + -1.16252 * tfactors.T913
                         + 0.267703 * tfactors.T9 + -0.0338411 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -16.696 * tfactors.T943i + (1.0/3.0) * -1.16252 * tfactors.T923i
                                  + 0.267703 + (5.0/3.0) * -0.0338411 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_o16__ne20(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // o16 + he4 --> ne20

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // co10r
    ln_set_rate =  9.50848 + -12.7643 * tfactors.T9i + -3.65925 * tfactors.T913
                         + 0.714224 * tfactors.T9 + -0.00107508 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  12.7643 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -3.65925 * tfactors.T923i
                                  + 0.714224 + (5.0/3.0) * -0.00107508 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // co10r
    ln_set_rate =  3.88571 + -10.3585 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  10.3585 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // co10n
    ln_set_rate =  23.903 + -39.7262 * tfactors.T913i + -0.210799 * tfactors.T913
                         + 0.442879 * tfactors.T9 + -0.0797753 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -39.7262 * tfactors.T943i + (1.0/3.0) * -0.210799 * tfactors.T923i
                                  + 0.442879 + (5.0/3.0) * -0.0797753 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_o17__f18(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // o17 + p --> f18

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10n
    ln_set_rate =  15.8929 + -16.4035 * tfactors.T913i + 4.31885 * tfactors.T913
                         + -0.709921 * tfactors.T9 + -2.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -16.4035 * tfactors.T943i + (1.0/3.0) * 4.31885 * tfactors.T923i
                                  + -0.709921 + (5.0/3.0) * -2.0 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  9.39048 + -6.22828 * tfactors.T9i + 2.31435 * tfactors.T913
                         + -0.302835 * tfactors.T9 + 0.020133 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  6.22828 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 2.31435 * tfactors.T923i
                                  + -0.302835 + (5.0/3.0) * 0.020133 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  -13.077 + -0.746296 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.746296 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_o18__f19(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // o18 + p --> f19

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10r
    ln_set_rate =  -35.0079 + -0.244743 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.244743 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10n
    ln_set_rate =  19.917 + -16.7246 * tfactors.T913i
                         + -3.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -16.7246 * tfactors.T943i
                                  + (5.0/3.0) * -3.0 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  7.26876 + -6.7253 * tfactors.T9i + 3.99059 * tfactors.T913
                         + -0.593127 * tfactors.T9 + 0.0877534 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  6.7253 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 3.99059 * tfactors.T923i
                                  + -0.593127 + (5.0/3.0) * 0.0877534 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  5.07648 + -1.65681 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.65681 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_f17__ne18(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // f17 + p --> ne18

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // cb09 
    ln_set_rate =  -7.84708 + -0.0323504 * tfactors.T9i + -14.2191 * tfactors.T913i + 34.0647 * tfactors.T913
                         + -16.5698 * tfactors.T9 + 2.48116 * tfactors.T953 + -2.13376 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.0323504 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -14.2191 * tfactors.T943i + (1.0/3.0) * 34.0647 * tfactors.T923i
                                  + -16.5698 + (5.0/3.0) * 2.48116 * tfactors.T923 + -2.13376 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // cb09 
    ln_set_rate =  27.5778 + -4.95969 * tfactors.T9i + -21.3249 * tfactors.T913i + -0.230774 * tfactors.T913
                         + 0.917931 * tfactors.T9 + -0.0440377 * tfactors.T953 + -7.36014 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  4.95969 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -21.3249 * tfactors.T943i + (1.0/3.0) * -0.230774 * tfactors.T923i
                                  + 0.917931 + (5.0/3.0) * -0.0440377 * tfactors.T923 + -7.36014 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_f18__ne19(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // f18 + p --> ne19

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10n
    ln_set_rate =  57.4084 + -21.4023 * tfactors.T913i + -93.766 * tfactors.T913
                         + 179.258 * tfactors.T9 + -202.561 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -21.4023 * tfactors.T943i + (1.0/3.0) * -93.766 * tfactors.T923i
                                  + 179.258 + (5.0/3.0) * -202.561 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  -5.85727 + -2.89147 * tfactors.T9i + 13.1683 * tfactors.T913
                         + -1.92023 * tfactors.T9 + 0.16901 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  2.89147 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 13.1683 * tfactors.T923i
                                  + -1.92023 + (5.0/3.0) * 0.16901 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  -29.449 + -0.39895 * tfactors.T9i + 22.4903 * tfactors.T913
                         + 0.307872 * tfactors.T9 + -0.296226 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.39895 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 22.4903 * tfactors.T923i
                                  + 0.307872 + (5.0/3.0) * -0.296226 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_f19__ne20(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // f19 + p --> ne20

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // nacrr
    ln_set_rate =  -5.63093 + -7.74414 * tfactors.T9i + 31.6442 * tfactors.T913i + -58.6563 * tfactors.T913
                         + 67.7365 * tfactors.T9 + -22.9721 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  7.74414 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 31.6442 * tfactors.T943i + (1.0/3.0) * -58.6563 * tfactors.T923i
                                  + 67.7365 + (5.0/3.0) * -22.9721 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // nacrr
    ln_set_rate =  12.3816 + -1.71383 * tfactors.T9i + -11.3832 * tfactors.T913i + 5.47872 * tfactors.T913
                         + -1.07203 * tfactors.T9 + 0.11196 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.71383 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -11.3832 * tfactors.T943i + (1.0/3.0) * 5.47872 * tfactors.T923i
                                  + -1.07203 + (5.0/3.0) * 0.11196 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // nacrn
    ln_set_rate =  18.2807 + -18.116 * tfactors.T913i + -1.4622 * tfactors.T913
                         + 6.95113 * tfactors.T9 + -2.90366 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -18.116 * tfactors.T943i + (1.0/3.0) * -1.4622 * tfactors.T923i
                                  + 6.95113 + (5.0/3.0) * -2.90366 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_ne18__mg22(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ne18 + he4 --> mg22

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  32.8865 + -46.4859 * tfactors.T913i + 0.956741 * tfactors.T913
                         + -0.914402 * tfactors.T9 + 0.0722478 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -46.4859 * tfactors.T943i + (1.0/3.0) * 0.956741 * tfactors.T923i
                                  + -0.914402 + (5.0/3.0) * 0.0722478 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_ne20__mg24(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ne20 + he4 --> mg24

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10r
    ln_set_rate =  -38.7055 + -2.50605 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  2.50605 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10n
    ln_set_rate =  24.5058 + -46.2525 * tfactors.T913i + 5.58901 * tfactors.T913
                         + 7.61843 * tfactors.T9 + -3.683 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -46.2525 * tfactors.T943i + (1.0/3.0) * 5.58901 * tfactors.T923i
                                  + 7.61843 + (5.0/3.0) * -3.683 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  -8.79827 + -12.7809 * tfactors.T9i + 16.9229 * tfactors.T913
                         + -2.57325 * tfactors.T9 + 0.208997 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  12.7809 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 16.9229 * tfactors.T923i
                                  + -2.57325 + (5.0/3.0) * 0.208997 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  1.98307 + -9.22026 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  9.22026 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_c12_c12__he4_ne20(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // c12 + c12 --> he4 + ne20

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // cf88r
    ln_set_rate =  61.2863 + -84.165 * tfactors.T913i + -1.56627 * tfactors.T913
                         + -0.0736084 * tfactors.T9 + -0.072797 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -84.165 * tfactors.T943i + (1.0/3.0) * -1.56627 * tfactors.T923i
                                  + -0.0736084 + (5.0/3.0) * -0.072797 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_n13__p_o16(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // n13 + he4 --> p + o16

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // cf88n
    ln_set_rate =  40.4644 + -35.829 * tfactors.T913i + -0.530275 * tfactors.T913
                         + -0.982462 * tfactors.T9 + 0.0808059 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -35.829 * tfactors.T943i + (1.0/3.0) * -0.530275 * tfactors.T923i
                                  + -0.982462 + (5.0/3.0) * 0.0808059 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_n15__he4_c12(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // n15 + p --> he4 + c12

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // nacrn
    ln_set_rate =  27.4764 + -15.253 * tfactors.T913i + 1.59318 * tfactors.T913
                         + 2.4479 * tfactors.T9 + -2.19708 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -15.253 * tfactors.T943i + (1.0/3.0) * 1.59318 * tfactors.T923i
                                  + 2.4479 + (5.0/3.0) * -2.19708 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // nacrr
    ln_set_rate =  -6.57522 + -1.1638 * tfactors.T9i + 22.7105 * tfactors.T913
                         + -2.90707 * tfactors.T9 + 0.205754 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.1638 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 22.7105 * tfactors.T923i
                                  + -2.90707 + (5.0/3.0) * 0.205754 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // nacrr
    ln_set_rate =  20.8972 + -7.406 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  7.406 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // nacrr
    ln_set_rate =  -4.87347 + -2.02117 * tfactors.T9i + 30.8497 * tfactors.T913
                         + -8.50433 * tfactors.T9 + -1.54426 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  2.02117 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 30.8497 * tfactors.T923i
                                  + -8.50433 + (5.0/3.0) * -1.54426 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_o14__p_f17(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // o14 + he4 --> p + f17

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // Ha96n
    ln_set_rate =  40.8358 + -39.388 * tfactors.T913i + -17.4673 * tfactors.T913
                         + 35.3029 * tfactors.T9 + -24.8162 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -39.388 * tfactors.T943i + (1.0/3.0) * -17.4673 * tfactors.T923i
                                  + 35.3029 + (5.0/3.0) * -24.8162 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // Ha96r
    ln_set_rate =  16.3087 + -22.51 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  22.51 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // Ha96r
    ln_set_rate =  11.1184 + -13.6 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  13.6 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // Ha96r
    ln_set_rate =  -106.091 + -0.453036 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.453036 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // Ha96r
    ln_set_rate =  12.1289 + -12.0223 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  12.0223 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // Ha96r
    ln_set_rate =  18.6518 + -26.0 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  26.0 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_c12_o16__he4_mg24(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // o16 + c12 --> he4 + mg24

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // cf88r
    ln_set_rate =  48.5341 + 0.37204 * tfactors.T9i + -133.413 * tfactors.T913i + 50.1572 * tfactors.T913
                         + -3.15987 * tfactors.T9 + 0.0178251 * tfactors.T953 + -23.7027 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  -0.37204 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -133.413 * tfactors.T943i + (1.0/3.0) * 50.1572 * tfactors.T923i
                                  + -3.15987 + (5.0/3.0) * 0.0178251 * tfactors.T923 + -23.7027 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_o17__he4_n14(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // o17 + p --> he4 + n14

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10r
    ln_set_rate =  5.5336 + -2.11477 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  2.11477 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  -7.20763 + -0.753395 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.753395 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10n
    ln_set_rate =  19.579 + -16.9078 * tfactors.T913i
                         + -2.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -16.9078 * tfactors.T943i
                                  + (5.0/3.0) * -2.0 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  10.174 + -4.95865 * tfactors.T9i + 5.10182 * tfactors.T913
                         + 0.379373 * tfactors.T9 + -0.0672515 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  4.95865 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 5.10182 * tfactors.T923i
                                  + 0.379373 + (5.0/3.0) * -0.0672515 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_o18__he4_n15(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // o18 + p --> he4 + n15

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10r
    ln_set_rate =  10.2725 + -1.663 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.663 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  -27.9044 + -0.245884 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.245884 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10n
    ln_set_rate =  26.9671 + -16.6979 * tfactors.T913i
                         + -3.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -16.6979 * tfactors.T943i
                                  + (5.0/3.0) * -3.0 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  8.94352 + -5.32335 * tfactors.T9i + 11.6568 * tfactors.T913
                         + -2.16303 * tfactors.T9 + 0.209965 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  5.32335 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 11.6568 * tfactors.T923i
                                  + -2.16303 + (5.0/3.0) * 0.209965 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_f18__he4_o15(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // f18 + p --> he4 + o15

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10n
    ln_set_rate =  62.0058 + -21.4023 * tfactors.T913i + -80.8891 * tfactors.T913
                         + 134.6 * tfactors.T9 + -126.504 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -21.4023 * tfactors.T943i + (1.0/3.0) * -80.8891 * tfactors.T923i
                                  + 134.6 + (5.0/3.0) * -126.504 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  1.75704 + -3.01675 * tfactors.T9i + 13.3223 * tfactors.T913
                         + -1.36696 * tfactors.T9 + 0.0757363 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  3.01675 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 13.3223 * tfactors.T923i
                                  + -1.36696 + (5.0/3.0) * 0.0757363 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  -31.7388 + -0.376432 * tfactors.T9i + 61.738 * tfactors.T913
                         + -108.29 * tfactors.T9 + -34.2365 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.376432 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 61.738 * tfactors.T923i
                                  + -108.29 + (5.0/3.0) * -34.2365 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_f19__he4_o16(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // f19 + p --> he4 + o16

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // nacr 
    ln_set_rate =  8.239 + -2.46828 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  2.46828 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // nacr 
    ln_set_rate =  -52.7043 + -0.12765 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.12765 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // nacr 
    ln_set_rate =  26.2916 + -18.116 * tfactors.T913i
                         + 1.86674 * tfactors.T9 + -7.5666 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -18.116 * tfactors.T943i
                                  + 1.86674 + (5.0/3.0) * -7.5666 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // nacrr
    ln_set_rate =  14.3586 + -3.286 * tfactors.T9i
                         + -0.21103 * tfactors.T9 + 2.87702 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  3.286 * tfactors.T9i * tfactors.T9i
                                  + -0.21103 + 2.87702 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // nacr 
    ln_set_rate =  15.1955 + -3.75185 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  3.75185 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_ne20__he4_f17(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ne20 + p --> he4 + f17

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // nacr 
    ln_set_rate =  41.563 + -47.9266 * tfactors.T9i + -43.18 * tfactors.T913i + 4.46827 * tfactors.T913
                         + -1.63915 * tfactors.T9 + 0.123483 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  47.9266 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -43.18 * tfactors.T943i + (1.0/3.0) * 4.46827 * tfactors.T923i
                                  + -1.63915 + (5.0/3.0) * 0.123483 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_he4_he4__c12(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // he4 + he4 + he4 --> c12

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // fy05r
    ln_set_rate =  -24.3505 + -4.12656 * tfactors.T9i + -13.49 * tfactors.T913i + 21.4259 * tfactors.T913
                         + -1.34769 * tfactors.T9 + 0.0879816 * tfactors.T953 + -13.1653 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  4.12656 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -13.49 * tfactors.T943i + (1.0/3.0) * 21.4259 * tfactors.T923i
                                  + -1.34769 + (5.0/3.0) * 0.0879816 * tfactors.T923 + -13.1653 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // fy05r
    ln_set_rate =  -11.7884 + -1.02446 * tfactors.T9i + -23.57 * tfactors.T913i + 20.4886 * tfactors.T913
                         + -12.9882 * tfactors.T9 + -20.0 * tfactors.T953 + -2.16667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.02446 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -23.57 * tfactors.T943i + (1.0/3.0) * 20.4886 * tfactors.T923i
                                  + -12.9882 + (5.0/3.0) * -20.0 * tfactors.T923 + -2.16667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // fy05n
    ln_set_rate =  -0.971052 + -37.06 * tfactors.T913i + 29.3493 * tfactors.T913
                         + -115.507 * tfactors.T9 + -10.0 * tfactors.T953 + -1.33333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -37.06 * tfactors.T943i + (1.0/3.0) * 29.3493 * tfactors.T923i
                                  + -115.507 + (5.0/3.0) * -10.0 * tfactors.T923 + -1.33333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}



template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
fill_reaclib_rates(const tf_t& tfactors, T& rate_eval)
{

    Real rate;
    Real drate_dT;

    rate_n13__c13__weak__wc12<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n13__c13__weak__wc12) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_n13__c13__weak__wc12) = drate_dT;

    }
    rate_o14__n14__weak__wc12<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_o14__n14__weak__wc12) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_o14__n14__weak__wc12) = drate_dT;

    }
    rate_o15__n15__weak__wc12<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_o15__n15__weak__wc12) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_o15__n15__weak__wc12) = drate_dT;

    }
    rate_f17__o17__weak__wc12<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_f17__o17__weak__wc12) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_f17__o17__weak__wc12) = drate_dT;

    }
    rate_f18__o18__weak__wc12<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_f18__o18__weak__wc12) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_f18__o18__weak__wc12) = drate_dT;

    }
    rate_ne18__f18__weak__wc12<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_ne18__f18__weak__wc12) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_ne18__f18__weak__wc12) = drate_dT;

    }
    rate_ne19__f19__weak__wc12<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_ne19__f19__weak__wc12) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_ne19__f19__weak__wc12) = drate_dT;

    }
    rate_p_c12__n13<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_c12__n13) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_c12__n13) = drate_dT;

    }
    rate_he4_c12__o16<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_c12__o16) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_c12__o16) = drate_dT;

    }
    rate_p_c13__n14<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_c13__n14) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_c13__n14) = drate_dT;

    }
    rate_p_n13__o14<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_n13__o14) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_n13__o14) = drate_dT;

    }
    rate_p_n14__o15<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_n14__o15) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_n14__o15) = drate_dT;

    }
    rate_he4_n14__f18<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_n14__f18) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_n14__f18) = drate_dT;

    }
    rate_p_n15__o16<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_n15__o16) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_n15__o16) = drate_dT;

    }
    rate_he4_n15__f19<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_n15__f19) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_n15__f19) = drate_dT;

    }
    rate_he4_o14__ne18<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_o14__ne18) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_o14__ne18) = drate_dT;

    }
    rate_he4_o15__ne19<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_o15__ne19) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_o15__ne19) = drate_dT;

    }
    rate_p_o16__f17<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_o16__f17) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_o16__f17) = drate_dT;

    }
    rate_he4_o16__ne20<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_o16__ne20) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_o16__ne20) = drate_dT;

    }
    rate_p_o17__f18<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_o17__f18) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_o17__f18) = drate_dT;

    }
    rate_p_o18__f19<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_o18__f19) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_o18__f19) = drate_dT;

    }
    rate_p_f17__ne18<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_f17__ne18) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_f17__ne18) = drate_dT;

    }
    rate_p_f18__ne19<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_f18__ne19) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_f18__ne19) = drate_dT;

    }
    rate_p_f19__ne20<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_f19__ne20) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_f19__ne20) = drate_dT;

    }
    rate_he4_ne18__mg22<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_ne18__mg22) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_ne18__mg22) = drate_dT;

    }
    rate_he4_ne20__mg24<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_ne20__mg24) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_ne20__mg24) = drate_dT;

    }
    rate_c12_c12__he4_ne20<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_c12_c12__he4_ne20) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_c12_c12__he4_ne20) = drate_dT;

    }
    rate_he4_n13__p_o16<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_n13__p_o16) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_n13__p_o16) = drate_dT;

    }
    rate_p_n15__he4_c12<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_n15__he4_c12) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_n15__he4_c12) = drate_dT;

    }
    rate_he4_o14__p_f17<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_o14__p_f17) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_o14__p_f17) = drate_dT;

    }
    rate_c12_o16__he4_mg24<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_c12_o16__he4_mg24) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_c12_o16__he4_mg24) = drate_dT;

    }
    rate_p_o17__he4_n14<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_o17__he4_n14) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_o17__he4_n14) = drate_dT;

    }
    rate_p_o18__he4_n15<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_o18__he4_n15) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_o18__he4_n15) = drate_dT;

    }
    rate_p_f18__he4_o15<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_f18__he4_o15) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_f18__he4_o15) = drate_dT;

    }
    rate_p_f19__he4_o16<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_f19__he4_o16) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_f19__he4_o16) = drate_dT;

    }
    rate_p_ne20__he4_f17<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_ne20__he4_f17) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_ne20__he4_f17) = drate_dT;

    }
    rate_he4_he4_he4__c12<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_he4_he4__c12) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_he4_he4__c12) = drate_dT;

    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
fill_approx_rates([[maybe_unused]] const tf_t& tfactors, T& rate_eval)
{

    Real rate;
    Real drate_dT;


}

#endif
