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
    amrex::Array1D<amrex::Real, 1, NumRates>  screened_rates;
    amrex::Real enuc_weak;
};

struct rate_derivs_t {
    amrex::Array1D<amrex::Real, 1, NumRates>  screened_rates;
    amrex::Array1D<amrex::Real, 1, NumRates>  dscreened_rates_dT;
    amrex::Real enuc_weak;
};


template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_N13_to_C13_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // N13 --> C13

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // wc12w
    ln_set_rate =  -6.7601;
    amrex::ignore_unused(tfactors);

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 = 0.0;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_O14_to_N14_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // O14 --> N14

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // wc12w
    ln_set_rate =  -4.62354;
    amrex::ignore_unused(tfactors);

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 = 0.0;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_O15_to_N15_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // O15 --> N15

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // wc12w
    ln_set_rate =  -5.17053;
    amrex::ignore_unused(tfactors);

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 = 0.0;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_C12_to_N13_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // C12 + p --> N13

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_C12_to_O16_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // C12 + He4 --> O16

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_C13_to_N14_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // C13 + p --> N14

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_N13_to_O14_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // N13 + p --> O14

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_N14_to_O15_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // N14 + p --> O15

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_N14_to_F18_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // N14 + He4 --> F18

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_N15_to_O16_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // N15 + p --> O16

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_N15_to_F19_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // N15 + He4 --> F19

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_O14_to_Ne18_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // O14 + He4 --> Ne18

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_O15_to_Ne19_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // O15 + He4 --> Ne19

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_O16_to_F17_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // O16 + p --> F17

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_O16_to_Ne20_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // O16 + He4 --> Ne20

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_O17_to_F18_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // O17 + p --> F18

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_O18_to_F19_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // O18 + p --> F19

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_F17_to_Ne18_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // F17 + p --> Ne18

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_F18_to_Ne19_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // F18 + p --> Ne19

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_F19_to_Ne20_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // F19 + p --> Ne20

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Ne18_to_Mg22_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ne18 + He4 --> Mg22

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Ne20_to_Mg24_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ne20 + He4 --> Mg24

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_C12_C12_to_He4_Ne20_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // C12 + C12 --> He4 + Ne20

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_N13_to_p_O16_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // N13 + He4 --> p + O16

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_N15_to_He4_C12_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // N15 + p --> He4 + C12

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_O14_to_p_F17_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // O14 + He4 --> p + F17

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_C12_O16_to_He4_Mg24_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // O16 + C12 --> He4 + Mg24

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_O17_to_He4_N14_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // O17 + p --> He4 + N14

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_O18_to_He4_N15_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // O18 + p --> He4 + N15

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_F18_to_He4_O15_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // F18 + p --> He4 + O15

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_F19_to_He4_O16_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // F19 + p --> He4 + O16

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Ne20_to_He4_F17_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ne20 + p --> He4 + F17

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_He4_He4_to_C12_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // 3 He4 --> C12

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_N13_to_p_C12_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // N13 --> p + C12

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.00028189664206047983;
    amrex::Real Q_kBT = -1.9434910639993015 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_C12_to_N13_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_C12_to_N13_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real C12_pf, dC12_pf_dT;
    // setting C12 partition function to 1.0 by default, independent of T
    C12_pf = 1.0_rt;
    dC12_pf_dT = 0.0_rt;

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real N13_pf, dN13_pf_dT;
    // setting N13 partition function to 1.0 by default, independent of T
    N13_pf = 1.0_rt;
    dN13_pf_dT = 0.0_rt;

    amrex::Real z_r = p_pf * C12_pf;
    amrex::Real z_p = N13_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = C12_pf * dp_pf_dT + p_pf * dC12_pf_dT;
        amrex::Real dz_p_dT = dN13_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_N14_to_p_C13_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // N14 --> p + C13

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0003796782582881305;
    amrex::Real Q_kBT = -7.550563563999526 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_C13_to_N14_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_C13_to_N14_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real C13_pf, dC13_pf_dT;
    // setting C13 partition function to 1.0 by default, independent of T
    C13_pf = 1.0_rt;
    dC13_pf_dT = 0.0_rt;

    amrex::Real N14_pf, dN14_pf_dT;
    // setting N14 partition function to 1.0 by default, independent of T
    N14_pf = 1.0_rt;
    dN14_pf_dT = 0.0_rt;

    amrex::Real z_r = p_pf * C13_pf;
    amrex::Real z_p = N14_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = C13_pf * dp_pf_dT + p_pf * dC13_pf_dT;
        amrex::Real dz_p_dT = dN14_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_O14_to_p_N13_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // O14 --> p + N13

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0011384340603262408;
    amrex::Real Q_kBT = -4.626670064000389 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_N13_to_O14_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_N13_to_O14_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real O14_pf, dO14_pf_dT;
    // interpolating O14 partition function
    get_partition_function_cached(O14, tfactors, pf_cache, O14_pf, dO14_pf_dT);

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real N13_pf, dN13_pf_dT;
    // setting N13 partition function to 1.0 by default, independent of T
    N13_pf = 1.0_rt;
    dN13_pf_dT = 0.0_rt;

    amrex::Real z_r = p_pf * N13_pf;
    amrex::Real z_p = O14_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = N13_pf * dp_pf_dT + p_pf * dN13_pf_dT;
        amrex::Real dz_p_dT = dO14_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_O15_to_p_N14_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // O15 --> p + N14

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0008608057335374369;
    amrex::Real Q_kBT = -7.2967878640004145 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_N14_to_O15_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_N14_to_O15_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real O15_pf, dO15_pf_dT;
    // interpolating O15 partition function
    get_partition_function_cached(O15, tfactors, pf_cache, O15_pf, dO15_pf_dT);

    amrex::Real N14_pf, dN14_pf_dT;
    // setting N14 partition function to 1.0 by default, independent of T
    N14_pf = 1.0_rt;
    dN14_pf_dT = 0.0_rt;

    amrex::Real z_r = p_pf * N14_pf;
    amrex::Real z_p = O15_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = N14_pf * dp_pf_dT + p_pf * dN14_pf_dT;
        amrex::Real dz_p_dT = dO15_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_O16_to_p_N15_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // O16 --> p + N15

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0011563309866673933;
    amrex::Real Q_kBT = -12.127411263999507 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_N15_to_O16_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_N15_to_O16_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real N15_pf, dN15_pf_dT;
    // setting N15 partition function to 1.0 by default, independent of T
    N15_pf = 1.0_rt;
    dN15_pf_dT = 0.0_rt;

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real O16_pf, dO16_pf_dT;
    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors, pf_cache, O16_pf, dO16_pf_dT);

    amrex::Real z_r = p_pf * N15_pf;
    amrex::Real z_p = O16_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = N15_pf * dp_pf_dT + p_pf * dN15_pf_dT;
        amrex::Real dz_p_dT = dO16_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_O16_to_He4_C12_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // O16 --> He4 + C12

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0016254843301234132;
    amrex::Real Q_kBT = -7.161917969999712 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_C12_to_O16_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_C12_to_O16_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real C12_pf, dC12_pf_dT;
    // setting C12 partition function to 1.0 by default, independent of T
    C12_pf = 1.0_rt;
    dC12_pf_dT = 0.0_rt;

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real O16_pf, dO16_pf_dT;
    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors, pf_cache, O16_pf, dO16_pf_dT);

    amrex::Real z_r = He4_pf * C12_pf;
    amrex::Real z_p = O16_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = C12_pf * dHe4_pf_dT + He4_pf * dC12_pf_dT;
        amrex::Real dz_p_dT = dO16_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_F17_to_p_O16_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // F17 --> p + O16

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 9.674279321716613e-05;
    amrex::Real Q_kBT = -0.6002689639985703 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_O16_to_F17_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_O16_to_F17_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real F17_pf, dF17_pf_dT;
    // interpolating F17 partition function
    get_partition_function_cached(F17, tfactors, pf_cache, F17_pf, dF17_pf_dT);

    amrex::Real O16_pf, dO16_pf_dT;
    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors, pf_cache, O16_pf, dO16_pf_dT);

    amrex::Real z_r = p_pf * O16_pf;
    amrex::Real z_p = F17_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = O16_pf * dp_pf_dT + p_pf * dO16_pf_dT;
        amrex::Real dz_p_dT = dF17_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_F18_to_p_O17_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // F18 --> p + O17

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0011679516397592826;
    amrex::Real Q_kBT = -5.607106863999662 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_O17_to_F18_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_O17_to_F18_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real F18_pf, dF18_pf_dT;
    // interpolating F18 partition function
    get_partition_function_cached(F18, tfactors, pf_cache, F18_pf, dF18_pf_dT);

    amrex::Real O17_pf, dO17_pf_dT;
    // interpolating O17 partition function
    get_partition_function_cached(O17, tfactors, pf_cache, O17_pf, dO17_pf_dT);

    amrex::Real z_r = p_pf * O17_pf;
    amrex::Real z_p = F18_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = O17_pf * dp_pf_dT + p_pf * dO17_pf_dT;
        amrex::Real dz_p_dT = dF18_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_F18_to_He4_N14_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // F18 --> He4 + N14

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0017159747049146441;
    amrex::Real Q_kBT = -4.415232669999798 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_N14_to_F18_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_N14_to_F18_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real F18_pf, dF18_pf_dT;
    // interpolating F18 partition function
    get_partition_function_cached(F18, tfactors, pf_cache, F18_pf, dF18_pf_dT);

    amrex::Real N14_pf, dN14_pf_dT;
    // setting N14 partition function to 1.0 by default, independent of T
    N14_pf = 1.0_rt;
    dN14_pf_dT = 0.0_rt;

    amrex::Real z_r = He4_pf * N14_pf;
    amrex::Real z_p = F18_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = N14_pf * dHe4_pf_dT + He4_pf * dN14_pf_dT;
        amrex::Real dz_p_dT = dF18_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_F19_to_p_O18_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // F19 --> p + O18

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.00029344799859633736;
    amrex::Real Q_kBT = -7.993599864000986 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_O18_to_F19_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_O18_to_F19_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real O18_pf, dO18_pf_dT;
    // interpolating O18 partition function
    get_partition_function_cached(O18, tfactors, pf_cache, O18_pf, dO18_pf_dT);

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real F19_pf, dF19_pf_dT;
    // interpolating F19 partition function
    get_partition_function_cached(F19, tfactors, pf_cache, F19_pf, dF19_pf_dT);

    amrex::Real z_r = p_pf * O18_pf;
    amrex::Real z_p = F19_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = O18_pf * dp_pf_dT + p_pf * dO18_pf_dT;
        amrex::Real dz_p_dT = dF19_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_F19_to_He4_N15_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // F19 --> He4 + N15

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.001754491925665686;
    amrex::Real Q_kBT = -4.0137990700009265 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_N15_to_F19_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_N15_to_F19_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real N15_pf, dN15_pf_dT;
    // setting N15 partition function to 1.0 by default, independent of T
    N15_pf = 1.0_rt;
    dN15_pf_dT = 0.0_rt;

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real F19_pf, dF19_pf_dT;
    // interpolating F19 partition function
    get_partition_function_cached(F19, tfactors, pf_cache, F19_pf, dF19_pf_dT);

    amrex::Real z_r = He4_pf * N15_pf;
    amrex::Real z_p = F19_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = N15_pf * dHe4_pf_dT + He4_pf * dN15_pf_dT;
        amrex::Real dz_p_dT = dF19_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ne18_to_p_F17_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ne18 --> p + F17

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0035030604144385565;
    amrex::Real Q_kBT = -3.9230710640041138 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_F17_to_Ne18_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_F17_to_Ne18_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real F17_pf, dF17_pf_dT;
    // interpolating F17 partition function
    get_partition_function_cached(F17, tfactors, pf_cache, F17_pf, dF17_pf_dT);

    amrex::Real Ne18_pf, dNe18_pf_dT;
    // interpolating Ne18 partition function
    get_partition_function_cached(Ne18, tfactors, pf_cache, Ne18_pf, dNe18_pf_dT);

    amrex::Real z_r = p_pf * F17_pf;
    amrex::Real z_p = Ne18_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = F17_pf * dp_pf_dT + p_pf * dF17_pf_dT;
        amrex::Real dz_p_dT = dNe18_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ne18_to_He4_O14_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ne18 --> He4 + O14

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0017165294291050784;
    amrex::Real Q_kBT = -5.115096870001259 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_O14_to_Ne18_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_O14_to_Ne18_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real O14_pf, dO14_pf_dT;
    // interpolating O14 partition function
    get_partition_function_cached(O14, tfactors, pf_cache, O14_pf, dO14_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Ne18_pf, dNe18_pf_dT;
    // interpolating Ne18 partition function
    get_partition_function_cached(Ne18, tfactors, pf_cache, Ne18_pf, dNe18_pf_dT);

    amrex::Real z_r = He4_pf * O14_pf;
    amrex::Real z_p = Ne18_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = O14_pf * dHe4_pf_dT + He4_pf * dO14_pf_dT;
        amrex::Real dz_p_dT = dNe18_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ne19_to_p_F18_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ne19 --> p + F18

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0008801585338347593;
    amrex::Real Q_kBT = -6.410021064003507 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_F18_to_Ne19_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_F18_to_Ne19_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real F18_pf, dF18_pf_dT;
    // interpolating F18 partition function
    get_partition_function_cached(F18, tfactors, pf_cache, F18_pf, dF18_pf_dT);

    amrex::Real Ne19_pf, dNe19_pf_dT;
    // interpolating Ne19 partition function
    get_partition_function_cached(Ne19, tfactors, pf_cache, Ne19_pf, dNe19_pf_dT);

    amrex::Real z_r = p_pf * F18_pf;
    amrex::Real z_p = Ne19_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = F18_pf * dp_pf_dT + p_pf * dF18_pf_dT;
        amrex::Real dz_p_dT = dNe19_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ne19_to_He4_O15_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ne19 --> He4 + O15

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0017545535787367311;
    amrex::Real Q_kBT = -3.5284658700012366 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_O15_to_Ne19_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_O15_to_Ne19_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real O15_pf, dO15_pf_dT;
    // interpolating O15 partition function
    get_partition_function_cached(O15, tfactors, pf_cache, O15_pf, dO15_pf_dT);

    amrex::Real Ne19_pf, dNe19_pf_dT;
    // interpolating Ne19 partition function
    get_partition_function_cached(Ne19, tfactors, pf_cache, Ne19_pf, dNe19_pf_dT);

    amrex::Real z_r = He4_pf * O15_pf;
    amrex::Real z_p = Ne19_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = O15_pf * dHe4_pf_dT + He4_pf * dO15_pf_dT;
        amrex::Real dz_p_dT = dNe19_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ne20_to_p_F19_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ne20 --> p + F19

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0011794427262325786;
    amrex::Real Q_kBT = -12.84345816399974 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_F19_to_Ne20_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_F19_to_Ne20_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real F19_pf, dF19_pf_dT;
    // interpolating F19 partition function
    get_partition_function_cached(F19, tfactors, pf_cache, F19_pf, dF19_pf_dT);

    amrex::Real Ne20_pf, dNe20_pf_dT;
    // interpolating Ne20 partition function
    get_partition_function_cached(Ne20, tfactors, pf_cache, Ne20_pf, dNe20_pf_dT);

    amrex::Real z_r = p_pf * F19_pf;
    amrex::Real z_p = Ne20_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = F19_pf * dp_pf_dT + p_pf * dF19_pf_dT;
        amrex::Real dz_p_dT = dNe20_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ne20_to_He4_O16_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ne20 --> He4 + O16

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0017895591866166972;
    amrex::Real Q_kBT = -4.729845970000538 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_O16_to_Ne20_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_O16_to_Ne20_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real O16_pf, dO16_pf_dT;
    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors, pf_cache, O16_pf, dO16_pf_dT);

    amrex::Real Ne20_pf, dNe20_pf_dT;
    // interpolating Ne20 partition function
    get_partition_function_cached(Ne20, tfactors, pf_cache, Ne20_pf, dNe20_pf_dT);

    amrex::Real z_r = He4_pf * O16_pf;
    amrex::Real z_p = Ne20_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = O16_pf * dHe4_pf_dT + He4_pf * dO16_pf_dT;
        amrex::Real dz_p_dT = dNe20_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Mg22_to_He4_Ne18_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mg22 --> He4 + Ne18

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.001852191832471393;
    amrex::Real Q_kBT = -8.142505869997956 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_Ne18_to_Mg22_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_Ne18_to_Mg22_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real Mg22_pf, dMg22_pf_dT;
    // interpolating Mg22 partition function
    get_partition_function_cached(Mg22, tfactors, pf_cache, Mg22_pf, dMg22_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Ne18_pf, dNe18_pf_dT;
    // interpolating Ne18 partition function
    get_partition_function_cached(Ne18, tfactors, pf_cache, Ne18_pf, dNe18_pf_dT);

    amrex::Real z_r = He4_pf * Ne18_pf;
    amrex::Real z_p = Mg22_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Ne18_pf * dHe4_pf_dT + He4_pf * dNe18_pf_dT;
        amrex::Real dz_p_dT = dMg22_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Mg24_to_He4_Ne20_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mg24 --> He4 + Ne20

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0019034449262058997;
    amrex::Real Q_kBT = -9.31656167000121 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_Ne20_to_Mg24_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_Ne20_to_Mg24_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Mg24_pf, dMg24_pf_dT;
    // interpolating Mg24 partition function
    get_partition_function_cached(Mg24, tfactors, pf_cache, Mg24_pf, dMg24_pf_dT);

    amrex::Real Ne20_pf, dNe20_pf_dT;
    // interpolating Ne20 partition function
    get_partition_function_cached(Ne20, tfactors, pf_cache, Ne20_pf, dNe20_pf_dT);

    amrex::Real z_r = He4_pf * Ne20_pf;
    amrex::Real z_p = Mg24_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Ne20_pf * dHe4_pf_dT + He4_pf * dNe20_pf_dT;
        amrex::Real dz_p_dT = dMg24_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_C12_to_He4_He4_He4_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // C12 --> 3 He4

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 2.0089297549344094e-07;
    amrex::Real Q_kBT = -7.274747609999003 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<6>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_He4_He4_to_C12_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (3.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_He4_He4_to_C12_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real C12_pf, dC12_pf_dT;
    // setting C12 partition function to 1.0 by default, independent of T
    C12_pf = 1.0_rt;
    dC12_pf_dT = 0.0_rt;

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real z_r = He4_pf * He4_pf * He4_pf;
    amrex::Real z_p = C12_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = dHe4_pf_dT + dHe4_pf_dT + dHe4_pf_dT;
        amrex::Real dz_p_dT = dC12_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_C12_to_p_N15_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // C12 + He4 --> p + N15

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.7113762742822628;
    amrex::Real Q_kBT = -4.965493293998968 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_N15_to_He4_C12_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_N15_to_He4_C12_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real N15_pf, dN15_pf_dT;
    // setting N15 partition function to 1.0 by default, independent of T
    N15_pf = 1.0_rt;
    dN15_pf_dT = 0.0_rt;

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real C12_pf, dC12_pf_dT;
    // setting C12 partition function to 1.0 by default, independent of T
    C12_pf = 1.0_rt;
    dC12_pf_dT = 0.0_rt;

    amrex::Real z_r = p_pf * N15_pf;
    amrex::Real z_p = He4_pf * C12_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = N15_pf * dp_pf_dT + p_pf * dN15_pf_dT;
        amrex::Real dz_p_dT = C12_pf * dHe4_pf_dT + He4_pf * dC12_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_N14_to_p_O17_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // N14 + He4 --> p + O17

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.6806345317411766;
    amrex::Real Q_kBT = -1.1918741939998632 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_O17_to_He4_N14_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_O17_to_He4_N14_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real O17_pf, dO17_pf_dT;
    // interpolating O17 partition function
    get_partition_function_cached(O17, tfactors, pf_cache, O17_pf, dO17_pf_dT);

    amrex::Real N14_pf, dN14_pf_dT;
    // setting N14 partition function to 1.0 by default, independent of T
    N14_pf = 1.0_rt;
    dN14_pf_dT = 0.0_rt;

    amrex::Real z_r = p_pf * O17_pf;
    amrex::Real z_p = He4_pf * N14_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = O17_pf * dp_pf_dT + p_pf * dO17_pf_dT;
        amrex::Real dz_p_dT = N14_pf * dHe4_pf_dT + He4_pf * dN14_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_N15_to_p_O18_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // N15 + He4 --> p + O18

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.16725525737885502;
    amrex::Real Q_kBT = -3.979800794000887 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_O18_to_He4_N15_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_O18_to_He4_N15_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real O18_pf, dO18_pf_dT;
    // interpolating O18 partition function
    get_partition_function_cached(O18, tfactors, pf_cache, O18_pf, dO18_pf_dT);

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real N15_pf, dN15_pf_dT;
    // setting N15 partition function to 1.0 by default, independent of T
    N15_pf = 1.0_rt;
    dN15_pf_dT = 0.0_rt;

    amrex::Real z_r = p_pf * O18_pf;
    amrex::Real z_p = He4_pf * N15_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = O18_pf * dp_pf_dT + p_pf * dO18_pf_dT;
        amrex::Real dz_p_dT = N15_pf * dHe4_pf_dT + He4_pf * dN15_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_O15_to_p_F18_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // O15 + He4 --> p + F18

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.5016424374275699;
    amrex::Real Q_kBT = -2.881555194002271 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_F18_to_He4_O15_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_F18_to_He4_O15_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real O15_pf, dO15_pf_dT;
    // interpolating O15 partition function
    get_partition_function_cached(O15, tfactors, pf_cache, O15_pf, dO15_pf_dT);

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real F18_pf, dF18_pf_dT;
    // interpolating F18 partition function
    get_partition_function_cached(F18, tfactors, pf_cache, F18_pf, dF18_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real z_r = p_pf * F18_pf;
    amrex::Real z_p = He4_pf * O15_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = F18_pf * dp_pf_dT + p_pf * dF18_pf_dT;
        amrex::Real dz_p_dT = O15_pf * dHe4_pf_dT + He4_pf * dO15_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_O16_to_He4_N13_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // O16 + p --> He4 + N13

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 5.766242259014467;
    amrex::Real Q_kBT = -5.218426906001445 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_N13_to_p_O16_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_N13_to_p_O16_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real O16_pf, dO16_pf_dT;
    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors, pf_cache, O16_pf, dO16_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real N13_pf, dN13_pf_dT;
    // setting N13 partition function to 1.0 by default, independent of T
    N13_pf = 1.0_rt;
    dN13_pf_dT = 0.0_rt;

    amrex::Real z_r = He4_pf * N13_pf;
    amrex::Real z_p = p_pf * O16_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = N13_pf * dHe4_pf_dT + He4_pf * dN13_pf_dT;
        amrex::Real dz_p_dT = O16_pf * dp_pf_dT + p_pf * dO16_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_O16_to_p_F19_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // O16 + He4 --> p + F19

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.659068856203861;
    amrex::Real Q_kBT = -8.113612194000856 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_F19_to_He4_O16_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_F19_to_He4_O16_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real O16_pf, dO16_pf_dT;
    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors, pf_cache, O16_pf, dO16_pf_dT);

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real F19_pf, dF19_pf_dT;
    // interpolating F19 partition function
    get_partition_function_cached(F19, tfactors, pf_cache, F19_pf, dF19_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real z_r = p_pf * F19_pf;
    amrex::Real z_p = He4_pf * O16_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = F19_pf * dp_pf_dT + p_pf * dF19_pf_dT;
        amrex::Real dz_p_dT = O16_pf * dHe4_pf_dT + He4_pf * dO16_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_F17_to_He4_O14_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // F17 + p --> He4 + O14

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.49000851427799036;
    amrex::Real Q_kBT = -1.192025805997145 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_O14_to_p_F17_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_O14_to_p_F17_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real O14_pf, dO14_pf_dT;
    // interpolating O14 partition function
    get_partition_function_cached(O14, tfactors, pf_cache, O14_pf, dO14_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real F17_pf, dF17_pf_dT;
    // interpolating F17 partition function
    get_partition_function_cached(F17, tfactors, pf_cache, F17_pf, dF17_pf_dT);

    amrex::Real z_r = He4_pf * O14_pf;
    amrex::Real z_p = p_pf * F17_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = O14_pf * dHe4_pf_dT + He4_pf * dO14_pf_dT;
        amrex::Real dz_p_dT = F17_pf * dp_pf_dT + p_pf * dF17_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_F17_to_p_Ne20_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // F17 + He4 --> p + Ne20

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.05405956614380885;
    amrex::Real Q_kBT = 4.1295770060003125 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_Ne20_to_He4_F17_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_Ne20_to_He4_F17_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real F17_pf, dF17_pf_dT;
    // interpolating F17 partition function
    get_partition_function_cached(F17, tfactors, pf_cache, F17_pf, dF17_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Ne20_pf, dNe20_pf_dT;
    // interpolating Ne20 partition function
    get_partition_function_cached(Ne20, tfactors, pf_cache, Ne20_pf, dNe20_pf_dT);

    amrex::Real z_r = p_pf * Ne20_pf;
    amrex::Real z_p = He4_pf * F17_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Ne20_pf * dp_pf_dT + p_pf * dNe20_pf_dT;
        amrex::Real dz_p_dT = F17_pf * dHe4_pf_dT + He4_pf * dF17_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Ne20_to_C12_C12_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ne20 + He4 --> C12 + C12

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 1.2066542763748098;
    amrex::Real Q_kBT = -4.617016329998765 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_C12_C12_to_He4_Ne20_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_C12_C12_to_He4_Ne20_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real C12_pf, dC12_pf_dT;
    // setting C12 partition function to 1.0 by default, independent of T
    C12_pf = 1.0_rt;
    dC12_pf_dT = 0.0_rt;

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Ne20_pf, dNe20_pf_dT;
    // interpolating Ne20 partition function
    get_partition_function_cached(Ne20, tfactors, pf_cache, Ne20_pf, dNe20_pf_dT);

    amrex::Real z_r = C12_pf * C12_pf;
    amrex::Real z_p = He4_pf * Ne20_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = dC12_pf_dT + dC12_pf_dT;
        amrex::Real dz_p_dT = Ne20_pf * dHe4_pf_dT + He4_pf * dNe20_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Mg24_to_C12_O16_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mg24 + He4 --> C12 + O16

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 2.8259884361677;
    amrex::Real Q_kBT = -6.771660029999435 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_C12_O16_to_He4_Mg24_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_C12_O16_to_He4_Mg24_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real C12_pf, dC12_pf_dT;
    // setting C12 partition function to 1.0 by default, independent of T
    C12_pf = 1.0_rt;
    dC12_pf_dT = 0.0_rt;

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Mg24_pf, dMg24_pf_dT;
    // interpolating Mg24 partition function
    get_partition_function_cached(Mg24, tfactors, pf_cache, Mg24_pf, dMg24_pf_dT);

    amrex::Real O16_pf, dO16_pf_dT;
    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors, pf_cache, O16_pf, dO16_pf_dT);

    amrex::Real z_r = C12_pf * O16_pf;
    amrex::Real z_p = He4_pf * Mg24_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = O16_pf * dC12_pf_dT + C12_pf * dO16_pf_dT;
        amrex::Real dz_p_dT = Mg24_pf * dHe4_pf_dT + He4_pf * dMg24_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}



template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
fill_reaclib_rates(const tf_t& tfactors, T& rate_eval)
{

    amrex::Real rate;
    amrex::Real drate_dT;

    rate_N13_to_C13_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_N13_to_C13_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_N13_to_C13_reaclib) = drate_dT;

    }
    rate_O14_to_N14_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_O14_to_N14_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_O14_to_N14_reaclib) = drate_dT;

    }
    rate_O15_to_N15_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_O15_to_N15_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_O15_to_N15_reaclib) = drate_dT;

    }
    rate_p_C12_to_N13_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_C12_to_N13_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_C12_to_N13_reaclib) = drate_dT;

    }
    rate_He4_C12_to_O16_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_C12_to_O16_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_C12_to_O16_reaclib) = drate_dT;

    }
    rate_p_C13_to_N14_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_C13_to_N14_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_C13_to_N14_reaclib) = drate_dT;

    }
    rate_p_N13_to_O14_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_N13_to_O14_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_N13_to_O14_reaclib) = drate_dT;

    }
    rate_p_N14_to_O15_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_N14_to_O15_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_N14_to_O15_reaclib) = drate_dT;

    }
    rate_He4_N14_to_F18_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_N14_to_F18_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_N14_to_F18_reaclib) = drate_dT;

    }
    rate_p_N15_to_O16_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_N15_to_O16_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_N15_to_O16_reaclib) = drate_dT;

    }
    rate_He4_N15_to_F19_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_N15_to_F19_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_N15_to_F19_reaclib) = drate_dT;

    }
    rate_He4_O14_to_Ne18_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_O14_to_Ne18_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_O14_to_Ne18_reaclib) = drate_dT;

    }
    rate_He4_O15_to_Ne19_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_O15_to_Ne19_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_O15_to_Ne19_reaclib) = drate_dT;

    }
    rate_p_O16_to_F17_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_O16_to_F17_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_O16_to_F17_reaclib) = drate_dT;

    }
    rate_He4_O16_to_Ne20_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_O16_to_Ne20_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_O16_to_Ne20_reaclib) = drate_dT;

    }
    rate_p_O17_to_F18_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_O17_to_F18_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_O17_to_F18_reaclib) = drate_dT;

    }
    rate_p_O18_to_F19_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_O18_to_F19_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_O18_to_F19_reaclib) = drate_dT;

    }
    rate_p_F17_to_Ne18_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_F17_to_Ne18_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_F17_to_Ne18_reaclib) = drate_dT;

    }
    rate_p_F18_to_Ne19_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_F18_to_Ne19_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_F18_to_Ne19_reaclib) = drate_dT;

    }
    rate_p_F19_to_Ne20_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_F19_to_Ne20_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_F19_to_Ne20_reaclib) = drate_dT;

    }
    rate_He4_Ne18_to_Mg22_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Ne18_to_Mg22_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Ne18_to_Mg22_reaclib) = drate_dT;

    }
    rate_He4_Ne20_to_Mg24_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Ne20_to_Mg24_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Ne20_to_Mg24_reaclib) = drate_dT;

    }
    rate_C12_C12_to_He4_Ne20_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_C12_C12_to_He4_Ne20_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_C12_C12_to_He4_Ne20_reaclib) = drate_dT;

    }
    rate_He4_N13_to_p_O16_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_N13_to_p_O16_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_N13_to_p_O16_reaclib) = drate_dT;

    }
    rate_p_N15_to_He4_C12_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_N15_to_He4_C12_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_N15_to_He4_C12_reaclib) = drate_dT;

    }
    rate_He4_O14_to_p_F17_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_O14_to_p_F17_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_O14_to_p_F17_reaclib) = drate_dT;

    }
    rate_C12_O16_to_He4_Mg24_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_C12_O16_to_He4_Mg24_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_C12_O16_to_He4_Mg24_reaclib) = drate_dT;

    }
    rate_p_O17_to_He4_N14_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_O17_to_He4_N14_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_O17_to_He4_N14_reaclib) = drate_dT;

    }
    rate_p_O18_to_He4_N15_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_O18_to_He4_N15_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_O18_to_He4_N15_reaclib) = drate_dT;

    }
    rate_p_F18_to_He4_O15_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_F18_to_He4_O15_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_F18_to_He4_O15_reaclib) = drate_dT;

    }
    rate_p_F19_to_He4_O16_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_F19_to_He4_O16_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_F19_to_He4_O16_reaclib) = drate_dT;

    }
    rate_p_Ne20_to_He4_F17_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Ne20_to_He4_F17_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Ne20_to_He4_F17_reaclib) = drate_dT;

    }
    rate_He4_He4_He4_to_C12_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_He4_He4_to_C12_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_He4_He4_to_C12_reaclib) = drate_dT;

    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
fill_derived_rates(const tf_t& tfactors, T& rate_eval)
{

    amrex::Real rate;
    amrex::Real drate_dT;

    part_fun::pf_cache_t pf_cache{};

    pf_cache.index_temp_array_1 = interp_net::find_index(tfactors.T9, part_fun::temp_array_1);

    rate_N13_to_p_C12_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_N13_to_p_C12_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_N13_to_p_C12_derived) = drate_dT;

    }
    rate_N14_to_p_C13_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_N14_to_p_C13_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_N14_to_p_C13_derived) = drate_dT;

    }
    rate_O14_to_p_N13_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_O14_to_p_N13_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_O14_to_p_N13_derived) = drate_dT;

    }
    rate_O15_to_p_N14_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_O15_to_p_N14_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_O15_to_p_N14_derived) = drate_dT;

    }
    rate_O16_to_p_N15_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_O16_to_p_N15_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_O16_to_p_N15_derived) = drate_dT;

    }
    rate_O16_to_He4_C12_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_O16_to_He4_C12_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_O16_to_He4_C12_derived) = drate_dT;

    }
    rate_F17_to_p_O16_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_F17_to_p_O16_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_F17_to_p_O16_derived) = drate_dT;

    }
    rate_F18_to_p_O17_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_F18_to_p_O17_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_F18_to_p_O17_derived) = drate_dT;

    }
    rate_F18_to_He4_N14_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_F18_to_He4_N14_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_F18_to_He4_N14_derived) = drate_dT;

    }
    rate_F19_to_p_O18_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_F19_to_p_O18_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_F19_to_p_O18_derived) = drate_dT;

    }
    rate_F19_to_He4_N15_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_F19_to_He4_N15_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_F19_to_He4_N15_derived) = drate_dT;

    }
    rate_Ne18_to_p_F17_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Ne18_to_p_F17_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ne18_to_p_F17_derived) = drate_dT;

    }
    rate_Ne18_to_He4_O14_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Ne18_to_He4_O14_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ne18_to_He4_O14_derived) = drate_dT;

    }
    rate_Ne19_to_p_F18_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Ne19_to_p_F18_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ne19_to_p_F18_derived) = drate_dT;

    }
    rate_Ne19_to_He4_O15_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Ne19_to_He4_O15_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ne19_to_He4_O15_derived) = drate_dT;

    }
    rate_Ne20_to_p_F19_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Ne20_to_p_F19_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ne20_to_p_F19_derived) = drate_dT;

    }
    rate_Ne20_to_He4_O16_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Ne20_to_He4_O16_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ne20_to_He4_O16_derived) = drate_dT;

    }
    rate_Mg22_to_He4_Ne18_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Mg22_to_He4_Ne18_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Mg22_to_He4_Ne18_derived) = drate_dT;

    }
    rate_Mg24_to_He4_Ne20_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Mg24_to_He4_Ne20_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Mg24_to_He4_Ne20_derived) = drate_dT;

    }
    rate_C12_to_He4_He4_He4_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_C12_to_He4_He4_He4_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_C12_to_He4_He4_He4_derived) = drate_dT;

    }
    rate_He4_C12_to_p_N15_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_He4_C12_to_p_N15_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_C12_to_p_N15_derived) = drate_dT;

    }
    rate_He4_N14_to_p_O17_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_He4_N14_to_p_O17_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_N14_to_p_O17_derived) = drate_dT;

    }
    rate_He4_N15_to_p_O18_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_He4_N15_to_p_O18_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_N15_to_p_O18_derived) = drate_dT;

    }
    rate_He4_O15_to_p_F18_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_He4_O15_to_p_F18_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_O15_to_p_F18_derived) = drate_dT;

    }
    rate_p_O16_to_He4_N13_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_p_O16_to_He4_N13_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_O16_to_He4_N13_derived) = drate_dT;

    }
    rate_He4_O16_to_p_F19_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_He4_O16_to_p_F19_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_O16_to_p_F19_derived) = drate_dT;

    }
    rate_p_F17_to_He4_O14_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_p_F17_to_He4_O14_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_F17_to_He4_O14_derived) = drate_dT;

    }
    rate_He4_F17_to_p_Ne20_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_He4_F17_to_p_Ne20_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_F17_to_p_Ne20_derived) = drate_dT;

    }
    rate_He4_Ne20_to_C12_C12_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_He4_Ne20_to_C12_C12_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Ne20_to_C12_C12_derived) = drate_dT;

    }
    rate_He4_Mg24_to_C12_O16_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_He4_Mg24_to_C12_O16_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Mg24_to_C12_O16_derived) = drate_dT;

    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
fill_approx_rates([[maybe_unused]] const tf_t& tfactors,
                  [[maybe_unused]] const amrex::Real rho,
                  [[maybe_unused]] const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
                  [[maybe_unused]] T& rate_eval)
{

    [[maybe_unused]] amrex::Real rate{};
    [[maybe_unused]] amrex::Real drate_dT{};


}

#endif
