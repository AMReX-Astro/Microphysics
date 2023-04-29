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
void rate_n__p__weak__wc12(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // n --> p

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // wc12w
    ln_set_rate =  -6.78161;

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
void rate_n_n13__n14(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // n13 + n --> n14

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // wiesr
    ln_set_rate =  -3.63074 + -2.99547 * tfactors.T9i + 9.44873e-10 * tfactors.T913i + -2.33713e-09 * tfactors.T913
                         + 1.97507e-10 * tfactors.T9 + -1.49747e-11 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  2.99547 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 9.44873e-10 * tfactors.T943i + (1.0/3.0) * -2.33713e-09 * tfactors.T923i
                                  + 1.97507e-10 + (5.0/3.0) * -1.49747e-11 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // wiesn
    ln_set_rate =  13.9377 + -0.0054652 * tfactors.T9i + 1.72241e-10 * tfactors.T913i + -5.62522e-10 * tfactors.T913
                         + 5.59212e-11 * tfactors.T9 + -4.6549e-12 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.0054652 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 1.72241e-10 * tfactors.T943i + (1.0/3.0) * -5.62522e-10 * tfactors.T923i
                                  + 5.59212e-11 + (5.0/3.0) * -4.6549e-12 * tfactors.T923;
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
void rate_he4_f18__na22(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // f18 + he4 --> na22

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // rpsmr
    ln_set_rate =  35.3786 + -1.82957 * tfactors.T9i + 18.8956 * tfactors.T913i + -65.6134 * tfactors.T913
                         + 1.71114 * tfactors.T9 + -0.0260999 * tfactors.T953 + 37.8396 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.82957 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 18.8956 * tfactors.T943i + (1.0/3.0) * -65.6134 * tfactors.T923i
                                  + 1.71114 + (5.0/3.0) * -0.0260999 * tfactors.T923 + 37.8396 * tfactors.T9i;
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
void rate_n_ne20__ne21(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ne20 + n --> ne21

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ka02r
    ln_set_rate =  12.7344 + -1.70393 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.70393 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // ka02n
    ln_set_rate =  8.57651;

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
void rate_p_ne21__na22(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ne21 + p --> na22

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10r
    ln_set_rate =  -47.6554 + -0.19618 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.19618 * tfactors.T9i * tfactors.T9i
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
    ln_set_rate =  19.0696 + -19.2096 * tfactors.T913i
                         + -1.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -19.2096 * tfactors.T943i
                                  + (5.0/3.0) * -1.0 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  -39.4862 + -4.21385 * tfactors.T9i + 21.1176 * tfactors.T913i + 34.0411 * tfactors.T913
                         + -4.45593 * tfactors.T9 + 0.328613 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  4.21385 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 21.1176 * tfactors.T943i + (1.0/3.0) * 34.0411 * tfactors.T923i
                                  + -4.45593 + (5.0/3.0) * 0.328613 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  1.75704 + -1.39957 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.39957 * tfactors.T9i * tfactors.T9i
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
void rate_n_na22__na23(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // na22 + n --> na23

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  12.8678 + 1.02148 * tfactors.T913
                         + -0.334638 * tfactors.T9 + 0.0258708 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * 1.02148 * tfactors.T923i
                                  + -0.334638 + (5.0/3.0) * 0.0258708 * tfactors.T923;
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
void rate_p_na23__mg24(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // na23 + p --> mg24

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10n
    ln_set_rate =  18.9075 + -20.6428 * tfactors.T913i + 1.52954 * tfactors.T913
                         + 2.7487 * tfactors.T9 + -1.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -20.6428 * tfactors.T943i + (1.0/3.0) * 1.52954 * tfactors.T923i
                                  + 2.7487 + (5.0/3.0) * -1.0 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  9.0594 + -3.28029 * tfactors.T9i + -0.360588 * tfactors.T913
                         + 1.4187 * tfactors.T9 + -0.184061 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  3.28029 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.360588 * tfactors.T923i
                                  + 1.4187 + (5.0/3.0) * -0.184061 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  -5.02585 + -1.61219 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.61219 * tfactors.T9i * tfactors.T9i
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
void rate_he4_mg24__si28(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // mg24 + he4 --> si28

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // st08r
    ln_set_rate =  -50.5494 + -12.8332 * tfactors.T9i + 21.3721 * tfactors.T913i + 37.7649 * tfactors.T913
                         + -4.10635 * tfactors.T9 + 0.249618 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  12.8332 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 21.3721 * tfactors.T943i + (1.0/3.0) * 37.7649 * tfactors.T923i
                                  + -4.10635 + (5.0/3.0) * 0.249618 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // st08r
    ln_set_rate =  8.03977 + -15.629 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  15.629 * tfactors.T9i * tfactors.T9i
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
void rate_p_al27__si28(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // al27 + p --> si28

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10r
    ln_set_rate =  -13.6664 + -1.90396 * tfactors.T9i + 23.8634 * tfactors.T913
                         + -3.70135 * tfactors.T9 + 0.28964 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.90396 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 23.8634 * tfactors.T923i
                                  + -3.70135 + (5.0/3.0) * 0.28964 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  86.0234 + -0.387313 * tfactors.T9i + -26.8327 * tfactors.T913i + -116.137 * tfactors.T913
                         + 0.00950567 * tfactors.T9 + 0.00999755 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.387313 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -26.8327 * tfactors.T943i + (1.0/3.0) * -116.137 * tfactors.T923i
                                  + 0.00950567 + (5.0/3.0) * 0.00999755 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10n
    ln_set_rate =  21.1065 + -23.2205 * tfactors.T913i
                         + -2.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -23.2205 * tfactors.T943i
                                  + (5.0/3.0) * -2.0 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_he4_si28__s32(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // si28 + he4 --> s32

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  47.9212 + -59.4896 * tfactors.T913i + 4.47205 * tfactors.T913
                         + -4.78989 * tfactors.T9 + 0.557201 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -59.4896 * tfactors.T943i + (1.0/3.0) * 4.47205 * tfactors.T923i
                                  + -4.78989 + (5.0/3.0) * 0.557201 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_p_p31__s32(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // p31 + p --> s32

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10r
    ln_set_rate =  0.821556 + -3.77704 * tfactors.T9i + 8.09341 * tfactors.T913
                         + -0.615971 * tfactors.T9 + 0.031159 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  3.77704 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 8.09341 * tfactors.T923i
                                  + -0.615971 + (5.0/3.0) * 0.031159 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  -2.66839 + -2.25958 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  2.25958 * tfactors.T9i * tfactors.T9i
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
    ln_set_rate =  19.2596 + -25.3278 * tfactors.T913i + 6.4931 * tfactors.T913
                         + -9.27513 * tfactors.T9 + -0.610439 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -25.3278 * tfactors.T943i + (1.0/3.0) * 6.4931 * tfactors.T923i
                                  + -9.27513 + (5.0/3.0) * -0.610439 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_c12_c12__p_na23(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // c12 + c12 --> p + na23

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // cf88r
    ln_set_rate =  60.9649 + -84.165 * tfactors.T913i + -1.4191 * tfactors.T913
                         + -0.114619 * tfactors.T9 + -0.070307 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -84.165 * tfactors.T943i + (1.0/3.0) * -1.4191 * tfactors.T923i
                                  + -0.114619 + (5.0/3.0) * -0.070307 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_c12_o16__p_al27(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // o16 + c12 --> p + al27

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // cf88r
    ln_set_rate =  68.5253 + 0.205134 * tfactors.T9i + -119.242 * tfactors.T913i + 13.3667 * tfactors.T913
                         + 0.295425 * tfactors.T9 + -0.267288 * tfactors.T953 + -9.91729 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  -0.205134 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -119.242 * tfactors.T943i + (1.0/3.0) * 13.3667 * tfactors.T923i
                                  + 0.295425 + (5.0/3.0) * -0.267288 * tfactors.T923 + -9.91729 * tfactors.T9i;
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
void rate_o16_o16__p_p31(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // o16 + o16 --> p + p31

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // cf88r
    ln_set_rate =  85.2628 + 0.223453 * tfactors.T9i + -145.844 * tfactors.T913i + 8.72612 * tfactors.T913
                         + -0.554035 * tfactors.T9 + -0.137562 * tfactors.T953 + -6.88807 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  -0.223453 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -145.844 * tfactors.T943i + (1.0/3.0) * 8.72612 * tfactors.T923i
                                  + -0.554035 + (5.0/3.0) * -0.137562 * tfactors.T923 + -6.88807 * tfactors.T9i;
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
void rate_o16_o16__he4_si28(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // o16 + o16 --> he4 + si28

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // cf88r
    ln_set_rate =  97.2435 + -0.268514 * tfactors.T9i + -119.324 * tfactors.T913i + -32.2497 * tfactors.T913
                         + 1.46214 * tfactors.T9 + -0.200893 * tfactors.T953 + 13.2148 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.268514 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -119.324 * tfactors.T943i + (1.0/3.0) * -32.2497 * tfactors.T923i
                                  + 1.46214 + (5.0/3.0) * -0.200893 * tfactors.T923 + 13.2148 * tfactors.T9i;
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
void rate_he4_f18__p_ne21(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // f18 + he4 --> p + ne21

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // rpsmr
    ln_set_rate =  49.7863 + -1.84559 * tfactors.T9i + 21.4461 * tfactors.T913i + -73.252 * tfactors.T913
                         + 2.42329 * tfactors.T9 + -0.077278 * tfactors.T953 + 40.7604 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.84559 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 21.4461 * tfactors.T943i + (1.0/3.0) * -73.252 * tfactors.T923i
                                  + 2.42329 + (5.0/3.0) * -0.077278 * tfactors.T923 + 40.7604 * tfactors.T9i;
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
void rate_he4_ne21__n_mg24(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ne21 + he4 --> n + mg24

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // nacrn
    ln_set_rate =  43.9762 + -46.88 * tfactors.T913i + -0.536629 * tfactors.T913
                         + 0.144715 * tfactors.T9 + -0.197624 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -46.88 * tfactors.T943i + (1.0/3.0) * -0.536629 * tfactors.T923i
                                  + 0.144715 + (5.0/3.0) * -0.197624 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // nacrr
    ln_set_rate =  -7.26831 + -13.2638 * tfactors.T9i + 18.0748 * tfactors.T913
                         + -0.981883 * tfactors.T9 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  13.2638 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 18.0748 * tfactors.T923i
                                  + -0.981883 + -1.5 * tfactors.T9i;
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
void rate_p_na23__he4_ne20(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // na23 + p --> he4 + ne20

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10r
    ln_set_rate =  -6.58736 + -2.31577 * tfactors.T9i + 19.7297 * tfactors.T913
                         + -2.20987 * tfactors.T9 + 0.153374 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  2.31577 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 19.7297 * tfactors.T923i
                                  + -2.20987 + (5.0/3.0) * 0.153374 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  0.0178295 + -1.86103 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.86103 * tfactors.T9i * tfactors.T9i
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
    ln_set_rate =  18.9756 + -20.0024 * tfactors.T913i + 11.5988 * tfactors.T913
                         + -1.37398 * tfactors.T9 + -1.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -20.0024 * tfactors.T943i + (1.0/3.0) * 11.5988 * tfactors.T923i
                                  + -1.37398 + (5.0/3.0) * -1.0 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_p_al27__he4_mg24(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // al27 + p --> he4 + mg24

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10r
    ln_set_rate =  -7.02789 + -4.2425 * tfactors.T9i + 18.0416 * tfactors.T913
                         + -1.54137 * tfactors.T9 + 0.0847506 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  4.2425 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 18.0416 * tfactors.T923i
                                  + -1.54137 + (5.0/3.0) * 0.0847506 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  -26.8683 + -0.963012 * tfactors.T9i + 5.18642 * tfactors.T913i + -34.7936 * tfactors.T913
                         + 168.225 * tfactors.T9 + -115.825 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.963012 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 5.18642 * tfactors.T943i + (1.0/3.0) * -34.7936 * tfactors.T923i
                                  + 168.225 + (5.0/3.0) * -115.825 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10n
    ln_set_rate =  29.4576 + -26.4162 * tfactors.T913i
                         + -2.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -26.4162 * tfactors.T943i
                                  + (5.0/3.0) * -2.0 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_p_p31__he4_si28(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // p31 + p --> he4 + si28

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10r
    ln_set_rate =  -10.893 + -3.42575 * tfactors.T9i + 21.521 * tfactors.T913
                         + -1.90355 * tfactors.T9 + 0.092724 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  3.42575 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 21.521 * tfactors.T923i
                                  + -1.90355 + (5.0/3.0) * 0.092724 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  -12.919 + -1.87716 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.87716 * tfactors.T9i * tfactors.T9i
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
    ln_set_rate =  60.8829 + -31.932 * tfactors.T913i + -77.0334 * tfactors.T913
                         + -43.6847 * tfactors.T9 + -4.28955 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -31.932 * tfactors.T943i + (1.0/3.0) * -77.0334 * tfactors.T923i
                                  + -43.6847 + (5.0/3.0) * -4.28955 * tfactors.T923 + -0.666667 * tfactors.T9i;
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

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_c12_c12__mg24__modified(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // c12 + c12 --> mg24

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // cf88r
    ln_set_rate =  -12.8056 + -30.1498 * tfactors.T9i + 11.4826 * tfactors.T913
                         + 1.82849 * tfactors.T9 + -0.34844 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  30.1498 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 11.4826 * tfactors.T923i
                                  + 1.82849 + (5.0/3.0) * -0.34844 * tfactors.T923;
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
void rate_o16_o16__s32__modified(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // o16 + o16 --> s32

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // cf88r
    ln_set_rate =  77.5491 + -0.373641 * tfactors.T9i + -120.83 * tfactors.T913i + -7.72334 * tfactors.T913
                         + -2.27939 * tfactors.T9 + 0.167655 * tfactors.T953 + 7.62001 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.373641 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -120.83 * tfactors.T943i + (1.0/3.0) * -7.72334 * tfactors.T923i
                                  + -2.27939 + (5.0/3.0) * 0.167655 * tfactors.T923 + 7.62001 * tfactors.T9i;
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
void rate_c12_o16__si28__modified(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // o16 + c12 --> si28

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // cf88r
    ln_set_rate =  -132.213 + -1.46479 * tfactors.T9i + -293.089 * tfactors.T913i + 414.404 * tfactors.T913
                         + -28.0562 * tfactors.T9 + 1.61807 * tfactors.T953 + -178.28 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.46479 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -293.089 * tfactors.T943i + (1.0/3.0) * 414.404 * tfactors.T923i
                                  + -28.0562 + (5.0/3.0) * 1.61807 * tfactors.T923 + -178.28 * tfactors.T9i;
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
void rate_he4_s32__ar36__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // s32 + he4 --> ar36

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  48.901 + -65.3709 * tfactors.T913i + 5.68294 * tfactors.T913
                         + -5.00388 * tfactors.T9 + 0.571407 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -65.3709 * tfactors.T943i + (1.0/3.0) * 5.68294 * tfactors.T923i
                                  + -5.00388 + (5.0/3.0) * 0.571407 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_p_cl35__ar36__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // cl35 + p --> ar36

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10r
    ln_set_rate =  -9.03294 + -2.00996 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  2.00996 * tfactors.T9i * tfactors.T9i
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
    ln_set_rate =  -42.5249 + -0.564651 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.564651 * tfactors.T9i * tfactors.T9i
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
    ln_set_rate =  35.6868 + -27.8971 * tfactors.T913i + -16.2304 * tfactors.T913
                         + 35.255 * tfactors.T9 + -25.8411 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -27.8971 * tfactors.T943i + (1.0/3.0) * -16.2304 * tfactors.T923i
                                  + 35.255 + (5.0/3.0) * -25.8411 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  -7.84699 + -3.65092 * tfactors.T9i + 18.0179 * tfactors.T913
                         + -2.86304 * tfactors.T9 + 0.250854 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  3.65092 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 18.0179 * tfactors.T923i
                                  + -2.86304 + (5.0/3.0) * 0.250854 * tfactors.T923 + -1.5 * tfactors.T9i;
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
void rate_p_cl35__he4_s32__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // cl35 + p --> he4 + s32

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10r
    ln_set_rate =  2.29121 + -6.00976 * tfactors.T9i + 5.33756 * tfactors.T913
                         + 1.64418 * tfactors.T9 + -0.246167 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  6.00976 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 5.33756 * tfactors.T923i
                                  + 1.64418 + (5.0/3.0) * -0.246167 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10r
    ln_set_rate =  -1.01202 + -3.93495 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  3.93495 * tfactors.T9i * tfactors.T9i
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
    ln_set_rate =  -57.5294 + -0.532931 * tfactors.T9i + 25.5338 * tfactors.T913
                         + 6.45824 * tfactors.T9 + -0.950294 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.532931 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 25.5338 * tfactors.T923i
                                  + 6.45824 + (5.0/3.0) * -0.950294 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10n
    ln_set_rate =  32.12 + -30.9147 * tfactors.T913i + -1.2345 * tfactors.T913
                         + 22.5118 * tfactors.T9 + -33.0589 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -30.9147 * tfactors.T943i + (1.0/3.0) * -1.2345 * tfactors.T923i
                                  + 22.5118 + (5.0/3.0) * -33.0589 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_he4_ar36__ca40__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ar36 + he4 --> ca40

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  52.3486 + -71.0046 * tfactors.T913i + 4.0656 * tfactors.T913
                         + -5.26509 * tfactors.T9 + 0.683546 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -71.0046 * tfactors.T943i + (1.0/3.0) * 4.0656 * tfactors.T923i
                                  + -5.26509 + (5.0/3.0) * 0.683546 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_p_k39__ca40__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // k39 + p --> ca40

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // lo18r
    ln_set_rate =  2761.38 + -5.22234 * tfactors.T9i + 802.18 * tfactors.T913i + -4010.27 * tfactors.T913
                         + 1136.19 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  5.22234 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 802.18 * tfactors.T943i + (1.0/3.0) * -4010.27 * tfactors.T923i
                                  + 1136.19 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // lo18r
    ln_set_rate =  588.099 + -12.5647 * tfactors.T9i + 641.844 * tfactors.T913i + -1248.49 * tfactors.T913
                         + 564.926 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  12.5647 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 641.844 * tfactors.T943i + (1.0/3.0) * -1248.49 * tfactors.T923i
                                  + 564.926 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // lo18r
    ln_set_rate =  102.252 + -1.66508 * tfactors.T9i + 41.1723 * tfactors.T913i + -149.299 * tfactors.T913
                         + 10.5229 * tfactors.T9 + -0.68208 * tfactors.T953 + 59.2367 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.66508 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 41.1723 * tfactors.T943i + (1.0/3.0) * -149.299 * tfactors.T923i
                                  + 10.5229 + (5.0/3.0) * -0.68208 * tfactors.T923 + 59.2367 * tfactors.T9i;
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
void rate_p_k39__he4_ar36__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // k39 + p --> he4 + ar36

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  20.5166 + -30.0732 * tfactors.T913i + 7.03263 * tfactors.T913
                         + -1.10085 * tfactors.T9 + 0.133768 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -30.0732 * tfactors.T943i + (1.0/3.0) * 7.03263 * tfactors.T923i
                                  + -1.10085 + (5.0/3.0) * 0.133768 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_he4_ca40__ti44__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ca40 + he4 --> ti44

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // chw0 
    ln_set_rate =  53.75 + -76.4273 * tfactors.T913i + 3.87451 * tfactors.T913
                         + -3.61477 * tfactors.T9 + 0.367451 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -76.4273 * tfactors.T943i + (1.0/3.0) * 3.87451 * tfactors.T923i
                                  + -3.61477 + (5.0/3.0) * 0.367451 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_p_sc43__ti44__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // sc43 + p --> ti44

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  36.8432 + -32.1734 * tfactors.T913i + -1.77078 * tfactors.T913
                         + -2.21706 * tfactors.T9 + 0.298499 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -32.1734 * tfactors.T943i + (1.0/3.0) * -1.77078 * tfactors.T923i
                                  + -2.21706 + (5.0/3.0) * 0.298499 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_p_sc43__he4_ca40__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // sc43 + p --> he4 + ca40

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  34.8559 + -32.1734 * tfactors.T913i + 0.0296879 * tfactors.T913
                         + -0.95232 * tfactors.T9 + 0.129022 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -32.1734 * tfactors.T943i + (1.0/3.0) * 0.0296879 * tfactors.T923i
                                  + -0.95232 + (5.0/3.0) * 0.129022 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_he4_ti44__cr48__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ti44 + he4 --> cr48

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  64.7958 + -81.667 * tfactors.T913i + -10.6333 * tfactors.T913
                         + -0.672613 * tfactors.T9 + 0.161209 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -81.667 * tfactors.T943i + (1.0/3.0) * -10.6333 * tfactors.T923i
                                  + -0.672613 + (5.0/3.0) * 0.161209 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_he4_ti44__p_v47__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ti44 + he4 --> p + v47

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // chw0r
    ln_set_rate =  -76.5154 + -10.7931 * tfactors.T9i + 70.2835 * tfactors.T913
                         + -7.99061 * tfactors.T9 + 0.486213 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  10.7931 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 70.2835 * tfactors.T923i
                                  + -7.99061 + (5.0/3.0) * 0.486213 * tfactors.T923 + -1.5 * tfactors.T9i;
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
void rate_p_v47__cr48__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // v47 + p --> cr48

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // nfisn
    ln_set_rate =  42.6798 + -6.0593 * tfactors.T9i + -34.0548 * tfactors.T913i + -3.41973 * tfactors.T913
                         + 1.16501 * tfactors.T9 + -0.105543 * tfactors.T953 + -7.70886 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  6.0593 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -34.0548 * tfactors.T943i + (1.0/3.0) * -3.41973 * tfactors.T923i
                                  + 1.16501 + (5.0/3.0) * -0.105543 * tfactors.T923 + -7.70886 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // nfisn
    ln_set_rate =  511.463 + -5.29491 * tfactors.T9i + 317.171 * tfactors.T913i + -911.679 * tfactors.T913
                         + 94.4245 * tfactors.T9 + -10.1973 * tfactors.T953 + 330.727 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  5.29491 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 317.171 * tfactors.T943i + (1.0/3.0) * -911.679 * tfactors.T923i
                                  + 94.4245 + (5.0/3.0) * -10.1973 * tfactors.T923 + 330.727 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // nfisn
    ln_set_rate =  23.8315 + 0.246665 * tfactors.T9i + -45.9868 * tfactors.T913i + 13.6822 * tfactors.T913
                         + -0.376902 * tfactors.T9 + -0.0194875 * tfactors.T953 + -8.42325 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  -0.246665 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -45.9868 * tfactors.T943i + (1.0/3.0) * 13.6822 * tfactors.T923i
                                  + -0.376902 + (5.0/3.0) * -0.0194875 * tfactors.T923 + -8.42325 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // nfisn
    ln_set_rate =  40.5626 + -0.514414 * tfactors.T9i + -110.655 * tfactors.T913i + 83.0232 * tfactors.T913
                         + -19.7762 * tfactors.T9 + 3.03961 * tfactors.T953 + -49.4742 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.514414 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -110.655 * tfactors.T943i + (1.0/3.0) * 83.0232 * tfactors.T923i
                                  + -19.7762 + (5.0/3.0) * 3.03961 * tfactors.T923 + -49.4742 * tfactors.T9i;
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
void rate_he4_cr48__fe52__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // cr48 + he4 --> fe52

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  65.1754 + -86.7459 * tfactors.T913i + -9.79373 * tfactors.T913
                         + -0.772169 * tfactors.T9 + 0.155883 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -86.7459 * tfactors.T943i + (1.0/3.0) * -9.79373 * tfactors.T923i
                                  + -0.772169 + (5.0/3.0) * 0.155883 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_he4_cr48__p_mn51__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // cr48 + he4 --> p + mn51

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  59.2276 + -86.7459 * tfactors.T913i + 1.05653 * tfactors.T913
                         + -1.15757 * tfactors.T9 + 0.0877546 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -86.7459 * tfactors.T943i + (1.0/3.0) * 1.05653 * tfactors.T923i
                                  + -1.15757 + (5.0/3.0) * 0.0877546 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_p_mn51__fe52__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // mn51 + p --> fe52

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  36.2596 + -36.1825 * tfactors.T913i + 0.873042 * tfactors.T913
                         + -2.89731 * tfactors.T9 + 0.364394 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -36.1825 * tfactors.T943i + (1.0/3.0) * 0.873042 * tfactors.T923i
                                  + -2.89731 + (5.0/3.0) * 0.364394 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_he4_fe52__ni56__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // fe52 + he4 --> ni56

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  66.6417 + -91.6819 * tfactors.T913i + -9.51885 * tfactors.T913
                         + -0.533014 * tfactors.T9 + 0.0892607 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -91.6819 * tfactors.T943i + (1.0/3.0) * -9.51885 * tfactors.T923i
                                  + -0.533014 + (5.0/3.0) * 0.0892607 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_he4_fe52__p_co55__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // fe52 + he4 --> p + co55

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  62.2207 + -91.6819 * tfactors.T913i + -0.329235 * tfactors.T913
                         + -0.780924 * tfactors.T9 + 0.0425179 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -91.6819 * tfactors.T943i + (1.0/3.0) * -0.329235 * tfactors.T923i
                                  + -0.780924 + (5.0/3.0) * 0.0425179 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_p_co55__ni56__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // co55 + p --> ni56

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  37.3736 + -38.1053 * tfactors.T913i + -0.210947 * tfactors.T913
                         + -2.68377 * tfactors.T9 + 0.355814 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -38.1053 * tfactors.T943i + (1.0/3.0) * -0.210947 * tfactors.T923i
                                  + -2.68377 + (5.0/3.0) * 0.355814 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_n13__p_c12__derived(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // n13 --> p + c12

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ls09c
    ln_set_rate =  40.040745583144826 + -22.547578710551672 * tfactors.T9i + -13.692 * tfactors.T913i + -0.230881 * tfactors.T913
                         + 4.44362 * tfactors.T9 + -3.15898 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  22.547578710551672 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -13.692 * tfactors.T943i + (1.0/3.0) * -0.230881 * tfactors.T923i
                                  + 4.44362 + (5.0/3.0) * -3.15898 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // ls09c
    ln_set_rate =  40.43534558314482 + -26.326068710551674 * tfactors.T9i + -5.10735 * tfactors.T913i + -2.24111 * tfactors.T913
                         + 0.148883 * tfactors.T9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  26.326068710551674 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -5.10735 * tfactors.T943i + (1.0/3.0) * -2.24111 * tfactors.T923i
                                  + 0.148883;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real c12_pf, dc12_pf_dT;
    // setting c12 partition function to 1.0 by default, independent of T
    c12_pf = 1.0_rt;
    dc12_pf_dT = 0.0_rt;

    Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    Real n13_pf, dn13_pf_dT;
    // setting n13 partition function to 1.0 by default, independent of T
    n13_pf = 1.0_rt;
    dn13_pf_dT = 0.0_rt;

    Real z_r = p_pf * c12_pf;
    Real z_p = n13_pf;

    Real dz_r_dT = c12_pf * dp_pf_dT + p_pf * dc12_pf_dT;
    Real dz_p_dT = dn13_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_n14__n_n13__derived(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // n14 --> n + n13

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // wiesc
    ln_set_rate =  19.55838975887632 + -125.4741960641516 * tfactors.T9i + 9.44873e-10 * tfactors.T913i + -2.33713e-09 * tfactors.T913
                         + 1.97507e-10 * tfactors.T9 + -1.49747e-11 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  125.4741960641516 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 9.44873e-10 * tfactors.T943i + (1.0/3.0) * -2.33713e-09 * tfactors.T923i
                                  + 1.97507e-10 + (5.0/3.0) * -1.49747e-11 * tfactors.T923;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // wiesc
    ln_set_rate =  37.12682975887632 + -122.4841912641516 * tfactors.T9i + 1.72241e-10 * tfactors.T913i + -5.62522e-10 * tfactors.T913
                         + 5.59212e-11 * tfactors.T9 + -4.6549e-12 * tfactors.T953 + 1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  122.4841912641516 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 1.72241e-10 * tfactors.T943i + (1.0/3.0) * -5.62522e-10 * tfactors.T923i
                                  + 5.59212e-11 + (5.0/3.0) * -4.6549e-12 * tfactors.T923 + 1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    Real n13_pf, dn13_pf_dT;
    // setting n13 partition function to 1.0 by default, independent of T
    n13_pf = 1.0_rt;
    dn13_pf_dT = 0.0_rt;

    Real n14_pf, dn14_pf_dT;
    // setting n14 partition function to 1.0 by default, independent of T
    n14_pf = 1.0_rt;
    dn14_pf_dT = 0.0_rt;

    Real z_r = n_pf * n13_pf;
    Real z_p = n14_pf;

    Real dz_r_dT = n13_pf * dn_pf_dT + n_pf * dn13_pf_dT;
    Real dz_p_dT = dn14_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_o16__he4_c12__derived(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // o16 --> he4 + c12

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // nac2 
    ln_set_rate =  279.2945280776573 + -84.95160042649215 * tfactors.T9i + 103.411 * tfactors.T913i + -420.567 * tfactors.T913
                         + 64.0874 * tfactors.T9 + -12.4624 * tfactors.T953 + 138.803 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  84.95160042649215 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 103.411 * tfactors.T943i + (1.0/3.0) * -420.567 * tfactors.T923i
                                  + 64.0874 + (5.0/3.0) * -12.4624 * tfactors.T923 + 138.803 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // nac2 
    ln_set_rate =  94.3131280776573 + -84.50317042649215 * tfactors.T9i + 58.9128 * tfactors.T913i + -148.273 * tfactors.T913
                         + 9.08324 * tfactors.T9 + -0.541041 * tfactors.T953 + 71.8554 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  84.50317042649215 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 58.9128 * tfactors.T943i + (1.0/3.0) * -148.273 * tfactors.T923i
                                  + 9.08324 + (5.0/3.0) * -0.541041 * tfactors.T923 + 71.8554 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real c12_pf, dc12_pf_dT;
    // setting c12 partition function to 1.0 by default, independent of T
    c12_pf = 1.0_rt;
    dc12_pf_dT = 0.0_rt;

    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real o16_pf, do16_pf_dT;
    // interpolating o16 partition function
    get_partition_function(O16, tfactors, o16_pf, do16_pf_dT);

    Real z_r = he4_pf * c12_pf;
    Real z_p = o16_pf;

    Real dz_r_dT = c12_pf * dhe4_pf_dT + he4_pf * dc12_pf_dT;
    Real dz_p_dT = do16_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_f18__he4_n14__derived(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // f18 --> he4 + n14

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10c
    ln_set_rate =  46.24897954391359 + -51.22930570025806 * tfactors.T9i + -36.2504 * tfactors.T913i
                         + -5.0 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  51.22930570025806 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -36.2504 * tfactors.T943i
                                  + (5.0/3.0) * -5.0 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  38.61457954391359 + -62.19490570025806 * tfactors.T9i + -5.6227 * tfactors.T913i;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  62.19490570025806 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -5.6227 * tfactors.T943i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  24.911917543913596 + -56.38964570025806 * tfactors.T9i;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  56.38964570025806 * tfactors.T9i * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real f18_pf, df18_pf_dT;
    // interpolating f18 partition function
    get_partition_function(F18, tfactors, f18_pf, df18_pf_dT);

    Real n14_pf, dn14_pf_dT;
    // setting n14 partition function to 1.0 by default, independent of T
    n14_pf = 1.0_rt;
    dn14_pf_dT = 0.0_rt;

    Real z_r = he4_pf * n14_pf;
    Real z_p = f18_pf;

    Real dz_r_dT = n14_pf * dhe4_pf_dT + he4_pf * dn14_pf_dT;
    Real dz_p_dT = df18_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_ne20__he4_o16__derived(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ne20 --> he4 + o16

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // co10c
    ln_set_rate =  34.265815859363634 + -67.65193003813836 * tfactors.T9i + -3.65925 * tfactors.T913
                         + 0.714224 * tfactors.T9 + -0.00107508 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  67.65193003813836 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -3.65925 * tfactors.T923i
                                  + 0.714224 + (5.0/3.0) * -0.00107508 * tfactors.T923;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // co10c
    ln_set_rate =  28.643045859363635 + -65.24613003813835 * tfactors.T9i;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  65.24613003813835 * tfactors.T9i * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // co10c
    ln_set_rate =  48.660335859363634 + -54.887630038138354 * tfactors.T9i + -39.7262 * tfactors.T913i + -0.210799 * tfactors.T913
                         + 0.442879 * tfactors.T9 + -0.0797753 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  54.887630038138354 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -39.7262 * tfactors.T943i + (1.0/3.0) * -0.210799 * tfactors.T923i
                                  + 0.442879 + (5.0/3.0) * -0.0797753 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real o16_pf, do16_pf_dT;
    // interpolating o16 partition function
    get_partition_function(O16, tfactors, o16_pf, do16_pf_dT);

    Real ne20_pf, dne20_pf_dT;
    // interpolating ne20 partition function
    get_partition_function(Ne20, tfactors, ne20_pf, dne20_pf_dT);

    Real z_r = he4_pf * o16_pf;
    Real z_p = ne20_pf;

    Real dz_r_dT = o16_pf * dhe4_pf_dT + he4_pf * do16_pf_dT;
    Real dz_p_dT = dne20_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_ne21__n_ne20__derived(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ne21 --> n + ne20

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ka02c
    ln_set_rate =  34.980677217841034 + -80.1620770211219 * tfactors.T9i;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  80.1620770211219 * tfactors.T9i * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // ka02c
    ln_set_rate =  30.822787217841032 + -78.4581470211219 * tfactors.T9i
                         + 1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  78.4581470211219 * tfactors.T9i * tfactors.T9i
                                  + 1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    Real ne21_pf, dne21_pf_dT;
    // interpolating ne21 partition function
    get_partition_function(Ne21, tfactors, ne21_pf, dne21_pf_dT);

    Real ne20_pf, dne20_pf_dT;
    // interpolating ne20 partition function
    get_partition_function(Ne20, tfactors, ne20_pf, dne20_pf_dT);

    Real z_r = n_pf * ne20_pf;
    Real z_p = ne21_pf;

    Real dz_r_dT = ne20_pf * dn_pf_dT + n_pf * dne20_pf_dT;
    Real dz_p_dT = dne21_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_na22__p_ne21__derived(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // na22 --> p + ne21

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10c
    ln_set_rate =  -24.57903898617269 + -78.40599033331655 * tfactors.T9i;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  78.40599033331655 * tfactors.T9i * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  42.14596101382731 + -78.20981033331655 * tfactors.T9i + -19.2096 * tfactors.T913i
                         + -1.0 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  78.20981033331655 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -19.2096 * tfactors.T943i
                                  + (5.0/3.0) * -1.0 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  -16.409838986172687 + -82.42366033331655 * tfactors.T9i + 21.1176 * tfactors.T913i + 34.0411 * tfactors.T913
                         + -4.45593 * tfactors.T9 + 0.328613 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  82.42366033331655 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 21.1176 * tfactors.T943i + (1.0/3.0) * 34.0411 * tfactors.T923i
                                  + -4.45593 + (5.0/3.0) * 0.328613 * tfactors.T923;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  24.83340101382731 + -79.60938033331655 * tfactors.T9i;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  79.60938033331655 * tfactors.T9i * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real ne21_pf, dne21_pf_dT;
    // interpolating ne21 partition function
    get_partition_function(Ne21, tfactors, ne21_pf, dne21_pf_dT);

    Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    Real na22_pf, dna22_pf_dT;
    // interpolating na22 partition function
    get_partition_function(Na22, tfactors, na22_pf, dna22_pf_dT);

    Real z_r = p_pf * ne21_pf;
    Real z_p = na22_pf;

    Real dz_r_dT = ne21_pf * dp_pf_dT + p_pf * dne21_pf_dT;
    Real dz_p_dT = dna22_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_na22__he4_f18__derived(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // na22 --> he4 + f18

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // rpsmr
    ln_set_rate =  59.32234728275453 + -100.23588367240258 * tfactors.T9i + 18.8956 * tfactors.T913i + -65.6134 * tfactors.T913
                         + 1.71114 * tfactors.T9 + -0.0260999 * tfactors.T953 + 39.3396 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  100.23588367240258 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 18.8956 * tfactors.T943i + (1.0/3.0) * -65.6134 * tfactors.T923i
                                  + 1.71114 + (5.0/3.0) * -0.0260999 * tfactors.T923 + 39.3396 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real f18_pf, df18_pf_dT;
    // interpolating f18 partition function
    get_partition_function(F18, tfactors, f18_pf, df18_pf_dT);

    Real na22_pf, dna22_pf_dT;
    // interpolating na22 partition function
    get_partition_function(Na22, tfactors, na22_pf, dna22_pf_dT);

    Real z_r = he4_pf * f18_pf;
    Real z_p = na22_pf;

    Real dz_r_dT = f18_pf * dhe4_pf_dT + he4_pf * df18_pf_dT;
    Real dz_p_dT = dna22_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_na23__n_na22__derived(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // na23 --> n + na22

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  37.06649496929424 + -144.1130291985219 * tfactors.T9i + 1.02148 * tfactors.T913
                         + -0.334638 * tfactors.T9 + 0.0258708 * tfactors.T953 + 1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  144.1130291985219 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 1.02148 * tfactors.T923i
                                  + -0.334638 + (5.0/3.0) * 0.0258708 * tfactors.T923 + 1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    Real na23_pf, dna23_pf_dT;
    // interpolating na23 partition function
    get_partition_function(Na23, tfactors, na23_pf, dna23_pf_dT);

    Real na22_pf, dna22_pf_dT;
    // interpolating na22 partition function
    get_partition_function(Na22, tfactors, na22_pf, dna22_pf_dT);

    Real z_r = n_pf * na22_pf;
    Real z_p = na23_pf;

    Real dz_r_dT = na22_pf * dn_pf_dT + n_pf * dna22_pf_dT;
    Real dz_p_dT = dna23_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_mg24__p_na23__derived(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // mg24 --> p + na23

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10c
    ln_set_rate =  43.93571176470677 + -135.6881490421346 * tfactors.T9i + -20.6428 * tfactors.T913i + 1.52954 * tfactors.T913
                         + 2.7487 * tfactors.T9 + -1.0 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  135.6881490421346 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -20.6428 * tfactors.T943i + (1.0/3.0) * 1.52954 * tfactors.T923i
                                  + 2.7487 + (5.0/3.0) * -1.0 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  34.087611764706764 + -138.9684390421346 * tfactors.T9i + -0.360588 * tfactors.T913
                         + 1.4187 * tfactors.T9 + -0.184061 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  138.9684390421346 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.360588 * tfactors.T923i
                                  + 1.4187 + (5.0/3.0) * -0.184061 * tfactors.T923;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  20.00236176470677 + -137.3003390421346 * tfactors.T9i;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  137.3003390421346 * tfactors.T9i * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real na23_pf, dna23_pf_dT;
    // interpolating na23 partition function
    get_partition_function(Na23, tfactors, na23_pf, dna23_pf_dT);

    Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    Real mg24_pf, dmg24_pf_dT;
    // interpolating mg24 partition function
    get_partition_function(Mg24, tfactors, mg24_pf, dmg24_pf_dT);

    Real z_r = p_pf * na23_pf;
    Real z_p = mg24_pf;

    Real dz_r_dT = na23_pf * dp_pf_dT + p_pf * dna23_pf_dT;
    Real dz_p_dT = dmg24_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_mg24__he4_ne20__derived(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // mg24 --> he4 + ne20

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10c
    ln_set_rate =  -13.886931148855972 + -110.62012330714883 * tfactors.T9i;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  110.62012330714883 * tfactors.T9i * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  49.32436885114403 + -108.11407330714883 * tfactors.T9i + -46.2525 * tfactors.T913i + 5.58901 * tfactors.T913
                         + 7.61843 * tfactors.T9 + -3.683 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  108.11407330714883 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -46.2525 * tfactors.T943i + (1.0/3.0) * 5.58901 * tfactors.T923i
                                  + 7.61843 + (5.0/3.0) * -3.683 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  16.020298851144027 + -120.89497330714883 * tfactors.T9i + 16.9229 * tfactors.T913
                         + -2.57325 * tfactors.T9 + 0.208997 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  120.89497330714883 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 16.9229 * tfactors.T923i
                                  + -2.57325 + (5.0/3.0) * 0.208997 * tfactors.T923;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  26.80163885114403 + -117.33433330714882 * tfactors.T9i;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  117.33433330714882 * tfactors.T9i * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real mg24_pf, dmg24_pf_dT;
    // interpolating mg24 partition function
    get_partition_function(Mg24, tfactors, mg24_pf, dmg24_pf_dT);

    Real ne20_pf, dne20_pf_dT;
    // interpolating ne20 partition function
    get_partition_function(Ne20, tfactors, ne20_pf, dne20_pf_dT);

    Real z_r = he4_pf * ne20_pf;
    Real z_p = mg24_pf;

    Real dz_r_dT = ne20_pf * dhe4_pf_dT + he4_pf * dne20_pf_dT;
    Real dz_p_dT = dmg24_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_si28__p_al27__derived(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // si28 --> p + al27

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10c
    ln_set_rate =  11.77656482818681 + -136.3434628922348 * tfactors.T9i + 23.8634 * tfactors.T913
                         + -3.70135 * tfactors.T9 + 0.28964 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  136.3434628922348 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 23.8634 * tfactors.T923i
                                  + -3.70135 + (5.0/3.0) * 0.28964 * tfactors.T923;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  111.46636482818681 + -134.8268158922348 * tfactors.T9i + -26.8327 * tfactors.T913i + -116.137 * tfactors.T913
                         + 0.00950567 * tfactors.T9 + 0.00999755 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  134.8268158922348 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -26.8327 * tfactors.T943i + (1.0/3.0) * -116.137 * tfactors.T923i
                                  + 0.00950567 + (5.0/3.0) * 0.00999755 * tfactors.T923;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  46.549464828186814 + -134.4395028922348 * tfactors.T9i + -23.2205 * tfactors.T913i
                         + -2.0 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  134.4395028922348 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -23.2205 * tfactors.T943i
                                  + (5.0/3.0) * -2.0 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    Real al27_pf, dal27_pf_dT;
    // interpolating al27 partition function
    get_partition_function(Al27, tfactors, al27_pf, dal27_pf_dT);

    Real si28_pf, dsi28_pf_dT;
    // interpolating si28 partition function
    get_partition_function(Si28, tfactors, si28_pf, dsi28_pf_dT);

    Real z_r = p_pf * al27_pf;
    Real z_p = si28_pf;

    Real dz_r_dT = al27_pf * dp_pf_dT + p_pf * dal27_pf_dT;
    Real dz_p_dT = dsi28_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_si28__he4_mg24__derived(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // si28 --> he4 + mg24

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // st08r
    ln_set_rate =  -25.68857483340593 + -128.69270892750794 * tfactors.T9i + 21.3721 * tfactors.T913i + 37.7649 * tfactors.T913
                         + -4.10635 * tfactors.T9 + 0.249618 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  128.69270892750794 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 21.3721 * tfactors.T943i + (1.0/3.0) * 37.7649 * tfactors.T923i
                                  + -4.10635 + (5.0/3.0) * 0.249618 * tfactors.T923;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // st08r
    ln_set_rate =  32.90059516659407 + -131.48850892750792 * tfactors.T9i;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  131.48850892750792 * tfactors.T9i * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real mg24_pf, dmg24_pf_dT;
    // interpolating mg24 partition function
    get_partition_function(Mg24, tfactors, mg24_pf, dmg24_pf_dT);

    Real si28_pf, dsi28_pf_dT;
    // interpolating si28 partition function
    get_partition_function(Si28, tfactors, si28_pf, dsi28_pf_dT);

    Real z_r = he4_pf * mg24_pf;
    Real z_p = si28_pf;

    Real dz_r_dT = mg24_pf * dhe4_pf_dT + he4_pf * dmg24_pf_dT;
    Real dz_p_dT = dsi28_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_s32__p_p31__derived(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // s32 --> p + p31

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10c
    ln_set_rate =  25.172836958303144 + -106.63693563716609 * tfactors.T9i + 8.09341 * tfactors.T913
                         + -0.615971 * tfactors.T9 + 0.031159 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  106.63693563716609 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 8.09341 * tfactors.T923i
                                  + -0.615971 + (5.0/3.0) * 0.031159 * tfactors.T923;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  21.682890958303144 + -105.11947563716609 * tfactors.T9i;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  105.11947563716609 * tfactors.T9i * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  43.61088095830314 + -102.85989563716609 * tfactors.T9i + -25.3278 * tfactors.T913i + 6.4931 * tfactors.T913
                         + -9.27513 * tfactors.T9 + -0.610439 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  102.85989563716609 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -25.3278 * tfactors.T943i + (1.0/3.0) * 6.4931 * tfactors.T923i
                                  + -9.27513 + (5.0/3.0) * -0.610439 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    Real s32_pf, ds32_pf_dT;
    // interpolating s32 partition function
    get_partition_function(S32, tfactors, s32_pf, ds32_pf_dT);

    Real p31_pf, dp31_pf_dT;
    // interpolating p31 partition function
    get_partition_function(P31, tfactors, p31_pf, dp31_pf_dT);

    Real z_r = p_pf * p31_pf;
    Real z_p = s32_pf;

    Real dz_r_dT = p31_pf * dp_pf_dT + p_pf * dp31_pf_dT;
    Real dz_p_dT = ds32_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_s32__he4_si28__derived(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // s32 --> he4 + si28

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  72.81295409739818 + -80.62610309662641 * tfactors.T9i + -59.4896 * tfactors.T913i + 4.47205 * tfactors.T913
                         + -4.78989 * tfactors.T9 + 0.557201 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  80.62610309662641 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -59.4896 * tfactors.T943i + (1.0/3.0) * 4.47205 * tfactors.T923i
                                  + -4.78989 + (5.0/3.0) * 0.557201 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real s32_pf, ds32_pf_dT;
    // interpolating s32 partition function
    get_partition_function(S32, tfactors, s32_pf, ds32_pf_dT);

    Real si28_pf, dsi28_pf_dT;
    // interpolating si28 partition function
    get_partition_function(Si28, tfactors, si28_pf, dsi28_pf_dT);

    Real z_r = he4_pf * si28_pf;
    Real z_p = s32_pf;

    Real dz_r_dT = si28_pf * dhe4_pf_dT + he4_pf * dsi28_pf_dT;
    Real dz_p_dT = ds32_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_c12__he4_he4_he4__derived(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // c12 --> he4 + he4 + he4

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // fy05c
    ln_set_rate =  22.39392447043974 + -88.54942933569913 * tfactors.T9i + -13.49 * tfactors.T913i + 21.4259 * tfactors.T913
                         + -1.34769 * tfactors.T9 + 0.0879816 * tfactors.T953 + -10.1653 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  88.54942933569913 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -13.49 * tfactors.T943i + (1.0/3.0) * 21.4259 * tfactors.T923i
                                  + -1.34769 + (5.0/3.0) * 0.0879816 * tfactors.T923 + -10.1653 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // fy05c
    ln_set_rate =  34.95602447043974 + -85.44732933569914 * tfactors.T9i + -23.57 * tfactors.T913i + 20.4886 * tfactors.T913
                         + -12.9882 * tfactors.T9 + -20.0 * tfactors.T953 + 0.8333300000000001 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  85.44732933569914 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -23.57 * tfactors.T943i + (1.0/3.0) * 20.4886 * tfactors.T923i
                                  + -12.9882 + (5.0/3.0) * -20.0 * tfactors.T923 + 0.8333300000000001 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // fy05c
    ln_set_rate =  45.77337247043974 + -84.42286933569913 * tfactors.T9i + -37.06 * tfactors.T913i + 29.3493 * tfactors.T913
                         + -115.507 * tfactors.T9 + -10.0 * tfactors.T953 + 1.66667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  84.42286933569913 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -37.06 * tfactors.T943i + (1.0/3.0) * 29.3493 * tfactors.T923i
                                  + -115.507 + (5.0/3.0) * -10.0 * tfactors.T923 + 1.66667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real c12_pf, dc12_pf_dT;
    // setting c12 partition function to 1.0 by default, independent of T
    c12_pf = 1.0_rt;
    dc12_pf_dT = 0.0_rt;

    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real z_r = he4_pf * he4_pf * he4_pf;
    Real z_p = c12_pf;

    Real dz_r_dT = dhe4_pf_dT + dhe4_pf_dT + dhe4_pf_dT;
    Real dz_p_dT = dc12_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_o16__he4_n13__derived(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // o16 + p --> he4 + n13

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // cf88n
    ln_set_rate =  42.232382494512464 + -60.55237555926846 * tfactors.T9i + -35.829 * tfactors.T913i + -0.530275 * tfactors.T913
                         + -0.982462 * tfactors.T9 + 0.0808059 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  60.55237555926846 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -35.829 * tfactors.T943i + (1.0/3.0) * -0.530275 * tfactors.T923i
                                  + -0.982462 + (5.0/3.0) * 0.0808059 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    Real o16_pf, do16_pf_dT;
    // interpolating o16 partition function
    get_partition_function(O16, tfactors, o16_pf, do16_pf_dT);

    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real n13_pf, dn13_pf_dT;
    // setting n13 partition function to 1.0 by default, independent of T
    n13_pf = 1.0_rt;
    dn13_pf_dT = 0.0_rt;

    Real z_r = he4_pf * n13_pf;
    Real z_p = p_pf * o16_pf;

    Real dz_r_dT = n13_pf * dhe4_pf_dT + he4_pf * dn13_pf_dT;
    Real dz_p_dT = o16_pf * dp_pf_dT + p_pf * do16_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_ne20__p_na23__derived(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ne20 + he4 --> p + na23

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10c
    ln_set_rate =  -6.377717086437262 + -29.88961364462334 * tfactors.T9i + 19.7297 * tfactors.T913
                         + -2.20987 * tfactors.T9 + 0.153374 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  29.88961364462334 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 19.7297 * tfactors.T923i
                                  + -2.20987 + (5.0/3.0) * 0.153374 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  0.227472413562738 + -29.434873644623337 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  29.434873644623337 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  19.18524291356274 + -27.573843644623338 * tfactors.T9i + -20.0024 * tfactors.T913i + 11.5988 * tfactors.T913
                         + -1.37398 * tfactors.T9 + -1.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  27.573843644623338 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -20.0024 * tfactors.T943i + (1.0/3.0) * 11.5988 * tfactors.T923i
                                  + -1.37398 + (5.0/3.0) * -1.0 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real na23_pf, dna23_pf_dT;
    // interpolating na23 partition function
    get_partition_function(Na23, tfactors, na23_pf, dna23_pf_dT);

    Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real ne20_pf, dne20_pf_dT;
    // interpolating ne20 partition function
    get_partition_function(Ne20, tfactors, ne20_pf, dne20_pf_dT);

    Real z_r = p_pf * na23_pf;
    Real z_p = he4_pf * ne20_pf;

    Real dz_r_dT = na23_pf * dp_pf_dT + p_pf * dna23_pf_dT;
    Real dz_p_dT = ne20_pf * dhe4_pf_dT + he4_pf * dne20_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_ne21__he4_f18__derived(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ne21 + p --> he4 + f18

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // rpsmr
    ln_set_rate =  50.65368626892722 + -22.049056049959066 * tfactors.T9i + 21.4461 * tfactors.T913i + -73.252 * tfactors.T913
                         + 2.42329 * tfactors.T9 + -0.077278 * tfactors.T953 + 40.7604 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  22.049056049959066 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 21.4461 * tfactors.T943i + (1.0/3.0) * -73.252 * tfactors.T923i
                                  + 2.42329 + (5.0/3.0) * -0.077278 * tfactors.T923 + 40.7604 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    Real ne21_pf, dne21_pf_dT;
    // interpolating ne21 partition function
    get_partition_function(Ne21, tfactors, ne21_pf, dne21_pf_dT);

    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real f18_pf, df18_pf_dT;
    // interpolating f18 partition function
    get_partition_function(F18, tfactors, f18_pf, df18_pf_dT);

    Real z_r = he4_pf * f18_pf;
    Real z_p = p_pf * ne21_pf;

    Real dz_r_dT = f18_pf * dhe4_pf_dT + he4_pf * df18_pf_dT;
    Real dz_p_dT = ne21_pf * dp_pf_dT + p_pf * dne21_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_n_mg24__he4_ne21__derived(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // mg24 + n --> he4 + ne21

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // nacrc
    ln_set_rate =  46.548491633302994 + -29.649543801059973 * tfactors.T9i + -46.88 * tfactors.T913i + -0.536629 * tfactors.T913
                         + 0.144715 * tfactors.T9 + -0.197624 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  29.649543801059973 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -46.88 * tfactors.T943i + (1.0/3.0) * -0.536629 * tfactors.T923i
                                  + 0.144715 + (5.0/3.0) * -0.197624 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // nacrc
    ln_set_rate =  -4.696018366697003 + -42.91334380105997 * tfactors.T9i + 18.0748 * tfactors.T913
                         + -0.981883 * tfactors.T9 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  42.91334380105997 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 18.0748 * tfactors.T923i
                                  + -0.981883 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    Real ne21_pf, dne21_pf_dT;
    // interpolating ne21 partition function
    get_partition_function(Ne21, tfactors, ne21_pf, dne21_pf_dT);

    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real mg24_pf, dmg24_pf_dT;
    // interpolating mg24 partition function
    get_partition_function(Mg24, tfactors, mg24_pf, dmg24_pf_dT);

    Real z_r = he4_pf * ne21_pf;
    Real z_p = n_pf * mg24_pf;

    Real dz_r_dT = ne21_pf * dhe4_pf_dT + he4_pf * dne21_pf_dT;
    Real dz_p_dT = mg24_pf * dn_pf_dT + n_pf * dmg24_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_mg24__p_al27__derived(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // mg24 + he4 --> p + al27

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10c
    ln_set_rate =  -6.445750338407262 + -22.82086933218981 * tfactors.T9i + 18.0416 * tfactors.T913
                         + -1.54137 * tfactors.T9 + 0.0847506 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  22.82086933218981 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 18.0416 * tfactors.T923i
                                  + -1.54137 + (5.0/3.0) * 0.0847506 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  -26.286160338407264 + -19.54138133218981 * tfactors.T9i + 5.18642 * tfactors.T913i + -34.7936 * tfactors.T913
                         + 168.225 * tfactors.T9 + -115.825 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  19.54138133218981 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 5.18642 * tfactors.T943i + (1.0/3.0) * -34.7936 * tfactors.T923i
                                  + 168.225 + (5.0/3.0) * -115.825 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  30.039739661592737 + -18.578369332189812 * tfactors.T9i + -26.4162 * tfactors.T913i
                         + -2.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  18.578369332189812 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -26.4162 * tfactors.T943i
                                  + (5.0/3.0) * -2.0 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    Real al27_pf, dal27_pf_dT;
    // interpolating al27 partition function
    get_partition_function(Al27, tfactors, al27_pf, dal27_pf_dT);

    Real mg24_pf, dmg24_pf_dT;
    // interpolating mg24 partition function
    get_partition_function(Mg24, tfactors, mg24_pf, dmg24_pf_dT);

    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real z_r = p_pf * al27_pf;
    Real z_p = he4_pf * mg24_pf;

    Real dz_r_dT = al27_pf * dp_pf_dT + p_pf * dal27_pf_dT;
    Real dz_p_dT = mg24_pf * dhe4_pf_dT + he4_pf * dmg24_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_si28__p_p31__derived(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // si28 + he4 --> p + p31

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10c
    ln_set_rate =  -11.433473139095032 + -25.65965858572089 * tfactors.T9i + 21.521 * tfactors.T913
                         + -1.90355 * tfactors.T9 + 0.092724 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  25.65965858572089 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 21.521 * tfactors.T923i
                                  + -1.90355 + (5.0/3.0) * 0.092724 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  -13.459473139095032 + -24.11106858572089 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  24.11106858572089 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  60.342426860904965 + -22.23390858572089 * tfactors.T9i + -31.932 * tfactors.T913i + -77.0334 * tfactors.T913
                         + -43.6847 * tfactors.T9 + -4.28955 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  22.23390858572089 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -31.932 * tfactors.T943i + (1.0/3.0) * -77.0334 * tfactors.T923i
                                  + -43.6847 + (5.0/3.0) * -4.28955 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real p31_pf, dp31_pf_dT;
    // interpolating p31 partition function
    get_partition_function(P31, tfactors, p31_pf, dp31_pf_dT);

    Real si28_pf, dsi28_pf_dT;
    // interpolating si28 partition function
    get_partition_function(Si28, tfactors, si28_pf, dsi28_pf_dT);

    Real z_r = p_pf * p31_pf;
    Real z_p = he4_pf * si28_pf;

    Real dz_r_dT = p31_pf * dp_pf_dT + p_pf * dp31_pf_dT;
    Real dz_p_dT = si28_pf * dhe4_pf_dT + he4_pf * dsi28_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_s32__p_cl35__derived__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // s32 + he4 --> p + cl35

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10c
    ln_set_rate =  2.425628238034531 + -27.66622776398283 * tfactors.T9i + 5.33756 * tfactors.T913
                         + 1.64418 * tfactors.T9 + -0.246167 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  27.66622776398283 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 5.33756 * tfactors.T923i
                                  + 1.64418 + (5.0/3.0) * -0.246167 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  -0.877601761965469 + -25.59141776398283 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  25.59141776398283 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  -57.39498176196547 + -22.18939876398283 * tfactors.T9i + 25.5338 * tfactors.T913
                         + 6.45824 * tfactors.T9 + -0.950294 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  22.18939876398283 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 25.5338 * tfactors.T923i
                                  + 6.45824 + (5.0/3.0) * -0.950294 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  32.25441823803453 + -21.65646776398283 * tfactors.T9i + -30.9147 * tfactors.T913i + -1.2345 * tfactors.T913
                         + 22.5118 * tfactors.T9 + -33.0589 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  21.65646776398283 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -30.9147 * tfactors.T943i + (1.0/3.0) * -1.2345 * tfactors.T923i
                                  + 22.5118 + (5.0/3.0) * -33.0589 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    Real cl35_pf, dcl35_pf_dT;
    // interpolating cl35 partition function
    get_partition_function(Cl35, tfactors, cl35_pf, dcl35_pf_dT);

    Real s32_pf, ds32_pf_dT;
    // interpolating s32 partition function
    get_partition_function(S32, tfactors, s32_pf, ds32_pf_dT);

    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real z_r = p_pf * cl35_pf;
    Real z_p = he4_pf * s32_pf;

    Real dz_r_dT = cl35_pf * dp_pf_dT + p_pf * dcl35_pf_dT;
    Real dz_p_dT = s32_pf * dhe4_pf_dT + he4_pf * ds32_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_ar36__he4_s32__derived__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ar36 --> he4 + s32

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  73.81637663285039 + -77.06281976216322 * tfactors.T9i + -65.3709 * tfactors.T913i + 5.68294 * tfactors.T913
                         + -5.00388 * tfactors.T9 + 0.571407 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  77.06281976216322 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -65.3709 * tfactors.T943i + (1.0/3.0) * 5.68294 * tfactors.T923i
                                  + -5.00388 + (5.0/3.0) * 0.571407 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real s32_pf, ds32_pf_dT;
    // interpolating s32 partition function
    get_partition_function(S32, tfactors, s32_pf, ds32_pf_dT);

    Real ar36_pf, dar36_pf_dT;
    // interpolating ar36 partition function
    get_partition_function(Ar36, tfactors, ar36_pf, dar36_pf_dT);

    Real z_r = he4_pf * s32_pf;
    Real z_p = ar36_pf;

    Real dz_r_dT = s32_pf * dhe4_pf_dT + he4_pf * ds32_pf_dT;
    Real dz_p_dT = dar36_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_ar36__p_cl35__derived__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ar36 --> p + cl35

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // il10c
    ln_set_rate =  16.016854870884917 + -100.72924752614605 * tfactors.T9i;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  100.72924752614605 * tfactors.T9i * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  -17.475105129115086 + -99.28393852614604 * tfactors.T9i;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  99.28393852614604 * tfactors.T9i * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  60.736594870884915 + -98.71928752614605 * tfactors.T9i + -27.8971 * tfactors.T913i + -16.2304 * tfactors.T913
                         + 35.255 * tfactors.T9 + -25.8411 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  98.71928752614605 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -27.8971 * tfactors.T943i + (1.0/3.0) * -16.2304 * tfactors.T923i
                                  + 35.255 + (5.0/3.0) * -25.8411 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // il10c
    ln_set_rate =  17.20280487088492 + -102.37020752614605 * tfactors.T9i + 18.0179 * tfactors.T913
                         + -2.86304 * tfactors.T9 + 0.250854 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  102.37020752614605 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 18.0179 * tfactors.T923i
                                  + -2.86304 + (5.0/3.0) * 0.250854 * tfactors.T923;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    Real cl35_pf, dcl35_pf_dT;
    // interpolating cl35 partition function
    get_partition_function(Cl35, tfactors, cl35_pf, dcl35_pf_dT);

    Real ar36_pf, dar36_pf_dT;
    // interpolating ar36 partition function
    get_partition_function(Ar36, tfactors, ar36_pf, dar36_pf_dT);

    Real z_r = p_pf * cl35_pf;
    Real z_p = ar36_pf;

    Real dz_r_dT = cl35_pf * dp_pf_dT + p_pf * dcl35_pf_dT;
    Real dz_p_dT = dar36_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_ar36__p_k39__derived__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ar36 + he4 --> p + k39

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  20.636664061510302 + -14.953349961318928 * tfactors.T9i + -30.0732 * tfactors.T913i + 7.03263 * tfactors.T913
                         + -1.10085 * tfactors.T9 + 0.133768 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  14.953349961318928 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -30.0732 * tfactors.T943i + (1.0/3.0) * 7.03263 * tfactors.T923i
                                  + -1.10085 + (5.0/3.0) * 0.133768 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    Real k39_pf, dk39_pf_dT;
    // interpolating k39 partition function
    get_partition_function(K39, tfactors, k39_pf, dk39_pf_dT);

    Real ar36_pf, dar36_pf_dT;
    // interpolating ar36 partition function
    get_partition_function(Ar36, tfactors, ar36_pf, dar36_pf_dT);

    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real z_r = p_pf * k39_pf;
    Real z_p = he4_pf * ar36_pf;

    Real dz_r_dT = k39_pf * dp_pf_dT + p_pf * dk39_pf_dT;
    Real dz_p_dT = ar36_pf * dhe4_pf_dT + he4_pf * dar36_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_ca40__he4_ar36__derived__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ca40 --> he4 + ar36

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  77.28261041284821 + -81.69174599574632 * tfactors.T9i + -71.0046 * tfactors.T913i + 4.0656 * tfactors.T913
                         + -5.26509 * tfactors.T9 + 0.683546 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  81.69174599574632 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -71.0046 * tfactors.T943i + (1.0/3.0) * 4.0656 * tfactors.T923i
                                  + -5.26509 + (5.0/3.0) * 0.683546 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real ca40_pf, dca40_pf_dT;
    // interpolating ca40 partition function
    get_partition_function(Ca40, tfactors, ca40_pf, dca40_pf_dT);

    Real ar36_pf, dar36_pf_dT;
    // interpolating ar36 partition function
    get_partition_function(Ar36, tfactors, ar36_pf, dar36_pf_dT);

    Real z_r = he4_pf * ar36_pf;
    Real z_p = ca40_pf;

    Real dz_r_dT = ar36_pf * dhe4_pf_dT + he4_pf * dar36_pf_dT;
    Real dz_p_dT = dca40_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_ca40__p_k39__derived__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ca40 --> p + k39

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // lo18r
    ln_set_rate =  2786.4340744743586 + -101.86662364079672 * tfactors.T9i + 802.18 * tfactors.T913i + -4010.27 * tfactors.T913
                         + 1137.69 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  101.86662364079672 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 802.18 * tfactors.T943i + (1.0/3.0) * -4010.27 * tfactors.T923i
                                  + 1137.69 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // lo18r
    ln_set_rate =  613.1530744743586 + -109.20898364079672 * tfactors.T9i + 641.844 * tfactors.T913i + -1248.49 * tfactors.T913
                         + 566.426 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  109.20898364079672 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 641.844 * tfactors.T943i + (1.0/3.0) * -1248.49 * tfactors.T923i
                                  + 566.426 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // lo18r
    ln_set_rate =  127.30607447435852 + -98.30936364079672 * tfactors.T9i + 41.1723 * tfactors.T913i + -149.299 * tfactors.T913
                         + 10.5229 * tfactors.T9 + -0.68208 * tfactors.T953 + 60.7367 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  98.30936364079672 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 41.1723 * tfactors.T943i + (1.0/3.0) * -149.299 * tfactors.T923i
                                  + 10.5229 + (5.0/3.0) * -0.68208 * tfactors.T923 + 60.7367 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    Real ca40_pf, dca40_pf_dT;
    // interpolating ca40 partition function
    get_partition_function(Ca40, tfactors, ca40_pf, dca40_pf_dT);

    Real k39_pf, dk39_pf_dT;
    // interpolating k39 partition function
    get_partition_function(K39, tfactors, k39_pf, dk39_pf_dT);

    Real z_r = p_pf * k39_pf;
    Real z_p = ca40_pf;

    Real dz_r_dT = k39_pf * dp_pf_dT + p_pf * dk39_pf_dT;
    Real dz_p_dT = dca40_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_he4_ca40__p_sc43__derived__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ca40 + he4 --> p + sc43

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  35.65752817292938 + -40.87575463203665 * tfactors.T9i + -32.1734 * tfactors.T913i + 0.0296879 * tfactors.T913
                         + -0.95232 * tfactors.T9 + 0.129022 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  40.87575463203665 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -32.1734 * tfactors.T943i + (1.0/3.0) * 0.0296879 * tfactors.T923i
                                  + -0.95232 + (5.0/3.0) * 0.129022 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    Real ca40_pf, dca40_pf_dT;
    // interpolating ca40 partition function
    get_partition_function(Ca40, tfactors, ca40_pf, dca40_pf_dT);

    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real sc43_pf, dsc43_pf_dT;
    // interpolating sc43 partition function
    get_partition_function(Sc43, tfactors, sc43_pf, dsc43_pf_dT);

    Real z_r = p_pf * sc43_pf;
    Real z_p = he4_pf * ca40_pf;

    Real dz_r_dT = sc43_pf * dp_pf_dT + p_pf * dsc43_pf_dT;
    Real dz_p_dT = ca40_pf * dhe4_pf_dT + he4_pf * dca40_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_ti44__he4_ca40__derived__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ti44 --> he4 + ca40

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // chw0 
    ln_set_rate =  78.69908591662846 + -59.4975248620018 * tfactors.T9i + -76.4273 * tfactors.T913i + 3.87451 * tfactors.T913
                         + -3.61477 * tfactors.T9 + 0.367451 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  59.4975248620018 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -76.4273 * tfactors.T943i + (1.0/3.0) * 3.87451 * tfactors.T923i
                                  + -3.61477 + (5.0/3.0) * 0.367451 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real ca40_pf, dca40_pf_dT;
    // interpolating ca40 partition function
    get_partition_function(Ca40, tfactors, ca40_pf, dca40_pf_dT);

    Real ti44_pf, dti44_pf_dT;
    // interpolating ti44 partition function
    get_partition_function(Ti44, tfactors, ti44_pf, dti44_pf_dT);

    Real z_r = he4_pf * ca40_pf;
    Real z_p = ti44_pf;

    Real dz_r_dT = ca40_pf * dhe4_pf_dT + he4_pf * dca40_pf_dT;
    Real dz_p_dT = dti44_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_ti44__p_sc43__derived__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ti44 --> p + sc43

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  62.59391408955786 + -100.37327949403844 * tfactors.T9i + -32.1734 * tfactors.T913i + -1.77078 * tfactors.T913
                         + -2.21706 * tfactors.T9 + 0.298499 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  100.37327949403844 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -32.1734 * tfactors.T943i + (1.0/3.0) * -1.77078 * tfactors.T923i
                                  + -2.21706 + (5.0/3.0) * 0.298499 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    Real ti44_pf, dti44_pf_dT;
    // interpolating ti44 partition function
    get_partition_function(Ti44, tfactors, ti44_pf, dti44_pf_dT);

    Real sc43_pf, dsc43_pf_dT;
    // interpolating sc43 partition function
    get_partition_function(Sc43, tfactors, sc43_pf, dsc43_pf_dT);

    Real z_r = p_pf * sc43_pf;
    Real z_p = ti44_pf;

    Real dz_r_dT = sc43_pf * dp_pf_dT + p_pf * dsc43_pf_dT;
    Real dz_p_dT = dti44_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_cr48__he4_ti44__derived__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // cr48 --> he4 + ti44

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  89.75733412085052 + -89.3041938384302 * tfactors.T9i + -81.667 * tfactors.T913i + -10.6333 * tfactors.T913
                         + -0.672613 * tfactors.T9 + 0.161209 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  89.3041938384302 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -81.667 * tfactors.T943i + (1.0/3.0) * -10.6333 * tfactors.T923i
                                  + -0.672613 + (5.0/3.0) * 0.161209 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real cr48_pf, dcr48_pf_dT;
    // interpolating cr48 partition function
    get_partition_function(Cr48, tfactors, cr48_pf, dcr48_pf_dT);

    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real ti44_pf, dti44_pf_dT;
    // interpolating ti44 partition function
    get_partition_function(Ti44, tfactors, ti44_pf, dti44_pf_dT);

    Real z_r = he4_pf * ti44_pf;
    Real z_p = cr48_pf;

    Real dz_r_dT = ti44_pf * dhe4_pf_dT + he4_pf * dti44_pf_dT;
    Real dz_p_dT = dcr48_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_cr48__p_v47__derived__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // cr48 --> p + v47

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // nfisn
    ln_set_rate =  67.7402710725382 + -100.12633621113824 * tfactors.T9i + -34.0548 * tfactors.T913i + -3.41973 * tfactors.T913
                         + 1.16501 * tfactors.T9 + -0.105543 * tfactors.T953 + -6.20886 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  100.12633621113824 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -34.0548 * tfactors.T943i + (1.0/3.0) * -3.41973 * tfactors.T923i
                                  + 1.16501 + (5.0/3.0) * -0.105543 * tfactors.T923 + -6.20886 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // nfisn
    ln_set_rate =  536.5234710725382 + -99.36194621113825 * tfactors.T9i + 317.171 * tfactors.T913i + -911.679 * tfactors.T913
                         + 94.4245 * tfactors.T9 + -10.1973 * tfactors.T953 + 332.227 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  99.36194621113825 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 317.171 * tfactors.T943i + (1.0/3.0) * -911.679 * tfactors.T923i
                                  + 94.4245 + (5.0/3.0) * -10.1973 * tfactors.T923 + 332.227 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // nfisn
    ln_set_rate =  48.89197107253821 + -93.82037121113825 * tfactors.T9i + -45.9868 * tfactors.T913i + 13.6822 * tfactors.T913
                         + -0.376902 * tfactors.T9 + -0.0194875 * tfactors.T953 + -6.9232499999999995 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  93.82037121113825 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -45.9868 * tfactors.T943i + (1.0/3.0) * 13.6822 * tfactors.T923i
                                  + -0.376902 + (5.0/3.0) * -0.0194875 * tfactors.T923 + -6.9232499999999995 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

    // nfisn
    ln_set_rate =  65.62307107253821 + -94.58145021113825 * tfactors.T9i + -110.655 * tfactors.T913i + 83.0232 * tfactors.T913
                         + -19.7762 * tfactors.T9 + 3.03961 * tfactors.T953 + -47.9742 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  94.58145021113825 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -110.655 * tfactors.T943i + (1.0/3.0) * 83.0232 * tfactors.T923i
                                  + -19.7762 + (5.0/3.0) * 3.03961 * tfactors.T923 + -47.9742 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real cr48_pf, dcr48_pf_dT;
    // interpolating cr48 partition function
    get_partition_function(Cr48, tfactors, cr48_pf, dcr48_pf_dT);

    Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    Real v47_pf, dv47_pf_dT;
    // interpolating v47 partition function
    get_partition_function(V47, tfactors, v47_pf, dv47_pf_dT);

    Real z_r = p_pf * v47_pf;
    Real z_p = cr48_pf;

    Real dz_r_dT = v47_pf * dp_pf_dT + p_pf * dv47_pf_dT;
    Real dz_p_dT = dcr48_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_v47__he4_ti44__derived__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // v47 + p --> he4 + ti44

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // chw0r
    ln_set_rate =  -76.61433695168769 + -6.029445311023438 * tfactors.T9i + 70.2835 * tfactors.T913
                         + -7.99061 * tfactors.T9 + 0.486213 * tfactors.T953 + -1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  6.029445311023438 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 70.2835 * tfactors.T923i
                                  + -7.99061 + (5.0/3.0) * 0.486213 * tfactors.T923 + -1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real v47_pf, dv47_pf_dT;
    // interpolating v47 partition function
    get_partition_function(V47, tfactors, v47_pf, dv47_pf_dT);

    Real ti44_pf, dti44_pf_dT;
    // interpolating ti44 partition function
    get_partition_function(Ti44, tfactors, ti44_pf, dti44_pf_dT);

    Real z_r = he4_pf * ti44_pf;
    Real z_p = p_pf * v47_pf;

    Real dz_r_dT = ti44_pf * dhe4_pf_dT + he4_pf * dti44_pf_dT;
    Real dz_p_dT = v47_pf * dp_pf_dT + p_pf * dv47_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_fe52__he4_cr48__derived__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // fe52 --> he4 + cr48

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  90.14738712482465 + -92.10912191363732 * tfactors.T9i + -86.7459 * tfactors.T913i + -9.79373 * tfactors.T913
                         + -0.772169 * tfactors.T9 + 0.155883 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  92.10912191363732 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -86.7459 * tfactors.T943i + (1.0/3.0) * -9.79373 * tfactors.T923i
                                  + -0.772169 + (5.0/3.0) * 0.155883 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real cr48_pf, dcr48_pf_dT;
    // interpolating cr48 partition function
    get_partition_function(Cr48, tfactors, cr48_pf, dcr48_pf_dT);

    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real fe52_pf, dfe52_pf_dT;
    // interpolating fe52 partition function
    get_partition_function(Fe52, tfactors, fe52_pf, dfe52_pf_dT);

    Real z_r = he4_pf * cr48_pf;
    Real z_p = fe52_pf;

    Real dz_r_dT = cr48_pf * dhe4_pf_dT + he4_pf * dcr48_pf_dT;
    Real dz_p_dT = dfe52_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_fe52__p_mn51__derived__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // fe52 --> p + mn51

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  61.72798916565747 + -85.6326403498911 * tfactors.T9i + -36.1825 * tfactors.T913i + 0.873042 * tfactors.T913
                         + -2.89731 * tfactors.T9 + 0.364394 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  85.6326403498911 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -36.1825 * tfactors.T943i + (1.0/3.0) * 0.873042 * tfactors.T923i
                                  + -2.89731 + (5.0/3.0) * 0.364394 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real mn51_pf, dmn51_pf_dT;
    // interpolating mn51 partition function
    get_partition_function(Mn51, tfactors, mn51_pf, dmn51_pf_dT);

    Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    Real fe52_pf, dfe52_pf_dT;
    // interpolating fe52 partition function
    get_partition_function(Fe52, tfactors, fe52_pf, dfe52_pf_dT);

    Real z_r = p_pf * mn51_pf;
    Real z_p = fe52_pf;

    Real dz_r_dT = mn51_pf * dp_pf_dT + p_pf * dmn51_pf_dT;
    Real dz_p_dT = dfe52_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_mn51__he4_cr48__derived__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // mn51 + p --> he4 + cr48

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  58.73119795916719 + -6.476551190854941 * tfactors.T9i + -86.7459 * tfactors.T913i + 1.05653 * tfactors.T913
                         + -1.15757 * tfactors.T9 + 0.0877546 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  6.476551190854941 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -86.7459 * tfactors.T943i + (1.0/3.0) * 1.05653 * tfactors.T923i
                                  + -1.15757 + (5.0/3.0) * 0.0877546 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    Real cr48_pf, dcr48_pf_dT;
    // interpolating cr48 partition function
    get_partition_function(Cr48, tfactors, cr48_pf, dcr48_pf_dT);

    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real mn51_pf, dmn51_pf_dT;
    // interpolating mn51 partition function
    get_partition_function(Mn51, tfactors, mn51_pf, dmn51_pf_dT);

    Real z_r = he4_pf * cr48_pf;
    Real z_p = p_pf * mn51_pf;

    Real dz_r_dT = cr48_pf * dhe4_pf_dT + he4_pf * dcr48_pf_dT;
    Real dz_p_dT = mn51_pf * dp_pf_dT + p_pf * dmn51_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_ni56__he4_fe52__derived__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ni56 --> he4 + fe52

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  91.62258922810437 + -92.801099329237 * tfactors.T9i + -91.6819 * tfactors.T913i + -9.51885 * tfactors.T913
                         + -0.533014 * tfactors.T9 + 0.0892607 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  92.801099329237 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -91.6819 * tfactors.T943i + (1.0/3.0) * -9.51885 * tfactors.T923i
                                  + -0.533014 + (5.0/3.0) * 0.0892607 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real ni56_pf, dni56_pf_dT;
    // interpolating ni56 partition function
    get_partition_function(Ni56, tfactors, ni56_pf, dni56_pf_dT);

    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real fe52_pf, dfe52_pf_dT;
    // interpolating fe52 partition function
    get_partition_function(Fe52, tfactors, fe52_pf, dfe52_pf_dT);

    Real z_r = he4_pf * fe52_pf;
    Real z_p = ni56_pf;

    Real dz_r_dT = fe52_pf * dhe4_pf_dT + he4_pf * dfe52_pf_dT;
    Real dz_p_dT = dni56_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_ni56__p_co55__derived__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // ni56 --> p + co55

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  63.13177060864089 + -83.14741674893808 * tfactors.T9i + -38.1053 * tfactors.T913i + -0.210947 * tfactors.T913
                         + -2.68377 * tfactors.T9 + 0.355814 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  83.14741674893808 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -38.1053 * tfactors.T943i + (1.0/3.0) * -0.210947 * tfactors.T923i
                                  + -2.68377 + (5.0/3.0) * 0.355814 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real ni56_pf, dni56_pf_dT;
    // interpolating ni56 partition function
    get_partition_function(Ni56, tfactors, ni56_pf, dni56_pf_dT);

    Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    Real co55_pf, dco55_pf_dT;
    // interpolating co55 partition function
    get_partition_function(Co55, tfactors, co55_pf, dco55_pf_dT);

    Real z_r = p_pf * co55_pf;
    Real z_p = ni56_pf;

    Real dz_r_dT = co55_pf * dp_pf_dT + p_pf * dco55_pf_dT;
    Real dz_p_dT = dni56_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_co55__he4_fe52__derived__removed(const tf_t& tfactors, Real& rate, Real& drate_dT) {

    // co55 + p --> he4 + fe52

    rate = 0.0;
    drate_dT = 0.0;

    Real ln_set_rate{0.0};
    Real dln_set_rate_dT9{0.0};
    Real set_rate{0.0};

    // ths8r
    ln_set_rate =  61.443418619463486 + -9.65364776674457 * tfactors.T9i + -91.6819 * tfactors.T913i + -0.329235 * tfactors.T913
                         + -0.780924 * tfactors.T9 + 0.0425179 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  9.65364776674457 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -91.6819 * tfactors.T943i + (1.0/3.0) * -0.329235 * tfactors.T923i
                                  + -0.780924 + (5.0/3.0) * 0.0425179 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    Real he4_pf, dhe4_pf_dT;
    // setting he4 partition function to 1.0 by default, independent of T
    he4_pf = 1.0_rt;
    dhe4_pf_dT = 0.0_rt;

    Real co55_pf, dco55_pf_dT;
    // interpolating co55 partition function
    get_partition_function(Co55, tfactors, co55_pf, dco55_pf_dT);

    Real fe52_pf, dfe52_pf_dT;
    // interpolating fe52 partition function
    get_partition_function(Fe52, tfactors, fe52_pf, dfe52_pf_dT);

    Real z_r = he4_pf * fe52_pf;
    Real z_p = p_pf * co55_pf;

    Real dz_r_dT = fe52_pf * dhe4_pf_dT + he4_pf * dfe52_pf_dT;
    Real dz_p_dT = co55_pf * dp_pf_dT + p_pf * dco55_pf_dT;

    Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_s32_he4__ar36__approx(const T& rate_eval, Real& rate, Real& drate_dT) {

    Real r_ag = rate_eval.screened_rates(k_he4_s32__ar36__removed);
    Real r_ap = rate_eval.screened_rates(k_he4_s32__p_cl35__derived__removed);
    Real r_pg = rate_eval.screened_rates(k_p_cl35__ar36__removed);
    Real r_pa = rate_eval.screened_rates(k_p_cl35__he4_s32__removed);
    Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ag + r_ap * r_pg * dd;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        Real drdT_ag = rate_eval.dscreened_rates_dT(k_he4_s32__ar36__removed);
        Real drdT_ap = rate_eval.dscreened_rates_dT(k_he4_s32__p_cl35__derived__removed);
        Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_cl35__ar36__removed);
        Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_cl35__he4_s32__removed);
        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_ar36__s32_he4__approx(const T& rate_eval, Real& rate, Real& drate_dT) {

    Real r_ga = rate_eval.screened_rates(k_ar36__he4_s32__derived__removed);
    Real r_pa = rate_eval.screened_rates(k_p_cl35__he4_s32__removed);
    Real r_gp = rate_eval.screened_rates(k_ar36__p_cl35__derived__removed);
    Real r_pg = rate_eval.screened_rates(k_p_cl35__ar36__removed);
    Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ga + r_gp * r_pa * dd;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        Real drdT_ga = rate_eval.dscreened_rates_dT(k_ar36__he4_s32__derived__removed);
        Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_cl35__he4_s32__removed);
        Real drdT_gp = rate_eval.dscreened_rates_dT(k_ar36__p_cl35__derived__removed);
        Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_cl35__ar36__removed);
        drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_ar36_he4__ca40__approx(const T& rate_eval, Real& rate, Real& drate_dT) {

    Real r_ag = rate_eval.screened_rates(k_he4_ar36__ca40__removed);
    Real r_ap = rate_eval.screened_rates(k_he4_ar36__p_k39__derived__removed);
    Real r_pg = rate_eval.screened_rates(k_p_k39__ca40__removed);
    Real r_pa = rate_eval.screened_rates(k_p_k39__he4_ar36__removed);
    Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ag + r_ap * r_pg * dd;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        Real drdT_ag = rate_eval.dscreened_rates_dT(k_he4_ar36__ca40__removed);
        Real drdT_ap = rate_eval.dscreened_rates_dT(k_he4_ar36__p_k39__derived__removed);
        Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_k39__ca40__removed);
        Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_k39__he4_ar36__removed);
        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_ca40__ar36_he4__approx(const T& rate_eval, Real& rate, Real& drate_dT) {

    Real r_ga = rate_eval.screened_rates(k_ca40__he4_ar36__derived__removed);
    Real r_pa = rate_eval.screened_rates(k_p_k39__he4_ar36__removed);
    Real r_gp = rate_eval.screened_rates(k_ca40__p_k39__derived__removed);
    Real r_pg = rate_eval.screened_rates(k_p_k39__ca40__removed);
    Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ga + r_gp * r_pa * dd;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        Real drdT_ga = rate_eval.dscreened_rates_dT(k_ca40__he4_ar36__derived__removed);
        Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_k39__he4_ar36__removed);
        Real drdT_gp = rate_eval.dscreened_rates_dT(k_ca40__p_k39__derived__removed);
        Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_k39__ca40__removed);
        drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_ca40_he4__ti44__approx(const T& rate_eval, Real& rate, Real& drate_dT) {

    Real r_ag = rate_eval.screened_rates(k_he4_ca40__ti44__removed);
    Real r_ap = rate_eval.screened_rates(k_he4_ca40__p_sc43__derived__removed);
    Real r_pg = rate_eval.screened_rates(k_p_sc43__ti44__removed);
    Real r_pa = rate_eval.screened_rates(k_p_sc43__he4_ca40__removed);
    Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ag + r_ap * r_pg * dd;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        Real drdT_ag = rate_eval.dscreened_rates_dT(k_he4_ca40__ti44__removed);
        Real drdT_ap = rate_eval.dscreened_rates_dT(k_he4_ca40__p_sc43__derived__removed);
        Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_sc43__ti44__removed);
        Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_sc43__he4_ca40__removed);
        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_ti44__ca40_he4__approx(const T& rate_eval, Real& rate, Real& drate_dT) {

    Real r_ga = rate_eval.screened_rates(k_ti44__he4_ca40__derived__removed);
    Real r_pa = rate_eval.screened_rates(k_p_sc43__he4_ca40__removed);
    Real r_gp = rate_eval.screened_rates(k_ti44__p_sc43__derived__removed);
    Real r_pg = rate_eval.screened_rates(k_p_sc43__ti44__removed);
    Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ga + r_gp * r_pa * dd;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        Real drdT_ga = rate_eval.dscreened_rates_dT(k_ti44__he4_ca40__derived__removed);
        Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_sc43__he4_ca40__removed);
        Real drdT_gp = rate_eval.dscreened_rates_dT(k_ti44__p_sc43__derived__removed);
        Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_sc43__ti44__removed);
        drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_ti44_he4__cr48__approx(const T& rate_eval, Real& rate, Real& drate_dT) {

    Real r_ag = rate_eval.screened_rates(k_he4_ti44__cr48__removed);
    Real r_ap = rate_eval.screened_rates(k_he4_ti44__p_v47__removed);
    Real r_pg = rate_eval.screened_rates(k_p_v47__cr48__removed);
    Real r_pa = rate_eval.screened_rates(k_p_v47__he4_ti44__derived__removed);
    Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ag + r_ap * r_pg * dd;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        Real drdT_ag = rate_eval.dscreened_rates_dT(k_he4_ti44__cr48__removed);
        Real drdT_ap = rate_eval.dscreened_rates_dT(k_he4_ti44__p_v47__removed);
        Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_v47__cr48__removed);
        Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_v47__he4_ti44__derived__removed);
        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_cr48__ti44_he4__approx(const T& rate_eval, Real& rate, Real& drate_dT) {

    Real r_ga = rate_eval.screened_rates(k_cr48__he4_ti44__derived__removed);
    Real r_pa = rate_eval.screened_rates(k_p_v47__he4_ti44__derived__removed);
    Real r_gp = rate_eval.screened_rates(k_cr48__p_v47__derived__removed);
    Real r_pg = rate_eval.screened_rates(k_p_v47__cr48__removed);
    Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ga + r_gp * r_pa * dd;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        Real drdT_ga = rate_eval.dscreened_rates_dT(k_cr48__he4_ti44__derived__removed);
        Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_v47__he4_ti44__derived__removed);
        Real drdT_gp = rate_eval.dscreened_rates_dT(k_cr48__p_v47__derived__removed);
        Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_v47__cr48__removed);
        drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_cr48_he4__fe52__approx(const T& rate_eval, Real& rate, Real& drate_dT) {

    Real r_ag = rate_eval.screened_rates(k_he4_cr48__fe52__removed);
    Real r_ap = rate_eval.screened_rates(k_he4_cr48__p_mn51__removed);
    Real r_pg = rate_eval.screened_rates(k_p_mn51__fe52__removed);
    Real r_pa = rate_eval.screened_rates(k_p_mn51__he4_cr48__derived__removed);
    Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ag + r_ap * r_pg * dd;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        Real drdT_ag = rate_eval.dscreened_rates_dT(k_he4_cr48__fe52__removed);
        Real drdT_ap = rate_eval.dscreened_rates_dT(k_he4_cr48__p_mn51__removed);
        Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_mn51__fe52__removed);
        Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_mn51__he4_cr48__derived__removed);
        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_fe52__cr48_he4__approx(const T& rate_eval, Real& rate, Real& drate_dT) {

    Real r_ga = rate_eval.screened_rates(k_fe52__he4_cr48__derived__removed);
    Real r_pa = rate_eval.screened_rates(k_p_mn51__he4_cr48__derived__removed);
    Real r_gp = rate_eval.screened_rates(k_fe52__p_mn51__derived__removed);
    Real r_pg = rate_eval.screened_rates(k_p_mn51__fe52__removed);
    Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ga + r_gp * r_pa * dd;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        Real drdT_ga = rate_eval.dscreened_rates_dT(k_fe52__he4_cr48__derived__removed);
        Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_mn51__he4_cr48__derived__removed);
        Real drdT_gp = rate_eval.dscreened_rates_dT(k_fe52__p_mn51__derived__removed);
        Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_mn51__fe52__removed);
        drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_fe52_he4__ni56__approx(const T& rate_eval, Real& rate, Real& drate_dT) {

    Real r_ag = rate_eval.screened_rates(k_he4_fe52__ni56__removed);
    Real r_ap = rate_eval.screened_rates(k_he4_fe52__p_co55__removed);
    Real r_pg = rate_eval.screened_rates(k_p_co55__ni56__removed);
    Real r_pa = rate_eval.screened_rates(k_p_co55__he4_fe52__derived__removed);
    Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ag + r_ap * r_pg * dd;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        Real drdT_ag = rate_eval.dscreened_rates_dT(k_he4_fe52__ni56__removed);
        Real drdT_ap = rate_eval.dscreened_rates_dT(k_he4_fe52__p_co55__removed);
        Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_co55__ni56__removed);
        Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_co55__he4_fe52__derived__removed);
        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_ni56__fe52_he4__approx(const T& rate_eval, Real& rate, Real& drate_dT) {

    Real r_ga = rate_eval.screened_rates(k_ni56__he4_fe52__derived__removed);
    Real r_pa = rate_eval.screened_rates(k_p_co55__he4_fe52__derived__removed);
    Real r_gp = rate_eval.screened_rates(k_ni56__p_co55__derived__removed);
    Real r_pg = rate_eval.screened_rates(k_p_co55__ni56__removed);
    Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ga + r_gp * r_pa * dd;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        Real drdT_ga = rate_eval.dscreened_rates_dT(k_ni56__he4_fe52__derived__removed);
        Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_co55__he4_fe52__derived__removed);
        Real drdT_gp = rate_eval.dscreened_rates_dT(k_ni56__p_co55__derived__removed);
        Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_co55__ni56__removed);
        drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa);
    }
}


template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
fill_reaclib_rates(const tf_t& tfactors, T& rate_eval)
{

    Real rate;
    Real drate_dT;

    rate_n__p__weak__wc12<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n__p__weak__wc12) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_n__p__weak__wc12) = drate_dT;

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
    rate_n_n13__n14<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_n13__n14) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_n_n13__n14) = drate_dT;

    }
    rate_he4_n14__f18<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_n14__f18) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_n14__f18) = drate_dT;

    }
    rate_he4_o16__ne20<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_o16__ne20) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_o16__ne20) = drate_dT;

    }
    rate_he4_f18__na22<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_f18__na22) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_f18__na22) = drate_dT;

    }
    rate_n_ne20__ne21<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_ne20__ne21) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_n_ne20__ne21) = drate_dT;

    }
    rate_he4_ne20__mg24<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_ne20__mg24) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_ne20__mg24) = drate_dT;

    }
    rate_p_ne21__na22<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_ne21__na22) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_ne21__na22) = drate_dT;

    }
    rate_n_na22__na23<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_na22__na23) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_n_na22__na23) = drate_dT;

    }
    rate_p_na23__mg24<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_na23__mg24) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_na23__mg24) = drate_dT;

    }
    rate_he4_mg24__si28<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_mg24__si28) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_mg24__si28) = drate_dT;

    }
    rate_p_al27__si28<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_al27__si28) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_al27__si28) = drate_dT;

    }
    rate_he4_si28__s32<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_si28__s32) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_si28__s32) = drate_dT;

    }
    rate_p_p31__s32<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_p31__s32) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_p31__s32) = drate_dT;

    }
    rate_c12_c12__p_na23<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_c12_c12__p_na23) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_c12_c12__p_na23) = drate_dT;

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
    rate_c12_o16__p_al27<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_c12_o16__p_al27) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_c12_o16__p_al27) = drate_dT;

    }
    rate_c12_o16__he4_mg24<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_c12_o16__he4_mg24) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_c12_o16__he4_mg24) = drate_dT;

    }
    rate_o16_o16__p_p31<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_o16_o16__p_p31) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_o16_o16__p_p31) = drate_dT;

    }
    rate_o16_o16__he4_si28<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_o16_o16__he4_si28) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_o16_o16__he4_si28) = drate_dT;

    }
    rate_he4_f18__p_ne21<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_f18__p_ne21) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_f18__p_ne21) = drate_dT;

    }
    rate_he4_ne21__n_mg24<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_ne21__n_mg24) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_ne21__n_mg24) = drate_dT;

    }
    rate_p_na23__he4_ne20<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_na23__he4_ne20) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_na23__he4_ne20) = drate_dT;

    }
    rate_p_al27__he4_mg24<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_al27__he4_mg24) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_al27__he4_mg24) = drate_dT;

    }
    rate_p_p31__he4_si28<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_p31__he4_si28) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_p31__he4_si28) = drate_dT;

    }
    rate_he4_he4_he4__c12<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_he4_he4__c12) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_he4_he4__c12) = drate_dT;

    }
    rate_c12_c12__mg24__modified<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_c12_c12__mg24__modified) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_c12_c12__mg24__modified) = drate_dT;

    }
    rate_o16_o16__s32__modified<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_o16_o16__s32__modified) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_o16_o16__s32__modified) = drate_dT;

    }
    rate_c12_o16__si28__modified<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_c12_o16__si28__modified) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_c12_o16__si28__modified) = drate_dT;

    }
    rate_he4_s32__ar36__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_s32__ar36__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_s32__ar36__removed) = drate_dT;

    }
    rate_p_cl35__ar36__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_cl35__ar36__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_cl35__ar36__removed) = drate_dT;

    }
    rate_p_cl35__he4_s32__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_cl35__he4_s32__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_cl35__he4_s32__removed) = drate_dT;

    }
    rate_he4_ar36__ca40__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_ar36__ca40__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_ar36__ca40__removed) = drate_dT;

    }
    rate_p_k39__ca40__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_k39__ca40__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_k39__ca40__removed) = drate_dT;

    }
    rate_p_k39__he4_ar36__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_k39__he4_ar36__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_k39__he4_ar36__removed) = drate_dT;

    }
    rate_he4_ca40__ti44__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_ca40__ti44__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_ca40__ti44__removed) = drate_dT;

    }
    rate_p_sc43__ti44__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_sc43__ti44__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_sc43__ti44__removed) = drate_dT;

    }
    rate_p_sc43__he4_ca40__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_sc43__he4_ca40__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_sc43__he4_ca40__removed) = drate_dT;

    }
    rate_he4_ti44__cr48__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_ti44__cr48__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_ti44__cr48__removed) = drate_dT;

    }
    rate_he4_ti44__p_v47__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_ti44__p_v47__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_ti44__p_v47__removed) = drate_dT;

    }
    rate_p_v47__cr48__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_v47__cr48__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_v47__cr48__removed) = drate_dT;

    }
    rate_he4_cr48__fe52__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_cr48__fe52__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_cr48__fe52__removed) = drate_dT;

    }
    rate_he4_cr48__p_mn51__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_cr48__p_mn51__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_cr48__p_mn51__removed) = drate_dT;

    }
    rate_p_mn51__fe52__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_mn51__fe52__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_mn51__fe52__removed) = drate_dT;

    }
    rate_he4_fe52__ni56__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_fe52__ni56__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_fe52__ni56__removed) = drate_dT;

    }
    rate_he4_fe52__p_co55__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_fe52__p_co55__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_fe52__p_co55__removed) = drate_dT;

    }
    rate_p_co55__ni56__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_co55__ni56__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_co55__ni56__removed) = drate_dT;

    }
    rate_n13__p_c12__derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n13__p_c12__derived) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_n13__p_c12__derived) = drate_dT;

    }
    rate_n14__n_n13__derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n14__n_n13__derived) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_n14__n_n13__derived) = drate_dT;

    }
    rate_o16__he4_c12__derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_o16__he4_c12__derived) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_o16__he4_c12__derived) = drate_dT;

    }
    rate_f18__he4_n14__derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_f18__he4_n14__derived) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_f18__he4_n14__derived) = drate_dT;

    }
    rate_ne20__he4_o16__derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_ne20__he4_o16__derived) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_ne20__he4_o16__derived) = drate_dT;

    }
    rate_ne21__n_ne20__derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_ne21__n_ne20__derived) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_ne21__n_ne20__derived) = drate_dT;

    }
    rate_na22__p_ne21__derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_na22__p_ne21__derived) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_na22__p_ne21__derived) = drate_dT;

    }
    rate_na22__he4_f18__derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_na22__he4_f18__derived) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_na22__he4_f18__derived) = drate_dT;

    }
    rate_na23__n_na22__derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_na23__n_na22__derived) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_na23__n_na22__derived) = drate_dT;

    }
    rate_mg24__p_na23__derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_mg24__p_na23__derived) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_mg24__p_na23__derived) = drate_dT;

    }
    rate_mg24__he4_ne20__derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_mg24__he4_ne20__derived) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_mg24__he4_ne20__derived) = drate_dT;

    }
    rate_si28__p_al27__derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_si28__p_al27__derived) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_si28__p_al27__derived) = drate_dT;

    }
    rate_si28__he4_mg24__derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_si28__he4_mg24__derived) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_si28__he4_mg24__derived) = drate_dT;

    }
    rate_s32__p_p31__derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_s32__p_p31__derived) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_s32__p_p31__derived) = drate_dT;

    }
    rate_s32__he4_si28__derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_s32__he4_si28__derived) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_s32__he4_si28__derived) = drate_dT;

    }
    rate_c12__he4_he4_he4__derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_c12__he4_he4_he4__derived) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_c12__he4_he4_he4__derived) = drate_dT;

    }
    rate_p_o16__he4_n13__derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_o16__he4_n13__derived) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_o16__he4_n13__derived) = drate_dT;

    }
    rate_he4_ne20__p_na23__derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_ne20__p_na23__derived) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_ne20__p_na23__derived) = drate_dT;

    }
    rate_p_ne21__he4_f18__derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_ne21__he4_f18__derived) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_ne21__he4_f18__derived) = drate_dT;

    }
    rate_n_mg24__he4_ne21__derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_mg24__he4_ne21__derived) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_n_mg24__he4_ne21__derived) = drate_dT;

    }
    rate_he4_mg24__p_al27__derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_mg24__p_al27__derived) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_mg24__p_al27__derived) = drate_dT;

    }
    rate_he4_si28__p_p31__derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_si28__p_p31__derived) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_si28__p_p31__derived) = drate_dT;

    }
    rate_he4_s32__p_cl35__derived__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_s32__p_cl35__derived__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_s32__p_cl35__derived__removed) = drate_dT;

    }
    rate_ar36__he4_s32__derived__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_ar36__he4_s32__derived__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_ar36__he4_s32__derived__removed) = drate_dT;

    }
    rate_ar36__p_cl35__derived__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_ar36__p_cl35__derived__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_ar36__p_cl35__derived__removed) = drate_dT;

    }
    rate_he4_ar36__p_k39__derived__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_ar36__p_k39__derived__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_ar36__p_k39__derived__removed) = drate_dT;

    }
    rate_ca40__he4_ar36__derived__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_ca40__he4_ar36__derived__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_ca40__he4_ar36__derived__removed) = drate_dT;

    }
    rate_ca40__p_k39__derived__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_ca40__p_k39__derived__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_ca40__p_k39__derived__removed) = drate_dT;

    }
    rate_he4_ca40__p_sc43__derived__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_he4_ca40__p_sc43__derived__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_he4_ca40__p_sc43__derived__removed) = drate_dT;

    }
    rate_ti44__he4_ca40__derived__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_ti44__he4_ca40__derived__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_ti44__he4_ca40__derived__removed) = drate_dT;

    }
    rate_ti44__p_sc43__derived__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_ti44__p_sc43__derived__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_ti44__p_sc43__derived__removed) = drate_dT;

    }
    rate_cr48__he4_ti44__derived__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_cr48__he4_ti44__derived__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_cr48__he4_ti44__derived__removed) = drate_dT;

    }
    rate_cr48__p_v47__derived__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_cr48__p_v47__derived__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_cr48__p_v47__derived__removed) = drate_dT;

    }
    rate_p_v47__he4_ti44__derived__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_v47__he4_ti44__derived__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_v47__he4_ti44__derived__removed) = drate_dT;

    }
    rate_fe52__he4_cr48__derived__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_fe52__he4_cr48__derived__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_fe52__he4_cr48__derived__removed) = drate_dT;

    }
    rate_fe52__p_mn51__derived__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_fe52__p_mn51__derived__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_fe52__p_mn51__derived__removed) = drate_dT;

    }
    rate_p_mn51__he4_cr48__derived__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_mn51__he4_cr48__derived__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_mn51__he4_cr48__derived__removed) = drate_dT;

    }
    rate_ni56__he4_fe52__derived__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_ni56__he4_fe52__derived__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_ni56__he4_fe52__derived__removed) = drate_dT;

    }
    rate_ni56__p_co55__derived__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_ni56__p_co55__derived__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_ni56__p_co55__derived__removed) = drate_dT;

    }
    rate_p_co55__he4_fe52__derived__removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_co55__he4_fe52__derived__removed) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_p_co55__he4_fe52__derived__removed) = drate_dT;

    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
fill_approx_rates([[maybe_unused]] const tf_t& tfactors, T& rate_eval)
{

    Real rate;
    Real drate_dT;

    rate_s32_he4__ar36__approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_s32_he4__ar36__approx) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_s32_he4__ar36__approx) = drate_dT;

    }
    rate_ar36__s32_he4__approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_ar36__s32_he4__approx) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_ar36__s32_he4__approx) = drate_dT;

    }
    rate_ar36_he4__ca40__approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_ar36_he4__ca40__approx) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_ar36_he4__ca40__approx) = drate_dT;

    }
    rate_ca40__ar36_he4__approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_ca40__ar36_he4__approx) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_ca40__ar36_he4__approx) = drate_dT;

    }
    rate_ca40_he4__ti44__approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_ca40_he4__ti44__approx) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_ca40_he4__ti44__approx) = drate_dT;

    }
    rate_ti44__ca40_he4__approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_ti44__ca40_he4__approx) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_ti44__ca40_he4__approx) = drate_dT;

    }
    rate_ti44_he4__cr48__approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_ti44_he4__cr48__approx) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_ti44_he4__cr48__approx) = drate_dT;

    }
    rate_cr48__ti44_he4__approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_cr48__ti44_he4__approx) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_cr48__ti44_he4__approx) = drate_dT;

    }
    rate_cr48_he4__fe52__approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_cr48_he4__fe52__approx) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_cr48_he4__fe52__approx) = drate_dT;

    }
    rate_fe52__cr48_he4__approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_fe52__cr48_he4__approx) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_fe52__cr48_he4__approx) = drate_dT;

    }
    rate_fe52_he4__ni56__approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_fe52_he4__ni56__approx) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_fe52_he4__ni56__approx) = drate_dT;

    }
    rate_ni56__fe52_he4__approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_ni56__fe52_he4__approx) = rate;
    if constexpr (std::is_same<T, rate_derivs_t>::value) {
        rate_eval.dscreened_rates_dT(k_ni56__fe52_he4__approx) = drate_dT;

    }

}

#endif
