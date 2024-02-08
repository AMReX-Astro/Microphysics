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
void rate_p_C12_to_N13(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
void rate_He4_C12_to_O16(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
void rate_He4_N14_to_F18(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // N14 + He4 --> F18

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_He4_O16_to_Ne20(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // O16 + He4 --> Ne20

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_He4_F18_to_Na22(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // F18 + He4 --> Na22

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_He4_Ne20_to_Mg24(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ne20 + He4 --> Mg24

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_p_Ne21_to_Na22(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ne21 + p --> Na22

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_p_Na23_to_Mg24(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Na23 + p --> Mg24

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_He4_Mg24_to_Si28(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Mg24 + He4 --> Si28

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_p_Al27_to_Si28(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Al27 + p --> Si28

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_He4_Si28_to_S32(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Si28 + He4 --> S32

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_p_P31_to_S32(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // P31 + p --> S32

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_He4_Cr48_to_Fe52(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Cr48 + He4 --> Fe52

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_p_Mn51_to_Fe52(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Mn51 + p --> Fe52

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_He4_Mn51_to_Co55(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Mn51 + He4 --> Co55

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  65.9219 + -89.274 * tfactors.T913i + -10.4373 * tfactors.T913
                         + 1.00492 * tfactors.T9 + -0.125548 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -89.274 * tfactors.T943i + (1.0/3.0) * -10.4373 * tfactors.T923i
                                  + 1.00492 + (5.0/3.0) * -0.125548 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_He4_Fe52_to_Ni56(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe52 + He4 --> Ni56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_p_Co55_to_Ni56(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co55 + p --> Ni56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_He4_Co55_to_Cu59(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co55 + He4 --> Cu59

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  65.2921 + -94.1404 * tfactors.T913i + -2.62786 * tfactors.T913
                         + -2.12066 * tfactors.T9 + 0.237999 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -94.1404 * tfactors.T943i + (1.0/3.0) * -2.62786 * tfactors.T923i
                                  + -2.12066 + (5.0/3.0) * 0.237999 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_He4_Ni56_to_Zn60(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni56 + He4 --> Zn60

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  61.0733 + -96.4898 * tfactors.T913i + 6.47209 * tfactors.T913
                         + -5.2029 * tfactors.T9 + 0.533391 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -96.4898 * tfactors.T943i + (1.0/3.0) * 6.47209 * tfactors.T923i
                                  + -5.2029 + (5.0/3.0) * 0.533391 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_p_Cu59_to_Zn60(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Cu59 + p --> Zn60

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  37.226 + -39.9808 * tfactors.T913i + 1.17311 * tfactors.T913
                         + -2.90486 * tfactors.T9 + 0.339644 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -39.9808 * tfactors.T943i + (1.0/3.0) * 1.17311 * tfactors.T923i
                                  + -2.90486 + (5.0/3.0) * 0.339644 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_C12_C12_to_p_Na23(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // C12 + C12 --> p + Na23

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_C12_C12_to_He4_Ne20(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_N13_to_p_O16(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_C12_O16_to_p_Al27(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // O16 + C12 --> p + Al27

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_C12_O16_to_He4_Mg24(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_O16_O16_to_p_P31(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // O16 + O16 --> p + P31

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_O16_O16_to_He4_Si28(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // O16 + O16 --> He4 + Si28

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_He4_F18_to_p_Ne21(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // F18 + He4 --> p + Ne21

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_p_Na23_to_He4_Ne20(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Na23 + p --> He4 + Ne20

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_p_Al27_to_He4_Mg24(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Al27 + p --> He4 + Mg24

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_p_P31_to_He4_Si28(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // P31 + p --> He4 + Si28

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_He4_Cr48_to_p_Mn51(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Cr48 + He4 --> p + Mn51

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_He4_Fe52_to_p_Co55(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe52 + He4 --> p + Co55

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_p_Cu59_to_He4_Ni56(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Cu59 + p --> He4 + Ni56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  20.7688 + -39.9808 * tfactors.T913i + 13.6773 * tfactors.T913
                         + -3.76429 * tfactors.T9 + 0.438096 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -39.9808 * tfactors.T943i + (1.0/3.0) * 13.6773 * tfactors.T923i
                                  + -3.76429 + (5.0/3.0) * 0.438096 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_He4_He4_He4_to_C12(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // He4 + He4 + He4 --> C12

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_C12_C12_to_Mg24_modified(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // C12 + C12 --> Mg24

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_O16_O16_to_S32_modified(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // O16 + O16 --> S32

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_C12_O16_to_Si28_modified(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // O16 + C12 --> Si28

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_p_Mn55_to_Fe56(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Mn55 + p --> Fe56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  36.1888 + -36.1995 * tfactors.T913i + 1.15346 * tfactors.T913
                         + -1.94437 * tfactors.T9 + 0.218429 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -36.1995 * tfactors.T943i + (1.0/3.0) * 1.15346 * tfactors.T923i
                                  + -1.94437 + (5.0/3.0) * 0.218429 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_n_Fe52_to_Fe53(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe52 + n --> Fe53

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  15.8885 + -0.344319 * tfactors.T913
                         + 0.178277 * tfactors.T9 + -0.0334326 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * -0.344319 * tfactors.T923i
                                  + 0.178277 + (5.0/3.0) * -0.0334326 * tfactors.T923;
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
void rate_n_Fe53_to_Fe54(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe53 + n --> Fe54

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  16.4534 + -1.10421 * tfactors.T913
                         + 0.379905 * tfactors.T9 + -0.0581878 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * -1.10421 * tfactors.T923i
                                  + 0.379905 + (5.0/3.0) * -0.0581878 * tfactors.T923;
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
void rate_He4_Fe53_to_Ni57(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe53 + He4 --> Ni57

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  69.0636 + -91.7231 * tfactors.T913i + -10.1816 * tfactors.T913
                         + -0.0406303 * tfactors.T9 + 0.0345056 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -91.7231 * tfactors.T943i + (1.0/3.0) * -10.1816 * tfactors.T923i
                                  + -0.0406303 + (5.0/3.0) * 0.0345056 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_n_Fe54_to_Fe55(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe54 + n --> Fe55

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ks03 
    ln_set_rate =  -0.80864 + 0.0591716 * tfactors.T9i + -8.66617 * tfactors.T913i + 26.4472 * tfactors.T913
                         + -1.9222 * tfactors.T9 + 0.0986404 * tfactors.T953 + -9.78317 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  -0.0591716 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -8.66617 * tfactors.T943i + (1.0/3.0) * 26.4472 * tfactors.T923i
                                  + -1.9222 + (5.0/3.0) * 0.0986404 * tfactors.T923 + -9.78317 * tfactors.T9i;
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
void rate_p_Fe54_to_Co55(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe54 + p --> Co55

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  36.2304 + -37.1544 * tfactors.T913i + 0.950364 * tfactors.T913
                         + -1.77529 * tfactors.T9 + 0.198562 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -37.1544 * tfactors.T943i + (1.0/3.0) * 0.950364 * tfactors.T923i
                                  + -1.77529 + (5.0/3.0) * 0.198562 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_He4_Fe54_to_Ni58(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe54 + He4 --> Ni58

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  60.2478 + -91.7628 * tfactors.T913i + 4.23027 * tfactors.T913
                         + -3.31305 * tfactors.T9 + 0.271293 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -91.7628 * tfactors.T943i + (1.0/3.0) * 4.23027 * tfactors.T923i
                                  + -3.31305 + (5.0/3.0) * 0.271293 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_n_Fe55_to_Fe56(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe55 + n --> Fe56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ks03 
    ln_set_rate =  21.7202 + -0.0955677 * tfactors.T9i + 8.06062 * tfactors.T913i + -14.4809 * tfactors.T913
                         + 0.94252 * tfactors.T9 + -0.0776007 * tfactors.T953 + 6.47093 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.0955677 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 8.06062 * tfactors.T943i + (1.0/3.0) * -14.4809 * tfactors.T923i
                                  + 0.94252 + (5.0/3.0) * -0.0776007 * tfactors.T923 + 6.47093 * tfactors.T9i;
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
void rate_p_Fe55_to_Co56(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe55 + p --> Co56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  35.7333 + -37.1585 * tfactors.T913i + 1.66198 * tfactors.T913
                         + -1.60842 * tfactors.T9 + 0.148916 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -37.1585 * tfactors.T943i + (1.0/3.0) * 1.66198 * tfactors.T923i
                                  + -1.60842 + (5.0/3.0) * 0.148916 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_p_Fe56_to_Co57(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe56 + p --> Co57

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  36.0665 + -37.1625 * tfactors.T913i + 1.06776 * tfactors.T913
                         + -1.31689 * tfactors.T9 + 0.122089 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -37.1625 * tfactors.T943i + (1.0/3.0) * 1.06776 * tfactors.T923i
                                  + -1.31689 + (5.0/3.0) * 0.122089 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_n_Co55_to_Co56(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co55 + n --> Co56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  17.065 + -1.86357 * tfactors.T913
                         + 0.616591 * tfactors.T9 + -0.0839313 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * -1.86357 * tfactors.T923i
                                  + 0.616591 + (5.0/3.0) * -0.0839313 * tfactors.T923;
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
void rate_n_Co56_to_Co57(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co56 + n --> Co57

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  17.3552 + -1.37855 * tfactors.T913
                         + 0.299896 * tfactors.T9 + -0.04382 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * -1.37855 * tfactors.T923i
                                  + 0.299896 + (5.0/3.0) * -0.04382 * tfactors.T923;
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
void rate_p_Co56_to_Ni57(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co56 + p --> Ni57

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  35.5937 + -38.1094 * tfactors.T913i + 2.57091 * tfactors.T913
                         + -2.07795 * tfactors.T9 + 0.20757 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -38.1094 * tfactors.T943i + (1.0/3.0) * 2.57091 * tfactors.T923i
                                  + -2.07795 + (5.0/3.0) * 0.20757 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_p_Co57_to_Ni58(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co57 + p --> Ni58

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  36.0159 + -38.1133 * tfactors.T913i + 1.77414 * tfactors.T913
                         + -1.48268 * tfactors.T9 + 0.121073 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -38.1133 * tfactors.T943i + (1.0/3.0) * 1.77414 * tfactors.T923i
                                  + -1.48268 + (5.0/3.0) * 0.121073 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_n_Ni56_to_Ni57(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni56 + n --> Ni57

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  16.0765 + -1.19665 * tfactors.T913
                         + 0.507179 * tfactors.T9 + -0.074604 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * -1.19665 * tfactors.T923i
                                  + 0.507179 + (5.0/3.0) * -0.074604 * tfactors.T923;
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
void rate_n_Ni57_to_Ni58(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni57 + n --> Ni58

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  17.2731 + -1.90814 * tfactors.T913
                         + 0.493188 * tfactors.T9 + -0.0684633 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * -1.90814 * tfactors.T923i
                                  + 0.493188 + (5.0/3.0) * -0.0684633 * tfactors.T923;
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
void rate_p_Ni58_to_Cu59(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni58 + p --> Cu59

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // si13n
    ln_set_rate =  36.6662 + -39.0526 * tfactors.T913i + 1.0436 * tfactors.T913
                         + -2.10834 * tfactors.T9 + 0.239559 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -39.0526 * tfactors.T943i + (1.0/3.0) * 1.0436 * tfactors.T923i
                                  + -2.10834 + (5.0/3.0) * 0.239559 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_He4_Mn51_to_p_Fe54(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Mn51 + He4 --> p + Fe54

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  62.2777 + -89.274 * tfactors.T913i + -0.862452 * tfactors.T913
                         + -0.635672 * tfactors.T9 + 0.0196464 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -89.274 * tfactors.T943i + (1.0/3.0) * -0.862452 * tfactors.T923i
                                  + -0.635672 + (5.0/3.0) * 0.0196464 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_He4_Fe53_to_p_Co56(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe53 + He4 --> p + Co56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  54.8614 + -91.7231 * tfactors.T913i + 10.052 * tfactors.T913
                         + -3.86332 * tfactors.T9 + 0.391105 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -91.7231 * tfactors.T943i + (1.0/3.0) * 10.052 * tfactors.T923i
                                  + -3.86332 + (5.0/3.0) * 0.391105 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_n_Fe55_to_p_Mn55(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe55 + n --> p + Mn55

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  4.70639 + 7.03052 * tfactors.T913
                         + 0.260509 * tfactors.T9 + -0.0813397 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * 7.03052 * tfactors.T923i
                                  + 0.260509 + (5.0/3.0) * -0.0813397 * tfactors.T923;
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
void rate_n_Co55_to_p_Fe55(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co55 + n --> p + Fe55

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  20.7398 + -1.62382 * tfactors.T913
                         + 0.58115 * tfactors.T9 + -0.0537057 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * -1.62382 * tfactors.T923i
                                  + 0.58115 + (5.0/3.0) * -0.0537057 * tfactors.T923;
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
void rate_He4_Co55_to_p_Ni58(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co55 + He4 --> p + Ni58

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  60.2281 + -94.1404 * tfactors.T913i + 3.39179 * tfactors.T913
                         + -1.71062 * tfactors.T9 + 0.133003 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -94.1404 * tfactors.T943i + (1.0/3.0) * 3.39179 * tfactors.T923i
                                  + -1.71062 + (5.0/3.0) * 0.133003 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_n_Co56_to_p_Fe56(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co56 + n --> p + Fe56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  20.4539 + -1.13331 * tfactors.T913
                         + 0.347185 * tfactors.T9 + -0.0328879 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * -1.13331 * tfactors.T923i
                                  + 0.347185 + (5.0/3.0) * -0.0328879 * tfactors.T923;
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
void rate_p_Co57_to_He4_Fe54(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co57 + p --> He4 + Fe54

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  -2.1962 + -38.1133 * tfactors.T913i + 29.3541 * tfactors.T913
                         + -4.75966 * tfactors.T9 + 0.40418 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -38.1133 * tfactors.T943i + (1.0/3.0) * 29.3541 * tfactors.T923i
                                  + -4.75966 + (5.0/3.0) * 0.40418 * tfactors.T923 + -0.666667 * tfactors.T9i;
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
void rate_n_Ni56_to_p_Co56(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni56 + n --> p + Co56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  15.5693 + 1.76846 * tfactors.T913
                         + 0.197992 * tfactors.T9 + -0.017494 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * 1.76846 * tfactors.T923i
                                  + 0.197992 + (5.0/3.0) * -0.017494 * tfactors.T923;
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
void rate_n_Ni56_to_He4_Fe53(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni56 + n --> He4 + Fe53

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  -7.84393 + 6.43259 * tfactors.T913
                         + 2.29591 * tfactors.T9 + -0.254724 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * 6.43259 * tfactors.T923i
                                  + 2.29591 + (5.0/3.0) * -0.254724 * tfactors.T923;
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
void rate_n_Ni57_to_p_Co57(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni57 + n --> p + Co57

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  19.9939 + -1.13177 * tfactors.T913
                         + 0.56875 * tfactors.T9 + -0.0579913 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * -1.13177 * tfactors.T923i
                                  + 0.56875 + (5.0/3.0) * -0.0579913 * tfactors.T923;
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
void rate_n_Ni57_to_He4_Fe54(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni57 + n --> He4 + Fe54

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  14.5327 + 0.48562 * tfactors.T913
                         + -0.0249302 * tfactors.T9 + 0.0509605 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * 0.48562 * tfactors.T923i
                                  + -0.0249302 + (5.0/3.0) * 0.0509605 * tfactors.T923;
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
void rate_n_Ni58_to_He4_Fe55(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni58 + n --> He4 + Fe55

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  -2.03901 + 3.44996 * tfactors.T913
                         + 2.98226 * tfactors.T9 + -0.387699 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * 3.44996 * tfactors.T923i
                                  + 2.98226 + (5.0/3.0) * -0.387699 * tfactors.T923;
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
void rate_n_Cu59_to_He4_Co56(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Cu59 + n --> He4 + Co56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  11.0108 + -1.00882 * tfactors.T913
                         + 1.88993 * tfactors.T9 + -0.202319 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * -1.00882 * tfactors.T923i
                                  + 1.88993 + (5.0/3.0) * -0.202319 * tfactors.T923;
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
void rate_n_Zn60_to_He4_Ni57(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Zn60 + n --> He4 + Ni57

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  12.4241 + 5.43552 * tfactors.T913
                         + -1.25169 * tfactors.T9 + 0.148383 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * 5.43552 * tfactors.T923i
                                  + -1.25169 + (5.0/3.0) * 0.148383 * tfactors.T923;
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
void rate_He4_S32_to_Ar36_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // S32 + He4 --> Ar36

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_p_Cl35_to_Ar36_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Cl35 + p --> Ar36

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_p_Cl35_to_He4_S32_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Cl35 + p --> He4 + S32

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_He4_Ar36_to_Ca40_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ar36 + He4 --> Ca40

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_p_K39_to_Ca40_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // K39 + p --> Ca40

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_p_K39_to_He4_Ar36_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // K39 + p --> He4 + Ar36

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_He4_Ca40_to_Ti44_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ca40 + He4 --> Ti44

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_p_Sc43_to_Ti44_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Sc43 + p --> Ti44

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_p_Sc43_to_He4_Ca40_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Sc43 + p --> He4 + Ca40

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_He4_Ti44_to_Cr48_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ti44 + He4 --> Cr48

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_He4_Ti44_to_p_V47_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ti44 + He4 --> p + V47

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_p_V47_to_Cr48_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // V47 + p --> Cr48

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
void rate_N13_to_p_C12_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // N13 --> p + C12

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ls09c
    ln_set_rate =  40.04074558314484 + -22.547578710551672 * tfactors.T9i + -13.692 * tfactors.T913i + -0.230881 * tfactors.T913
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
    ln_set_rate =  40.43534558314484 + -26.326068710551674 * tfactors.T9i + -5.10735 * tfactors.T913i + -2.24111 * tfactors.T913
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

    amrex::Real dz_r_dT = C12_pf * dp_pf_dT + p_pf * dC12_pf_dT;
    amrex::Real dz_p_dT = dN13_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_O16_to_He4_C12_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // O16 --> He4 + C12

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
    ln_set_rate =  94.31312807765731 + -84.50317042649215 * tfactors.T9i + 58.9128 * tfactors.T913i + -148.273 * tfactors.T913
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
    get_partition_function(O16, tfactors, O16_pf, dO16_pf_dT);

    amrex::Real z_r = He4_pf * C12_pf;
    amrex::Real z_p = O16_pf;

    amrex::Real dz_r_dT = C12_pf * dHe4_pf_dT + He4_pf * dC12_pf_dT;
    amrex::Real dz_p_dT = dO16_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_F18_to_He4_N14_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // F18 --> He4 + N14

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10c
    ln_set_rate =  46.248979543913606 + -51.22930570025806 * tfactors.T9i + -36.2504 * tfactors.T913i
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
    ln_set_rate =  38.61457954391361 + -62.19490570025806 * tfactors.T9i + -5.6227 * tfactors.T913i;

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
    ln_set_rate =  24.91191754391361 + -56.38964570025806 * tfactors.T9i;

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


    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real F18_pf, dF18_pf_dT;
    // interpolating F18 partition function
    get_partition_function(F18, tfactors, F18_pf, dF18_pf_dT);

    amrex::Real N14_pf, dN14_pf_dT;
    // setting N14 partition function to 1.0 by default, independent of T
    N14_pf = 1.0_rt;
    dN14_pf_dT = 0.0_rt;

    amrex::Real z_r = He4_pf * N14_pf;
    amrex::Real z_p = F18_pf;

    amrex::Real dz_r_dT = N14_pf * dHe4_pf_dT + He4_pf * dN14_pf_dT;
    amrex::Real dz_p_dT = dF18_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ne20_to_He4_O16_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ne20 --> He4 + O16

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // co10c
    ln_set_rate =  34.26581585936365 + -67.65193003813836 * tfactors.T9i + -3.65925 * tfactors.T913
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
    ln_set_rate =  28.64304585936365 + -65.24613003813835 * tfactors.T9i;

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
    ln_set_rate =  48.66033585936365 + -54.887630038138354 * tfactors.T9i + -39.7262 * tfactors.T913i + -0.210799 * tfactors.T913
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


    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real O16_pf, dO16_pf_dT;
    // interpolating O16 partition function
    get_partition_function(O16, tfactors, O16_pf, dO16_pf_dT);

    amrex::Real Ne20_pf, dNe20_pf_dT;
    // interpolating Ne20 partition function
    get_partition_function(Ne20, tfactors, Ne20_pf, dNe20_pf_dT);

    amrex::Real z_r = He4_pf * O16_pf;
    amrex::Real z_p = Ne20_pf;

    amrex::Real dz_r_dT = O16_pf * dHe4_pf_dT + He4_pf * dO16_pf_dT;
    amrex::Real dz_p_dT = dNe20_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Na22_to_p_Ne21_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Na22 --> p + Ne21

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10c
    ln_set_rate =  -24.579038986172677 + -78.40599033331655 * tfactors.T9i;

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
    ln_set_rate =  42.145961013827325 + -78.20981033331655 * tfactors.T9i + -19.2096 * tfactors.T913i
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
    ln_set_rate =  -16.409838986172673 + -82.42366033331655 * tfactors.T9i + 21.1176 * tfactors.T913i + 34.0411 * tfactors.T913
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
    ln_set_rate =  24.833401013827324 + -79.60938033331655 * tfactors.T9i;

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


    amrex::Real Ne21_pf, dNe21_pf_dT;
    // interpolating Ne21 partition function
    get_partition_function(Ne21, tfactors, Ne21_pf, dNe21_pf_dT);

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Na22_pf, dNa22_pf_dT;
    // interpolating Na22 partition function
    get_partition_function(Na22, tfactors, Na22_pf, dNa22_pf_dT);

    amrex::Real z_r = p_pf * Ne21_pf;
    amrex::Real z_p = Na22_pf;

    amrex::Real dz_r_dT = Ne21_pf * dp_pf_dT + p_pf * dNe21_pf_dT;
    amrex::Real dz_p_dT = dNa22_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Na22_to_He4_F18_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Na22 --> He4 + F18

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // rpsmr
    ln_set_rate =  59.32234728275454 + -100.23588367240258 * tfactors.T9i + 18.8956 * tfactors.T913i + -65.6134 * tfactors.T913
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


    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real F18_pf, dF18_pf_dT;
    // interpolating F18 partition function
    get_partition_function(F18, tfactors, F18_pf, dF18_pf_dT);

    amrex::Real Na22_pf, dNa22_pf_dT;
    // interpolating Na22 partition function
    get_partition_function(Na22, tfactors, Na22_pf, dNa22_pf_dT);

    amrex::Real z_r = He4_pf * F18_pf;
    amrex::Real z_p = Na22_pf;

    amrex::Real dz_r_dT = F18_pf * dHe4_pf_dT + He4_pf * dF18_pf_dT;
    amrex::Real dz_p_dT = dNa22_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Mg24_to_p_Na23_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Mg24 --> p + Na23

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10c
    ln_set_rate =  43.93571176470678 + -135.6881490421346 * tfactors.T9i + -20.6428 * tfactors.T913i + 1.52954 * tfactors.T913
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
    ln_set_rate =  34.08761176470678 + -138.9684390421346 * tfactors.T9i + -0.360588 * tfactors.T913
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
    ln_set_rate =  20.002361764706784 + -137.3003390421346 * tfactors.T9i;

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


    amrex::Real Na23_pf, dNa23_pf_dT;
    // interpolating Na23 partition function
    get_partition_function(Na23, tfactors, Na23_pf, dNa23_pf_dT);

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Mg24_pf, dMg24_pf_dT;
    // interpolating Mg24 partition function
    get_partition_function(Mg24, tfactors, Mg24_pf, dMg24_pf_dT);

    amrex::Real z_r = p_pf * Na23_pf;
    amrex::Real z_p = Mg24_pf;

    amrex::Real dz_r_dT = Na23_pf * dp_pf_dT + p_pf * dNa23_pf_dT;
    amrex::Real dz_p_dT = dMg24_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Mg24_to_He4_Ne20_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Mg24 --> He4 + Ne20

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10c
    ln_set_rate =  -13.886931148855957 + -110.62012330714883 * tfactors.T9i;

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
    ln_set_rate =  49.324368851144044 + -108.11407330714883 * tfactors.T9i + -46.2525 * tfactors.T913i + 5.58901 * tfactors.T913
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
    ln_set_rate =  16.02029885114404 + -120.89497330714883 * tfactors.T9i + 16.9229 * tfactors.T913
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
    ln_set_rate =  26.801638851144045 + -117.33433330714882 * tfactors.T9i;

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


    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Mg24_pf, dMg24_pf_dT;
    // interpolating Mg24 partition function
    get_partition_function(Mg24, tfactors, Mg24_pf, dMg24_pf_dT);

    amrex::Real Ne20_pf, dNe20_pf_dT;
    // interpolating Ne20 partition function
    get_partition_function(Ne20, tfactors, Ne20_pf, dNe20_pf_dT);

    amrex::Real z_r = He4_pf * Ne20_pf;
    amrex::Real z_p = Mg24_pf;

    amrex::Real dz_r_dT = Ne20_pf * dHe4_pf_dT + He4_pf * dNe20_pf_dT;
    amrex::Real dz_p_dT = dMg24_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Si28_to_p_Al27_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Si28 --> p + Al27

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10c
    ln_set_rate =  11.776564828186824 + -136.3434628922348 * tfactors.T9i + 23.8634 * tfactors.T913
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
    ln_set_rate =  46.54946482818683 + -134.4395028922348 * tfactors.T9i + -23.2205 * tfactors.T913i
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


    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Al27_pf, dAl27_pf_dT;
    // interpolating Al27 partition function
    get_partition_function(Al27, tfactors, Al27_pf, dAl27_pf_dT);

    amrex::Real Si28_pf, dSi28_pf_dT;
    // interpolating Si28 partition function
    get_partition_function(Si28, tfactors, Si28_pf, dSi28_pf_dT);

    amrex::Real z_r = p_pf * Al27_pf;
    amrex::Real z_p = Si28_pf;

    amrex::Real dz_r_dT = Al27_pf * dp_pf_dT + p_pf * dAl27_pf_dT;
    amrex::Real dz_p_dT = dSi28_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Si28_to_He4_Mg24_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Si28 --> He4 + Mg24

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // st08r
    ln_set_rate =  -25.688574833405916 + -128.69270892750794 * tfactors.T9i + 21.3721 * tfactors.T913i + 37.7649 * tfactors.T913
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
    ln_set_rate =  32.90059516659409 + -131.48850892750792 * tfactors.T9i;

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


    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Mg24_pf, dMg24_pf_dT;
    // interpolating Mg24 partition function
    get_partition_function(Mg24, tfactors, Mg24_pf, dMg24_pf_dT);

    amrex::Real Si28_pf, dSi28_pf_dT;
    // interpolating Si28 partition function
    get_partition_function(Si28, tfactors, Si28_pf, dSi28_pf_dT);

    amrex::Real z_r = He4_pf * Mg24_pf;
    amrex::Real z_p = Si28_pf;

    amrex::Real dz_r_dT = Mg24_pf * dHe4_pf_dT + He4_pf * dMg24_pf_dT;
    amrex::Real dz_p_dT = dSi28_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_S32_to_p_P31_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // S32 --> p + P31

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10c
    ln_set_rate =  25.172836958303158 + -106.63693563716609 * tfactors.T9i + 8.09341 * tfactors.T913
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
    ln_set_rate =  21.682890958303158 + -105.11947563716609 * tfactors.T9i;

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
    ln_set_rate =  43.610880958303156 + -102.85989563716609 * tfactors.T9i + -25.3278 * tfactors.T913i + 6.4931 * tfactors.T913
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


    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real S32_pf, dS32_pf_dT;
    // interpolating S32 partition function
    get_partition_function(S32, tfactors, S32_pf, dS32_pf_dT);

    amrex::Real P31_pf, dP31_pf_dT;
    // interpolating P31 partition function
    get_partition_function(P31, tfactors, P31_pf, dP31_pf_dT);

    amrex::Real z_r = p_pf * P31_pf;
    amrex::Real z_p = S32_pf;

    amrex::Real dz_r_dT = P31_pf * dp_pf_dT + p_pf * dP31_pf_dT;
    amrex::Real dz_p_dT = dS32_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_S32_to_He4_Si28_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // S32 --> He4 + Si28

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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


    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real S32_pf, dS32_pf_dT;
    // interpolating S32 partition function
    get_partition_function(S32, tfactors, S32_pf, dS32_pf_dT);

    amrex::Real Si28_pf, dSi28_pf_dT;
    // interpolating Si28 partition function
    get_partition_function(Si28, tfactors, Si28_pf, dSi28_pf_dT);

    amrex::Real z_r = He4_pf * Si28_pf;
    amrex::Real z_p = S32_pf;

    amrex::Real dz_r_dT = Si28_pf * dHe4_pf_dT + He4_pf * dSi28_pf_dT;
    amrex::Real dz_p_dT = dS32_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Fe52_to_p_Mn51_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe52 --> p + Mn51

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  61.72798916565748 + -85.6326403498911 * tfactors.T9i + -36.1825 * tfactors.T913i + 0.873042 * tfactors.T913
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


    amrex::Real Mn51_pf, dMn51_pf_dT;
    // interpolating Mn51 partition function
    get_partition_function(Mn51, tfactors, Mn51_pf, dMn51_pf_dT);

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Fe52_pf, dFe52_pf_dT;
    // interpolating Fe52 partition function
    get_partition_function(Fe52, tfactors, Fe52_pf, dFe52_pf_dT);

    amrex::Real z_r = p_pf * Mn51_pf;
    amrex::Real z_p = Fe52_pf;

    amrex::Real dz_r_dT = Mn51_pf * dp_pf_dT + p_pf * dMn51_pf_dT;
    amrex::Real dz_p_dT = dFe52_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Fe52_to_He4_Cr48_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe52 --> He4 + Cr48

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  90.14738712482466 + -92.10912191363732 * tfactors.T9i + -86.7459 * tfactors.T913i + -9.79373 * tfactors.T913
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


    amrex::Real Cr48_pf, dCr48_pf_dT;
    // interpolating Cr48 partition function
    get_partition_function(Cr48, tfactors, Cr48_pf, dCr48_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Fe52_pf, dFe52_pf_dT;
    // interpolating Fe52 partition function
    get_partition_function(Fe52, tfactors, Fe52_pf, dFe52_pf_dT);

    amrex::Real z_r = He4_pf * Cr48_pf;
    amrex::Real z_p = Fe52_pf;

    amrex::Real dz_r_dT = Cr48_pf * dHe4_pf_dT + He4_pf * dCr48_pf_dT;
    amrex::Real dz_p_dT = dFe52_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Co55_to_He4_Mn51_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co55 --> He4 + Mn51

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  90.61300778512097 + -95.28620688500884 * tfactors.T9i + -89.274 * tfactors.T913i + -10.4373 * tfactors.T913
                         + 1.00492 * tfactors.T9 + -0.125548 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  95.28620688500884 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -89.274 * tfactors.T943i + (1.0/3.0) * -10.4373 * tfactors.T923i
                                  + 1.00492 + (5.0/3.0) * -0.125548 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real Mn51_pf, dMn51_pf_dT;
    // interpolating Mn51 partition function
    get_partition_function(Mn51, tfactors, Mn51_pf, dMn51_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Co55_pf, dCo55_pf_dT;
    // interpolating Co55 partition function
    get_partition_function(Co55, tfactors, Co55_pf, dCo55_pf_dT);

    amrex::Real z_r = He4_pf * Mn51_pf;
    amrex::Real z_p = Co55_pf;

    amrex::Real dz_r_dT = Mn51_pf * dHe4_pf_dT + He4_pf * dMn51_pf_dT;
    amrex::Real dz_p_dT = dCo55_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ni56_to_p_Co55_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni56 --> p + Co55

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  63.131770608640906 + -83.14741674893808 * tfactors.T9i + -38.1053 * tfactors.T913i + -0.210947 * tfactors.T913
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


    amrex::Real Ni56_pf, dNi56_pf_dT;
    // interpolating Ni56 partition function
    get_partition_function(Ni56, tfactors, Ni56_pf, dNi56_pf_dT);

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Co55_pf, dCo55_pf_dT;
    // interpolating Co55 partition function
    get_partition_function(Co55, tfactors, Co55_pf, dCo55_pf_dT);

    amrex::Real z_r = p_pf * Co55_pf;
    amrex::Real z_p = Ni56_pf;

    amrex::Real dz_r_dT = Co55_pf * dp_pf_dT + p_pf * dCo55_pf_dT;
    amrex::Real dz_p_dT = dNi56_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ni56_to_He4_Fe52_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni56 --> He4 + Fe52

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  91.62258922810439 + -92.801099329237 * tfactors.T9i + -91.6819 * tfactors.T913i + -9.51885 * tfactors.T913
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


    amrex::Real Ni56_pf, dNi56_pf_dT;
    // interpolating Ni56 partition function
    get_partition_function(Ni56, tfactors, Ni56_pf, dNi56_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Fe52_pf, dFe52_pf_dT;
    // interpolating Fe52 partition function
    get_partition_function(Fe52, tfactors, Fe52_pf, dFe52_pf_dT);

    amrex::Real z_r = He4_pf * Fe52_pf;
    amrex::Real z_p = Ni56_pf;

    amrex::Real dz_r_dT = Fe52_pf * dHe4_pf_dT + He4_pf * dFe52_pf_dT;
    amrex::Real dz_p_dT = dNi56_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Cu59_to_He4_Co55_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Cu59 --> He4 + Co55

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  90.97199197888506 + -55.17460977128912 * tfactors.T9i + -94.1404 * tfactors.T913i + -2.62786 * tfactors.T913
                         + -2.12066 * tfactors.T9 + 0.237999 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  55.17460977128912 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -94.1404 * tfactors.T943i + (1.0/3.0) * -2.62786 * tfactors.T923i
                                  + -2.12066 + (5.0/3.0) * 0.237999 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real Cu59_pf, dCu59_pf_dT;
    // interpolating Cu59 partition function
    get_partition_function(Cu59, tfactors, Cu59_pf, dCu59_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Co55_pf, dCo55_pf_dT;
    // interpolating Co55 partition function
    get_partition_function(Co55, tfactors, Co55_pf, dCo55_pf_dT);

    amrex::Real z_r = He4_pf * Co55_pf;
    amrex::Real z_p = Cu59_pf;

    amrex::Real dz_r_dT = Co55_pf * dHe4_pf_dT + He4_pf * dCo55_pf_dT;
    amrex::Real dz_p_dT = dCu59_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Zn60_to_p_Cu59_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Zn60 --> p + Cu59

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  62.2928405088604 + -59.40944656945774 * tfactors.T9i + -39.9808 * tfactors.T913i + 1.17311 * tfactors.T913
                         + -2.90486 * tfactors.T9 + 0.339644 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  59.40944656945774 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -39.9808 * tfactors.T943i + (1.0/3.0) * 1.17311 * tfactors.T923i
                                  + -2.90486 + (5.0/3.0) * 0.339644 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real Cu59_pf, dCu59_pf_dT;
    // interpolating Cu59 partition function
    get_partition_function(Cu59, tfactors, Cu59_pf, dCu59_pf_dT);

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Zn60_pf, dZn60_pf_dT;
    // interpolating Zn60 partition function
    get_partition_function(Zn60, tfactors, Zn60_pf, dZn60_pf_dT);

    amrex::Real z_r = p_pf * Cu59_pf;
    amrex::Real z_p = Zn60_pf;

    amrex::Real dz_r_dT = Cu59_pf * dp_pf_dT + p_pf * dCu59_pf_dT;
    amrex::Real dz_p_dT = dZn60_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Zn60_to_He4_Ni56_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Zn60 --> He4 + Ni56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  86.06186187910455 + -31.43675563699001 * tfactors.T9i + -96.4898 * tfactors.T913i + 6.47209 * tfactors.T913
                         + -5.2029 * tfactors.T9 + 0.533391 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  31.43675563699001 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -96.4898 * tfactors.T943i + (1.0/3.0) * 6.47209 * tfactors.T923i
                                  + -5.2029 + (5.0/3.0) * 0.533391 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real Ni56_pf, dNi56_pf_dT;
    // interpolating Ni56 partition function
    get_partition_function(Ni56, tfactors, Ni56_pf, dNi56_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Zn60_pf, dZn60_pf_dT;
    // interpolating Zn60 partition function
    get_partition_function(Zn60, tfactors, Zn60_pf, dZn60_pf_dT);

    amrex::Real z_r = He4_pf * Ni56_pf;
    amrex::Real z_p = Zn60_pf;

    amrex::Real dz_r_dT = Ni56_pf * dHe4_pf_dT + He4_pf * dNi56_pf_dT;
    amrex::Real dz_p_dT = dZn60_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_C12_to_He4_He4_He4_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // C12 --> He4 + He4 + He4

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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

    amrex::Real dz_r_dT = dHe4_pf_dT + dHe4_pf_dT + dHe4_pf_dT;
    amrex::Real dz_p_dT = dC12_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_O16_to_He4_N13_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // O16 + p --> He4 + N13

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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


    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real O16_pf, dO16_pf_dT;
    // interpolating O16 partition function
    get_partition_function(O16, tfactors, O16_pf, dO16_pf_dT);

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

    amrex::Real dz_r_dT = N13_pf * dHe4_pf_dT + He4_pf * dN13_pf_dT;
    amrex::Real dz_p_dT = O16_pf * dp_pf_dT + p_pf * dO16_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Ne20_to_p_Na23_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ne20 + He4 --> p + Na23

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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


    amrex::Real Na23_pf, dNa23_pf_dT;
    // interpolating Na23 partition function
    get_partition_function(Na23, tfactors, Na23_pf, dNa23_pf_dT);

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Ne20_pf, dNe20_pf_dT;
    // interpolating Ne20 partition function
    get_partition_function(Ne20, tfactors, Ne20_pf, dNe20_pf_dT);

    amrex::Real z_r = p_pf * Na23_pf;
    amrex::Real z_p = He4_pf * Ne20_pf;

    amrex::Real dz_r_dT = Na23_pf * dp_pf_dT + p_pf * dNa23_pf_dT;
    amrex::Real dz_p_dT = Ne20_pf * dHe4_pf_dT + He4_pf * dNe20_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Ne21_to_He4_F18_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ne21 + p --> He4 + F18

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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


    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Ne21_pf, dNe21_pf_dT;
    // interpolating Ne21 partition function
    get_partition_function(Ne21, tfactors, Ne21_pf, dNe21_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real F18_pf, dF18_pf_dT;
    // interpolating F18 partition function
    get_partition_function(F18, tfactors, F18_pf, dF18_pf_dT);

    amrex::Real z_r = He4_pf * F18_pf;
    amrex::Real z_p = p_pf * Ne21_pf;

    amrex::Real dz_r_dT = F18_pf * dHe4_pf_dT + He4_pf * dF18_pf_dT;
    amrex::Real dz_p_dT = Ne21_pf * dp_pf_dT + p_pf * dNe21_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Mg24_to_p_Al27_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Mg24 + He4 --> p + Al27

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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


    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Al27_pf, dAl27_pf_dT;
    // interpolating Al27 partition function
    get_partition_function(Al27, tfactors, Al27_pf, dAl27_pf_dT);

    amrex::Real Mg24_pf, dMg24_pf_dT;
    // interpolating Mg24 partition function
    get_partition_function(Mg24, tfactors, Mg24_pf, dMg24_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real z_r = p_pf * Al27_pf;
    amrex::Real z_p = He4_pf * Mg24_pf;

    amrex::Real dz_r_dT = Al27_pf * dp_pf_dT + p_pf * dAl27_pf_dT;
    amrex::Real dz_p_dT = Mg24_pf * dHe4_pf_dT + He4_pf * dMg24_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Si28_to_p_P31_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Si28 + He4 --> p + P31

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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


    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real P31_pf, dP31_pf_dT;
    // interpolating P31 partition function
    get_partition_function(P31, tfactors, P31_pf, dP31_pf_dT);

    amrex::Real Si28_pf, dSi28_pf_dT;
    // interpolating Si28 partition function
    get_partition_function(Si28, tfactors, Si28_pf, dSi28_pf_dT);

    amrex::Real z_r = p_pf * P31_pf;
    amrex::Real z_p = He4_pf * Si28_pf;

    amrex::Real dz_r_dT = P31_pf * dp_pf_dT + p_pf * dP31_pf_dT;
    amrex::Real dz_p_dT = Si28_pf * dHe4_pf_dT + He4_pf * dSi28_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Mn51_to_He4_Cr48_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Mn51 + p --> He4 + Cr48

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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


    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Cr48_pf, dCr48_pf_dT;
    // interpolating Cr48 partition function
    get_partition_function(Cr48, tfactors, Cr48_pf, dCr48_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Mn51_pf, dMn51_pf_dT;
    // interpolating Mn51 partition function
    get_partition_function(Mn51, tfactors, Mn51_pf, dMn51_pf_dT);

    amrex::Real z_r = He4_pf * Cr48_pf;
    amrex::Real z_p = p_pf * Mn51_pf;

    amrex::Real dz_r_dT = Cr48_pf * dHe4_pf_dT + He4_pf * dCr48_pf_dT;
    amrex::Real dz_p_dT = Mn51_pf * dp_pf_dT + p_pf * dMn51_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Co55_to_He4_Fe52_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co55 + p --> He4 + Fe52

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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


    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Co55_pf, dCo55_pf_dT;
    // interpolating Co55 partition function
    get_partition_function(Co55, tfactors, Co55_pf, dCo55_pf_dT);

    amrex::Real Fe52_pf, dFe52_pf_dT;
    // interpolating Fe52 partition function
    get_partition_function(Fe52, tfactors, Fe52_pf, dFe52_pf_dT);

    amrex::Real z_r = He4_pf * Fe52_pf;
    amrex::Real z_p = p_pf * Co55_pf;

    amrex::Real dz_r_dT = Fe52_pf * dHe4_pf_dT + He4_pf * dFe52_pf_dT;
    amrex::Real dz_p_dT = Co55_pf * dp_pf_dT + p_pf * dCo55_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Ni56_to_p_Cu59_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni56 + He4 --> p + Cu59

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  20.847078629755856 + -27.97280697764895 * tfactors.T9i + -39.9808 * tfactors.T913i + 13.6773 * tfactors.T913
                         + -3.76429 * tfactors.T9 + 0.438096 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  27.97280697764895 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -39.9808 * tfactors.T943i + (1.0/3.0) * 13.6773 * tfactors.T923i
                                  + -3.76429 + (5.0/3.0) * 0.438096 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real Cu59_pf, dCu59_pf_dT;
    // interpolating Cu59 partition function
    get_partition_function(Cu59, tfactors, Cu59_pf, dCu59_pf_dT);

    amrex::Real Ni56_pf, dNi56_pf_dT;
    // interpolating Ni56 partition function
    get_partition_function(Ni56, tfactors, Ni56_pf, dNi56_pf_dT);

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real z_r = p_pf * Cu59_pf;
    amrex::Real z_p = He4_pf * Ni56_pf;

    amrex::Real dz_r_dT = Cu59_pf * dp_pf_dT + p_pf * dCu59_pf_dT;
    amrex::Real dz_p_dT = Ni56_pf * dHe4_pf_dT + He4_pf * dNi56_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Fe53_to_n_Fe52_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe53 --> n + Fe52

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  37.48624299107921 + -123.99427613085159 * tfactors.T9i + -0.344319 * tfactors.T913
                         + 0.178277 * tfactors.T9 + -0.0334326 * tfactors.T953 + 1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  123.99427613085159 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.344319 * tfactors.T923i
                                  + 0.178277 + (5.0/3.0) * -0.0334326 * tfactors.T923 + 1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real Fe52_pf, dFe52_pf_dT;
    // interpolating Fe52 partition function
    get_partition_function(Fe52, tfactors, Fe52_pf, dFe52_pf_dT);

    amrex::Real Fe53_pf, dFe53_pf_dT;
    // interpolating Fe53 partition function
    get_partition_function(Fe53, tfactors, Fe53_pf, dFe53_pf_dT);

    amrex::Real z_r = n_pf * Fe52_pf;
    amrex::Real z_p = Fe53_pf;

    amrex::Real dz_r_dT = Fe52_pf * dn_pf_dT + n_pf * dFe52_pf_dT;
    amrex::Real dz_p_dT = dFe53_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Fe54_to_n_Fe53_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe54 --> n + Fe53

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  42.210560167376684 + -155.25104569177333 * tfactors.T9i + -1.10421 * tfactors.T913
                         + 0.379905 * tfactors.T9 + -0.0581878 * tfactors.T953 + 1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  155.25104569177333 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.10421 * tfactors.T923i
                                  + 0.379905 + (5.0/3.0) * -0.0581878 * tfactors.T923 + 1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real Fe53_pf, dFe53_pf_dT;
    // interpolating Fe53 partition function
    get_partition_function(Fe53, tfactors, Fe53_pf, dFe53_pf_dT);

    amrex::Real Fe54_pf, dFe54_pf_dT;
    // interpolating Fe54 partition function
    get_partition_function(Fe54, tfactors, Fe54_pf, dFe54_pf_dT);

    amrex::Real z_r = n_pf * Fe53_pf;
    amrex::Real z_p = Fe54_pf;

    amrex::Real dz_r_dT = Fe53_pf * dn_pf_dT + n_pf * dFe53_pf_dT;
    amrex::Real dz_p_dT = dFe54_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Fe55_to_n_Fe54_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe55 --> n + Fe54

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ks03 
    ln_set_rate =  21.483298756092896 + -107.84311925142697 * tfactors.T9i + -8.66617 * tfactors.T913i + 26.4472 * tfactors.T913
                         + -1.9222 * tfactors.T9 + 0.0986404 * tfactors.T953 + -8.28317 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  107.84311925142697 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -8.66617 * tfactors.T943i + (1.0/3.0) * 26.4472 * tfactors.T923i
                                  + -1.9222 + (5.0/3.0) * 0.0986404 * tfactors.T923 + -8.28317 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real Fe55_pf, dFe55_pf_dT;
    // interpolating Fe55 partition function
    get_partition_function(Fe55, tfactors, Fe55_pf, dFe55_pf_dT);

    amrex::Real Fe54_pf, dFe54_pf_dT;
    // interpolating Fe54 partition function
    get_partition_function(Fe54, tfactors, Fe54_pf, dFe54_pf_dT);

    amrex::Real z_r = n_pf * Fe54_pf;
    amrex::Real z_p = Fe55_pf;

    amrex::Real dz_r_dT = Fe54_pf * dn_pf_dT + n_pf * dFe54_pf_dT;
    amrex::Real dz_p_dT = dFe55_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Fe56_to_n_Fe55_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe56 --> n + Fe55

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ks03 
    ln_set_rate =  46.785223428080954 + -130.0232339465001 * tfactors.T9i + 8.06062 * tfactors.T913i + -14.4809 * tfactors.T913
                         + 0.94252 * tfactors.T9 + -0.0776007 * tfactors.T953 + 7.97093 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  130.0232339465001 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 8.06062 * tfactors.T943i + (1.0/3.0) * -14.4809 * tfactors.T923i
                                  + 0.94252 + (5.0/3.0) * -0.0776007 * tfactors.T923 + 7.97093 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real Fe56_pf, dFe56_pf_dT;
    // interpolating Fe56 partition function
    get_partition_function(Fe56, tfactors, Fe56_pf, dFe56_pf_dT);

    amrex::Real Fe55_pf, dFe55_pf_dT;
    // interpolating Fe55 partition function
    get_partition_function(Fe55, tfactors, Fe55_pf, dFe55_pf_dT);

    amrex::Real z_r = n_pf * Fe55_pf;
    amrex::Real z_p = Fe56_pf;

    amrex::Real dz_r_dT = Fe55_pf * dn_pf_dT + n_pf * dFe55_pf_dT;
    amrex::Real dz_p_dT = dFe56_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Fe56_to_p_Mn55_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe56 --> p + Mn55

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  61.659288536189116 + -118.17693119642051 * tfactors.T9i + -36.1995 * tfactors.T913i + 1.15346 * tfactors.T913
                         + -1.94437 * tfactors.T9 + 0.218429 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  118.17693119642051 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -36.1995 * tfactors.T943i + (1.0/3.0) * 1.15346 * tfactors.T923i
                                  + -1.94437 + (5.0/3.0) * 0.218429 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Fe56_pf, dFe56_pf_dT;
    // interpolating Fe56 partition function
    get_partition_function(Fe56, tfactors, Fe56_pf, dFe56_pf_dT);

    amrex::Real Mn55_pf, dMn55_pf_dT;
    // interpolating Mn55 partition function
    get_partition_function(Mn55, tfactors, Mn55_pf, dMn55_pf_dT);

    amrex::Real z_r = p_pf * Mn55_pf;
    amrex::Real z_p = Fe56_pf;

    amrex::Real dz_r_dT = Mn55_pf * dp_pf_dT + p_pf * dMn55_pf_dT;
    amrex::Real dz_p_dT = dFe56_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Co55_to_p_Fe54_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co55 --> p + Fe54

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  57.82919157553295 + -58.76609208478817 * tfactors.T9i + -37.1544 * tfactors.T913i + 0.950364 * tfactors.T913
                         + -1.77529 * tfactors.T9 + 0.198562 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  58.76609208478817 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -37.1544 * tfactors.T943i + (1.0/3.0) * 0.950364 * tfactors.T923i
                                  + -1.77529 + (5.0/3.0) * 0.198562 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Co55_pf, dCo55_pf_dT;
    // interpolating Co55 partition function
    get_partition_function(Co55, tfactors, Co55_pf, dCo55_pf_dT);

    amrex::Real Fe54_pf, dFe54_pf_dT;
    // interpolating Fe54 partition function
    get_partition_function(Fe54, tfactors, Fe54_pf, dFe54_pf_dT);

    amrex::Real z_r = p_pf * Fe54_pf;
    amrex::Real z_p = Co55_pf;

    amrex::Real dz_r_dT = Fe54_pf * dp_pf_dT + p_pf * dFe54_pf_dT;
    amrex::Real dz_p_dT = dCo55_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Co56_to_n_Co55_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co56 --> n + Co55

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  40.62594603130468 + -117.00951667337291 * tfactors.T9i + -1.86357 * tfactors.T913
                         + 0.616591 * tfactors.T9 + -0.0839313 * tfactors.T953 + 1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  117.00951667337291 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.86357 * tfactors.T923i
                                  + 0.616591 + (5.0/3.0) * -0.0839313 * tfactors.T923 + 1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real Co55_pf, dCo55_pf_dT;
    // interpolating Co55 partition function
    get_partition_function(Co55, tfactors, Co55_pf, dCo55_pf_dT);

    amrex::Real Co56_pf, dCo56_pf_dT;
    // interpolating Co56 partition function
    get_partition_function(Co56, tfactors, Co56_pf, dCo56_pf_dT);

    amrex::Real z_r = n_pf * Co55_pf;
    amrex::Real z_p = Co56_pf;

    amrex::Real dz_r_dT = Co55_pf * dn_pf_dT + n_pf * dCo55_pf_dT;
    amrex::Real dz_p_dT = dCo56_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Co56_to_p_Fe55_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co56 --> p + Fe55

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  58.601098850744734 + -67.87436231336505 * tfactors.T9i + -37.1585 * tfactors.T913i + 1.66198 * tfactors.T913
                         + -1.60842 * tfactors.T9 + 0.148916 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  67.87436231336505 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -37.1585 * tfactors.T943i + (1.0/3.0) * 1.66198 * tfactors.T923i
                                  + -1.60842 + (5.0/3.0) * 0.148916 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real Fe55_pf, dFe55_pf_dT;
    // interpolating Fe55 partition function
    get_partition_function(Fe55, tfactors, Fe55_pf, dFe55_pf_dT);

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Co56_pf, dCo56_pf_dT;
    // interpolating Co56 partition function
    get_partition_function(Co56, tfactors, Co56_pf, dCo56_pf_dT);

    amrex::Real z_r = p_pf * Fe55_pf;
    amrex::Real z_p = Co56_pf;

    amrex::Real dz_r_dT = Fe55_pf * dp_pf_dT + p_pf * dFe55_pf_dT;
    amrex::Real dz_p_dT = dCo56_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Co57_to_n_Co56_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co57 --> n + Co56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  41.15219049522237 + -132.01531905660212 * tfactors.T9i + -1.37855 * tfactors.T913
                         + 0.299896 * tfactors.T9 + -0.04382 * tfactors.T953 + 1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  132.01531905660212 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.37855 * tfactors.T923i
                                  + 0.299896 + (5.0/3.0) * -0.04382 * tfactors.T923 + 1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real Co56_pf, dCo56_pf_dT;
    // interpolating Co56 partition function
    get_partition_function(Co56, tfactors, Co56_pf, dCo56_pf_dT);

    amrex::Real Co57_pf, dCo57_pf_dT;
    // interpolating Co57 partition function
    get_partition_function(Co57, tfactors, Co57_pf, dCo57_pf_dT);

    amrex::Real z_r = n_pf * Co56_pf;
    amrex::Real z_p = Co57_pf;

    amrex::Real dz_r_dT = Co56_pf * dn_pf_dT + n_pf * dCo56_pf_dT;
    amrex::Real dz_p_dT = dCo57_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Co57_to_p_Fe56_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co57 --> p + Fe56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  57.66626591788614 + -69.94994642462048 * tfactors.T9i + -37.1625 * tfactors.T913i + 1.06776 * tfactors.T913
                         + -1.31689 * tfactors.T9 + 0.122089 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  69.94994642462048 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -37.1625 * tfactors.T943i + (1.0/3.0) * 1.06776 * tfactors.T923i
                                  + -1.31689 + (5.0/3.0) * 0.122089 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Fe56_pf, dFe56_pf_dT;
    // interpolating Fe56 partition function
    get_partition_function(Fe56, tfactors, Fe56_pf, dFe56_pf_dT);

    amrex::Real Co57_pf, dCo57_pf_dT;
    // interpolating Co57 partition function
    get_partition_function(Co57, tfactors, Co57_pf, dCo57_pf_dT);

    amrex::Real z_r = p_pf * Fe56_pf;
    amrex::Real z_p = Co57_pf;

    amrex::Real dz_r_dT = Fe56_pf * dp_pf_dT + p_pf * dFe56_pf_dT;
    amrex::Real dz_p_dT = dCo57_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ni57_to_n_Ni56_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni57 --> n + Ni56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  38.36941309844609 + -118.94166894064354 * tfactors.T9i + -1.19665 * tfactors.T913
                         + 0.507179 * tfactors.T9 + -0.074604 * tfactors.T953 + 1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  118.94166894064354 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.19665 * tfactors.T923i
                                  + 0.507179 + (5.0/3.0) * -0.074604 * tfactors.T923 + 1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real Ni56_pf, dNi56_pf_dT;
    // interpolating Ni56 partition function
    get_partition_function(Ni56, tfactors, Ni56_pf, dNi56_pf_dT);

    amrex::Real Ni57_pf, dNi57_pf_dT;
    // interpolating Ni57 partition function
    get_partition_function(Ni57, tfactors, Ni57_pf, dNi57_pf_dT);

    amrex::Real z_r = n_pf * Ni56_pf;
    amrex::Real z_p = Ni57_pf;

    amrex::Real dz_r_dT = Ni56_pf * dn_pf_dT + n_pf * dNi56_pf_dT;
    amrex::Real dz_p_dT = dNi57_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ni57_to_p_Co56_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni57 --> p + Co56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  60.08383767578231 + -85.07956901620872 * tfactors.T9i + -38.1094 * tfactors.T913i + 2.57091 * tfactors.T913
                         + -2.07795 * tfactors.T9 + 0.20757 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  85.07956901620872 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -38.1094 * tfactors.T943i + (1.0/3.0) * 2.57091 * tfactors.T923i
                                  + -2.07795 + (5.0/3.0) * 0.20757 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real Ni57_pf, dNi57_pf_dT;
    // interpolating Ni57 partition function
    get_partition_function(Ni57, tfactors, Ni57_pf, dNi57_pf_dT);

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Co56_pf, dCo56_pf_dT;
    // interpolating Co56 partition function
    get_partition_function(Co56, tfactors, Co56_pf, dCo56_pf_dT);

    amrex::Real z_r = p_pf * Co56_pf;
    amrex::Real z_p = Ni57_pf;

    amrex::Real dz_r_dT = Co56_pf * dp_pf_dT + p_pf * dCo56_pf_dT;
    amrex::Real dz_p_dT = dNi57_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ni57_to_He4_Fe53_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni57 --> He4 + Fe53

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  94.73965933547126 + -87.74826004866655 * tfactors.T9i + -91.7231 * tfactors.T913i + -10.1816 * tfactors.T913
                         + -0.0406303 * tfactors.T9 + 0.0345056 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  87.74826004866655 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -91.7231 * tfactors.T943i + (1.0/3.0) * -10.1816 * tfactors.T923i
                                  + -0.0406303 + (5.0/3.0) * 0.0345056 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real Ni57_pf, dNi57_pf_dT;
    // interpolating Ni57 partition function
    get_partition_function(Ni57, tfactors, Ni57_pf, dNi57_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Fe53_pf, dFe53_pf_dT;
    // interpolating Fe53 partition function
    get_partition_function(Fe53, tfactors, Fe53_pf, dFe53_pf_dT);

    amrex::Real z_r = He4_pf * Fe53_pf;
    amrex::Real z_p = Ni57_pf;

    amrex::Real dz_r_dT = Fe53_pf * dHe4_pf_dT + He4_pf * dFe53_pf_dT;
    amrex::Real dz_p_dT = dNi57_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ni58_to_n_Ni57_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni58 --> n + Ni57

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  42.339063572267165 + -141.77239789336582 * tfactors.T9i + -1.90814 * tfactors.T913
                         + 0.493188 * tfactors.T9 + -0.0684633 * tfactors.T953 + 1.5 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  141.77239789336582 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.90814 * tfactors.T923i
                                  + 0.493188 + (5.0/3.0) * -0.0684633 * tfactors.T923 + 1.5 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real Ni57_pf, dNi57_pf_dT;
    // interpolating Ni57 partition function
    get_partition_function(Ni57, tfactors, Ni57_pf, dNi57_pf_dT);

    amrex::Real Ni58_pf, dNi58_pf_dT;
    // interpolating Ni58 partition function
    get_partition_function(Ni58, tfactors, Ni58_pf, dNi58_pf_dT);

    amrex::Real z_r = n_pf * Ni57_pf;
    amrex::Real z_p = Ni58_pf;

    amrex::Real dz_r_dT = Ni57_pf * dn_pf_dT + n_pf * dNi57_pf_dT;
    amrex::Real dz_p_dT = dNi58_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ni58_to_p_Co57_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni58 --> p + Co57

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  61.775010752827114 + -94.83746016924091 * tfactors.T9i + -38.1133 * tfactors.T913i + 1.77414 * tfactors.T913
                         + -1.48268 * tfactors.T9 + 0.121073 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  94.83746016924091 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -38.1133 * tfactors.T943i + (1.0/3.0) * 1.77414 * tfactors.T923i
                                  + -1.48268 + (5.0/3.0) * 0.121073 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Ni58_pf, dNi58_pf_dT;
    // interpolating Ni58 partition function
    get_partition_function(Ni58, tfactors, Ni58_pf, dNi58_pf_dT);

    amrex::Real Co57_pf, dCo57_pf_dT;
    // interpolating Co57 partition function
    get_partition_function(Co57, tfactors, Co57_pf, dCo57_pf_dT);

    amrex::Real z_r = p_pf * Co57_pf;
    amrex::Real z_p = Ni58_pf;

    amrex::Real dz_r_dT = Co57_pf * dp_pf_dT + p_pf * dCo57_pf_dT;
    amrex::Real dz_p_dT = dNi58_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ni58_to_He4_Fe54_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni58 --> He4 + Fe54

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  85.23266274036175 + -74.27065665689001 * tfactors.T9i + -91.7628 * tfactors.T913i + 4.23027 * tfactors.T913
                         + -3.31305 * tfactors.T9 + 0.271293 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  74.27065665689001 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -91.7628 * tfactors.T943i + (1.0/3.0) * 4.23027 * tfactors.T923i
                                  + -3.31305 + (5.0/3.0) * 0.271293 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Ni58_pf, dNi58_pf_dT;
    // interpolating Ni58 partition function
    get_partition_function(Ni58, tfactors, Ni58_pf, dNi58_pf_dT);

    amrex::Real Fe54_pf, dFe54_pf_dT;
    // interpolating Fe54 partition function
    get_partition_function(Fe54, tfactors, Fe54_pf, dFe54_pf_dT);

    amrex::Real z_r = He4_pf * Fe54_pf;
    amrex::Real z_p = Ni58_pf;

    amrex::Real dz_r_dT = Fe54_pf * dHe4_pf_dT + He4_pf * dFe54_pf_dT;
    amrex::Real dz_p_dT = dNi58_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Cu59_to_p_Ni58_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Cu59 --> p + Ni58

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // si13n
    ln_set_rate =  58.96002081405625 + -39.670045199187285 * tfactors.T9i + -39.0526 * tfactors.T913i + 1.0436 * tfactors.T913
                         + -2.10834 * tfactors.T9 + 0.239559 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  39.670045199187285 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -39.0526 * tfactors.T943i + (1.0/3.0) * 1.0436 * tfactors.T923i
                                  + -2.10834 + (5.0/3.0) * 0.239559 * tfactors.T923 + 0.833333 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real Cu59_pf, dCu59_pf_dT;
    // interpolating Cu59 partition function
    get_partition_function(Cu59, tfactors, Cu59_pf, dCu59_pf_dT);

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Ni58_pf, dNi58_pf_dT;
    // interpolating Ni58 partition function
    get_partition_function(Ni58, tfactors, Ni58_pf, dNi58_pf_dT);

    amrex::Real z_r = p_pf * Ni58_pf;
    amrex::Real z_p = Cu59_pf;

    amrex::Real dz_r_dT = Ni58_pf * dp_pf_dT + p_pf * dNi58_pf_dT;
    amrex::Real dz_p_dT = dCu59_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Mn55_to_n_Fe55_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Mn55 + p --> n + Fe55

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  4.300924891891836 + -11.761875387476456 * tfactors.T9i + 7.03052 * tfactors.T913
                         + 0.260509 * tfactors.T9 + -0.0813397 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  11.761875387476456 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 7.03052 * tfactors.T923i
                                  + 0.260509 + (5.0/3.0) * -0.0813397 * tfactors.T923;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Fe55_pf, dFe55_pf_dT;
    // interpolating Fe55 partition function
    get_partition_function(Fe55, tfactors, Fe55_pf, dFe55_pf_dT);

    amrex::Real Mn55_pf, dMn55_pf_dT;
    // interpolating Mn55 partition function
    get_partition_function(Mn55, tfactors, Mn55_pf, dMn55_pf_dT);

    amrex::Real z_r = n_pf * Fe55_pf;
    amrex::Real z_p = p_pf * Mn55_pf;

    amrex::Real dz_r_dT = Fe55_pf * dn_pf_dT + n_pf * dFe55_pf_dT;
    amrex::Real dz_p_dT = Mn55_pf * dp_pf_dT + p_pf * dMn55_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Fe53_to_n_Ni56_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe53 + He4 --> n + Ni56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  -11.227076237025186 + -31.193524937158223 * tfactors.T9i + 6.43259 * tfactors.T913
                         + 2.29591 * tfactors.T9 + -0.254724 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  31.193524937158223 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 6.43259 * tfactors.T923i
                                  + 2.29591 + (5.0/3.0) * -0.254724 * tfactors.T923;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real Ni56_pf, dNi56_pf_dT;
    // interpolating Ni56 partition function
    get_partition_function(Ni56, tfactors, Ni56_pf, dNi56_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Fe53_pf, dFe53_pf_dT;
    // interpolating Fe53 partition function
    get_partition_function(Fe53, tfactors, Fe53_pf, dFe53_pf_dT);

    amrex::Real z_r = n_pf * Ni56_pf;
    amrex::Real z_p = He4_pf * Fe53_pf;

    amrex::Real dz_r_dT = Ni56_pf * dn_pf_dT + n_pf * dNi56_pf_dT;
    amrex::Real dz_p_dT = Fe53_pf * dHe4_pf_dT + He4_pf * dFe53_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Fe54_to_He4_Mn51_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe54 + p --> He4 + Mn51

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  65.37001620958803 + -36.52011480022066 * tfactors.T9i + -89.274 * tfactors.T913i + -0.862452 * tfactors.T913
                         + -0.635672 * tfactors.T9 + 0.0196464 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  36.52011480022066 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -89.274 * tfactors.T943i + (1.0/3.0) * -0.862452 * tfactors.T923i
                                  + -0.635672 + (5.0/3.0) * 0.0196464 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Mn51_pf, dMn51_pf_dT;
    // interpolating Mn51 partition function
    get_partition_function(Mn51, tfactors, Mn51_pf, dMn51_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Fe54_pf, dFe54_pf_dT;
    // interpolating Fe54 partition function
    get_partition_function(Fe54, tfactors, Fe54_pf, dFe54_pf_dT);

    amrex::Real z_r = He4_pf * Mn51_pf;
    amrex::Real z_p = p_pf * Fe54_pf;

    amrex::Real dz_r_dT = Mn51_pf * dHe4_pf_dT + He4_pf * dMn51_pf_dT;
    amrex::Real dz_p_dT = Fe54_pf * dp_pf_dT + p_pf * dFe54_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Fe54_to_n_Ni57_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe54 + He4 --> n + Ni57

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  14.613800831905415 + -67.50220541720067 * tfactors.T9i + 0.48562 * tfactors.T913
                         + -0.0249302 * tfactors.T9 + 0.0509605 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  67.50220541720067 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 0.48562 * tfactors.T923i
                                  + -0.0249302 + (5.0/3.0) * 0.0509605 * tfactors.T923;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real Ni57_pf, dNi57_pf_dT;
    // interpolating Ni57 partition function
    get_partition_function(Ni57, tfactors, Ni57_pf, dNi57_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Fe54_pf, dFe54_pf_dT;
    // interpolating Fe54 partition function
    get_partition_function(Fe54, tfactors, Fe54_pf, dFe54_pf_dT);

    amrex::Real z_r = n_pf * Ni57_pf;
    amrex::Real z_p = He4_pf * Fe54_pf;

    amrex::Real dz_r_dT = Ni57_pf * dn_pf_dT + n_pf * dNi57_pf_dT;
    amrex::Real dz_p_dT = Fe54_pf * dHe4_pf_dT + He4_pf * dFe54_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Fe54_to_p_Co57_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe54 + He4 --> p + Co57

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  -1.4219519875346407 + -20.56680351235092 * tfactors.T9i + -38.1133 * tfactors.T913i + 29.3541 * tfactors.T913
                         + -4.75966 * tfactors.T9 + 0.40418 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  20.56680351235092 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -38.1133 * tfactors.T943i + (1.0/3.0) * 29.3541 * tfactors.T923i
                                  + -4.75966 + (5.0/3.0) * 0.40418 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Co57_pf, dCo57_pf_dT;
    // interpolating Co57 partition function
    get_partition_function(Co57, tfactors, Co57_pf, dCo57_pf_dT);

    amrex::Real Fe54_pf, dFe54_pf_dT;
    // interpolating Fe54 partition function
    get_partition_function(Fe54, tfactors, Fe54_pf, dFe54_pf_dT);

    amrex::Real z_r = p_pf * Co57_pf;
    amrex::Real z_p = He4_pf * Fe54_pf;

    amrex::Real dz_r_dT = Co57_pf * dp_pf_dT + p_pf * dCo57_pf_dT;
    amrex::Real dz_p_dT = Fe54_pf * dHe4_pf_dT + He4_pf * dFe54_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Fe55_to_n_Co55_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe55 + p --> n + Co55

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  21.432947180559943 + -49.13538645037029 * tfactors.T9i + -1.62382 * tfactors.T913
                         + 0.58115 * tfactors.T9 + -0.0537057 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  49.13538645037029 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.62382 * tfactors.T923i
                                  + 0.58115 + (5.0/3.0) * -0.0537057 * tfactors.T923;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Co55_pf, dCo55_pf_dT;
    // interpolating Co55 partition function
    get_partition_function(Co55, tfactors, Co55_pf, dCo55_pf_dT);

    amrex::Real Fe55_pf, dFe55_pf_dT;
    // interpolating Fe55 partition function
    get_partition_function(Fe55, tfactors, Fe55_pf, dFe55_pf_dT);

    amrex::Real z_r = n_pf * Co55_pf;
    amrex::Real z_p = p_pf * Fe55_pf;

    amrex::Real dz_r_dT = Co55_pf * dn_pf_dT + n_pf * dCo55_pf_dT;
    amrex::Real dz_p_dT = Fe55_pf * dp_pf_dT + p_pf * dFe55_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Fe55_to_n_Ni58_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe55 + He4 --> n + Ni58

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  -4.731933984268859 + -33.63082187826845 * tfactors.T9i + 3.44996 * tfactors.T913
                         + 2.98226 * tfactors.T9 + -0.387699 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  33.63082187826845 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 3.44996 * tfactors.T923i
                                  + 2.98226 + (5.0/3.0) * -0.387699 * tfactors.T923;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Ni58_pf, dNi58_pf_dT;
    // interpolating Ni58 partition function
    get_partition_function(Ni58, tfactors, Ni58_pf, dNi58_pf_dT);

    amrex::Real Fe55_pf, dFe55_pf_dT;
    // interpolating Fe55 partition function
    get_partition_function(Fe55, tfactors, Fe55_pf, dFe55_pf_dT);

    amrex::Real z_r = n_pf * Ni58_pf;
    amrex::Real z_p = He4_pf * Fe55_pf;

    amrex::Real dz_r_dT = Ni58_pf * dn_pf_dT + n_pf * dNi58_pf_dT;
    amrex::Real dz_p_dT = Fe55_pf * dHe4_pf_dT + He4_pf * dFe55_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Fe56_to_n_Co56_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe56 + p --> n + Co56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  22.65112457733622 + -62.065024496438 * tfactors.T9i + -1.13331 * tfactors.T913
                         + 0.347185 * tfactors.T9 + -0.0328879 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  62.065024496438 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.13331 * tfactors.T923i
                                  + 0.347185 + (5.0/3.0) * -0.0328879 * tfactors.T923;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Fe56_pf, dFe56_pf_dT;
    // interpolating Fe56 partition function
    get_partition_function(Fe56, tfactors, Fe56_pf, dFe56_pf_dT);

    amrex::Real Co56_pf, dCo56_pf_dT;
    // interpolating Co56 partition function
    get_partition_function(Co56, tfactors, Co56_pf, dCo56_pf_dT);

    amrex::Real z_r = n_pf * Co56_pf;
    amrex::Real z_p = p_pf * Fe56_pf;

    amrex::Real dz_r_dT = Co56_pf * dn_pf_dT + n_pf * dCo56_pf_dT;
    amrex::Real dz_p_dT = Fe56_pf * dp_pf_dT + p_pf * dFe56_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Co56_to_n_Ni56_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co56 + p --> n + Ni56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  13.37207542266378 + -33.86221596961605 * tfactors.T9i + 1.76846 * tfactors.T913
                         + 0.197992 * tfactors.T9 + -0.017494 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  33.86221596961605 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 1.76846 * tfactors.T923i
                                  + 0.197992 + (5.0/3.0) * -0.017494 * tfactors.T923;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real Ni56_pf, dNi56_pf_dT;
    // interpolating Ni56 partition function
    get_partition_function(Ni56, tfactors, Ni56_pf, dNi56_pf_dT);

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Co56_pf, dCo56_pf_dT;
    // interpolating Co56 partition function
    get_partition_function(Co56, tfactors, Co56_pf, dCo56_pf_dT);

    amrex::Real z_r = n_pf * Ni56_pf;
    amrex::Real z_p = p_pf * Co56_pf;

    amrex::Real dz_r_dT = Ni56_pf * dn_pf_dT + n_pf * dNi56_pf_dT;
    amrex::Real dz_p_dT = Co56_pf * dp_pf_dT + p_pf * dCo56_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Co56_to_He4_Fe53_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co56 + p --> He4 + Fe53

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  56.04732165968897 + -2.6687374505303194 * tfactors.T9i + -91.7231 * tfactors.T913i + 10.052 * tfactors.T913
                         + -3.86332 * tfactors.T9 + 0.391105 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  2.6687374505303194 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -91.7231 * tfactors.T943i + (1.0/3.0) * 10.052 * tfactors.T923i
                                  + -3.86332 + (5.0/3.0) * 0.391105 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Fe53_pf, dFe53_pf_dT;
    // interpolating Fe53 partition function
    get_partition_function(Fe53, tfactors, Fe53_pf, dFe53_pf_dT);

    amrex::Real Co56_pf, dCo56_pf_dT;
    // interpolating Co56 partition function
    get_partition_function(Co56, tfactors, Co56_pf, dCo56_pf_dT);

    amrex::Real z_r = He4_pf * Fe53_pf;
    amrex::Real z_p = p_pf * Co56_pf;

    amrex::Real dz_r_dT = Fe53_pf * dHe4_pf_dT + He4_pf * dFe53_pf_dT;
    amrex::Real dz_p_dT = Co56_pf * dp_pf_dT + p_pf * dCo56_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Co56_to_n_Cu59_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co56 + He4 --> n + Cu59

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  8.891854052419635 + -61.835022947265 * tfactors.T9i + -1.00882 * tfactors.T913
                         + 1.88993 * tfactors.T9 + -0.202319 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  61.835022947265 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.00882 * tfactors.T923i
                                  + 1.88993 + (5.0/3.0) * -0.202319 * tfactors.T923;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real Cu59_pf, dCu59_pf_dT;
    // interpolating Cu59 partition function
    get_partition_function(Cu59, tfactors, Cu59_pf, dCu59_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Co56_pf, dCo56_pf_dT;
    // interpolating Co56 partition function
    get_partition_function(Co56, tfactors, Co56_pf, dCo56_pf_dT);

    amrex::Real z_r = n_pf * Cu59_pf;
    amrex::Real z_p = He4_pf * Co56_pf;

    amrex::Real dz_r_dT = Cu59_pf * dn_pf_dT + n_pf * dCu59_pf_dT;
    amrex::Real dz_p_dT = Co56_pf * dHe4_pf_dT + He4_pf * dCo56_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Co57_to_n_Ni57_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co57 + p --> n + Ni57

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  19.300752819440056 + -46.935401904849755 * tfactors.T9i + -1.13177 * tfactors.T913
                         + 0.56875 * tfactors.T9 + -0.0579913 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  46.935401904849755 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.13177 * tfactors.T923i
                                  + 0.56875 + (5.0/3.0) * -0.0579913 * tfactors.T923;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real Ni57_pf, dNi57_pf_dT;
    // interpolating Ni57 partition function
    get_partition_function(Ni57, tfactors, Ni57_pf, dNi57_pf_dT);

    amrex::Real Co57_pf, dCo57_pf_dT;
    // interpolating Co57 partition function
    get_partition_function(Co57, tfactors, Co57_pf, dCo57_pf_dT);

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real z_r = n_pf * Ni57_pf;
    amrex::Real z_p = p_pf * Co57_pf;

    amrex::Real dz_r_dT = Ni57_pf * dn_pf_dT + n_pf * dNi57_pf_dT;
    amrex::Real dz_p_dT = Co57_pf * dp_pf_dT + p_pf * dCo57_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Ni57_to_n_Zn60_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni57 + He4 --> n + Zn60

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  9.728451219341542 + -87.50502934883477 * tfactors.T9i + 5.43552 * tfactors.T913
                         + -1.25169 * tfactors.T9 + 0.148383 * tfactors.T953;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  87.50502934883477 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 5.43552 * tfactors.T923i
                                  + -1.25169 + (5.0/3.0) * 0.148383 * tfactors.T923;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real Ni57_pf, dNi57_pf_dT;
    // interpolating Ni57 partition function
    get_partition_function(Ni57, tfactors, Ni57_pf, dNi57_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Zn60_pf, dZn60_pf_dT;
    // interpolating Zn60 partition function
    get_partition_function(Zn60, tfactors, Zn60_pf, dZn60_pf_dT);

    amrex::Real z_r = n_pf * Zn60_pf;
    amrex::Real z_p = He4_pf * Ni57_pf;

    amrex::Real dz_r_dT = Zn60_pf * dn_pf_dT + n_pf * dZn60_pf_dT;
    amrex::Real dz_p_dT = Ni57_pf * dHe4_pf_dT + He4_pf * dNi57_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Ni58_to_He4_Co55_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni58 + p --> He4 + Co55

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  63.614171164828804 + -15.504564572101842 * tfactors.T9i + -94.1404 * tfactors.T913i + 3.39179 * tfactors.T913
                         + -1.71062 * tfactors.T9 + 0.133003 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  15.504564572101842 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -94.1404 * tfactors.T943i + (1.0/3.0) * 3.39179 * tfactors.T923i
                                  + -1.71062 + (5.0/3.0) * 0.133003 * tfactors.T923 + -0.666667 * tfactors.T9i;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 / 1.0e9;
    }


    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Co55_pf, dCo55_pf_dT;
    // interpolating Co55 partition function
    get_partition_function(Co55, tfactors, Co55_pf, dCo55_pf_dT);

    amrex::Real Ni58_pf, dNi58_pf_dT;
    // interpolating Ni58 partition function
    get_partition_function(Ni58, tfactors, Ni58_pf, dNi58_pf_dT);

    amrex::Real z_r = He4_pf * Co55_pf;
    amrex::Real z_p = p_pf * Ni58_pf;

    amrex::Real dz_r_dT = Co55_pf * dHe4_pf_dT + He4_pf * dCo55_pf_dT;
    amrex::Real dz_p_dT = Ni58_pf * dp_pf_dT + p_pf * dNi58_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_S32_to_p_Cl35_derived_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // S32 + He4 --> p + Cl35

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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


    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Cl35_pf, dCl35_pf_dT;
    // interpolating Cl35 partition function
    get_partition_function(Cl35, tfactors, Cl35_pf, dCl35_pf_dT);

    amrex::Real S32_pf, dS32_pf_dT;
    // interpolating S32 partition function
    get_partition_function(S32, tfactors, S32_pf, dS32_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real z_r = p_pf * Cl35_pf;
    amrex::Real z_p = He4_pf * S32_pf;

    amrex::Real dz_r_dT = Cl35_pf * dp_pf_dT + p_pf * dCl35_pf_dT;
    amrex::Real dz_p_dT = S32_pf * dHe4_pf_dT + He4_pf * dS32_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ar36_to_He4_S32_derived_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ar36 --> He4 + S32

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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


    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real S32_pf, dS32_pf_dT;
    // interpolating S32 partition function
    get_partition_function(S32, tfactors, S32_pf, dS32_pf_dT);

    amrex::Real Ar36_pf, dAr36_pf_dT;
    // interpolating Ar36 partition function
    get_partition_function(Ar36, tfactors, Ar36_pf, dAr36_pf_dT);

    amrex::Real z_r = He4_pf * S32_pf;
    amrex::Real z_p = Ar36_pf;

    amrex::Real dz_r_dT = S32_pf * dHe4_pf_dT + He4_pf * dS32_pf_dT;
    amrex::Real dz_p_dT = dAr36_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ar36_to_p_Cl35_derived_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ar36 --> p + Cl35

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10c
    ln_set_rate =  16.01685487088493 + -100.72924752614605 * tfactors.T9i;

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
    ln_set_rate =  -17.47510512911507 + -99.28393852614604 * tfactors.T9i;

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
    ln_set_rate =  60.73659487088493 + -98.71928752614605 * tfactors.T9i + -27.8971 * tfactors.T913i + -16.2304 * tfactors.T913
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
    ln_set_rate =  17.202804870884933 + -102.37020752614605 * tfactors.T9i + 18.0179 * tfactors.T913
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


    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Cl35_pf, dCl35_pf_dT;
    // interpolating Cl35 partition function
    get_partition_function(Cl35, tfactors, Cl35_pf, dCl35_pf_dT);

    amrex::Real Ar36_pf, dAr36_pf_dT;
    // interpolating Ar36 partition function
    get_partition_function(Ar36, tfactors, Ar36_pf, dAr36_pf_dT);

    amrex::Real z_r = p_pf * Cl35_pf;
    amrex::Real z_p = Ar36_pf;

    amrex::Real dz_r_dT = Cl35_pf * dp_pf_dT + p_pf * dCl35_pf_dT;
    amrex::Real dz_p_dT = dAr36_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Ar36_to_p_K39_derived_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ar36 + He4 --> p + K39

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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


    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real K39_pf, dK39_pf_dT;
    // interpolating K39 partition function
    get_partition_function(K39, tfactors, K39_pf, dK39_pf_dT);

    amrex::Real Ar36_pf, dAr36_pf_dT;
    // interpolating Ar36 partition function
    get_partition_function(Ar36, tfactors, Ar36_pf, dAr36_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real z_r = p_pf * K39_pf;
    amrex::Real z_p = He4_pf * Ar36_pf;

    amrex::Real dz_r_dT = K39_pf * dp_pf_dT + p_pf * dK39_pf_dT;
    amrex::Real dz_p_dT = Ar36_pf * dHe4_pf_dT + He4_pf * dAr36_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ca40_to_He4_Ar36_derived_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ca40 --> He4 + Ar36

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  77.28261041284823 + -81.69174599574632 * tfactors.T9i + -71.0046 * tfactors.T913i + 4.0656 * tfactors.T913
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


    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Ca40_pf, dCa40_pf_dT;
    // interpolating Ca40 partition function
    get_partition_function(Ca40, tfactors, Ca40_pf, dCa40_pf_dT);

    amrex::Real Ar36_pf, dAr36_pf_dT;
    // interpolating Ar36 partition function
    get_partition_function(Ar36, tfactors, Ar36_pf, dAr36_pf_dT);

    amrex::Real z_r = He4_pf * Ar36_pf;
    amrex::Real z_p = Ca40_pf;

    amrex::Real dz_r_dT = Ar36_pf * dHe4_pf_dT + He4_pf * dAr36_pf_dT;
    amrex::Real dz_p_dT = dCa40_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ca40_to_p_K39_derived_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ca40 --> p + K39

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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


    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Ca40_pf, dCa40_pf_dT;
    // interpolating Ca40 partition function
    get_partition_function(Ca40, tfactors, Ca40_pf, dCa40_pf_dT);

    amrex::Real K39_pf, dK39_pf_dT;
    // interpolating K39 partition function
    get_partition_function(K39, tfactors, K39_pf, dK39_pf_dT);

    amrex::Real z_r = p_pf * K39_pf;
    amrex::Real z_p = Ca40_pf;

    amrex::Real dz_r_dT = K39_pf * dp_pf_dT + p_pf * dK39_pf_dT;
    amrex::Real dz_p_dT = dCa40_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Ca40_to_p_Sc43_derived_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ca40 + He4 --> p + Sc43

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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


    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Ca40_pf, dCa40_pf_dT;
    // interpolating Ca40 partition function
    get_partition_function(Ca40, tfactors, Ca40_pf, dCa40_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Sc43_pf, dSc43_pf_dT;
    // interpolating Sc43 partition function
    get_partition_function(Sc43, tfactors, Sc43_pf, dSc43_pf_dT);

    amrex::Real z_r = p_pf * Sc43_pf;
    amrex::Real z_p = He4_pf * Ca40_pf;

    amrex::Real dz_r_dT = Sc43_pf * dp_pf_dT + p_pf * dSc43_pf_dT;
    amrex::Real dz_p_dT = Ca40_pf * dHe4_pf_dT + He4_pf * dCa40_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ti44_to_He4_Ca40_derived_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ti44 --> He4 + Ca40

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // chw0 
    ln_set_rate =  78.69908591662849 + -59.4975248620018 * tfactors.T9i + -76.4273 * tfactors.T913i + 3.87451 * tfactors.T913
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


    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Ca40_pf, dCa40_pf_dT;
    // interpolating Ca40 partition function
    get_partition_function(Ca40, tfactors, Ca40_pf, dCa40_pf_dT);

    amrex::Real Ti44_pf, dTi44_pf_dT;
    // interpolating Ti44 partition function
    get_partition_function(Ti44, tfactors, Ti44_pf, dTi44_pf_dT);

    amrex::Real z_r = He4_pf * Ca40_pf;
    amrex::Real z_p = Ti44_pf;

    amrex::Real dz_r_dT = Ca40_pf * dHe4_pf_dT + He4_pf * dCa40_pf_dT;
    amrex::Real dz_p_dT = dTi44_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ti44_to_p_Sc43_derived_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ti44 --> p + Sc43

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  62.59391408955788 + -100.37327949403844 * tfactors.T9i + -32.1734 * tfactors.T913i + -1.77078 * tfactors.T913
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


    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Ti44_pf, dTi44_pf_dT;
    // interpolating Ti44 partition function
    get_partition_function(Ti44, tfactors, Ti44_pf, dTi44_pf_dT);

    amrex::Real Sc43_pf, dSc43_pf_dT;
    // interpolating Sc43 partition function
    get_partition_function(Sc43, tfactors, Sc43_pf, dSc43_pf_dT);

    amrex::Real z_r = p_pf * Sc43_pf;
    amrex::Real z_p = Ti44_pf;

    amrex::Real dz_r_dT = Sc43_pf * dp_pf_dT + p_pf * dSc43_pf_dT;
    amrex::Real dz_p_dT = dTi44_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Cr48_to_He4_Ti44_derived_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Cr48 --> He4 + Ti44

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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


    amrex::Real Cr48_pf, dCr48_pf_dT;
    // interpolating Cr48 partition function
    get_partition_function(Cr48, tfactors, Cr48_pf, dCr48_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Ti44_pf, dTi44_pf_dT;
    // interpolating Ti44 partition function
    get_partition_function(Ti44, tfactors, Ti44_pf, dTi44_pf_dT);

    amrex::Real z_r = He4_pf * Ti44_pf;
    amrex::Real z_p = Cr48_pf;

    amrex::Real dz_r_dT = Ti44_pf * dHe4_pf_dT + He4_pf * dTi44_pf_dT;
    amrex::Real dz_p_dT = dCr48_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Cr48_to_p_V47_derived_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Cr48 --> p + V47

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // nfisn
    ln_set_rate =  67.74027107253823 + -100.12633621113824 * tfactors.T9i + -34.0548 * tfactors.T913i + -3.41973 * tfactors.T913
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
    ln_set_rate =  48.891971072538226 + -93.82037121113825 * tfactors.T9i + -45.9868 * tfactors.T913i + 13.6822 * tfactors.T913
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
    ln_set_rate =  65.62307107253824 + -94.58145021113825 * tfactors.T9i + -110.655 * tfactors.T913i + 83.0232 * tfactors.T913
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


    amrex::Real Cr48_pf, dCr48_pf_dT;
    // interpolating Cr48 partition function
    get_partition_function(Cr48, tfactors, Cr48_pf, dCr48_pf_dT);

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real V47_pf, dV47_pf_dT;
    // interpolating V47 partition function
    get_partition_function(V47, tfactors, V47_pf, dV47_pf_dT);

    amrex::Real z_r = p_pf * V47_pf;
    amrex::Real z_p = Cr48_pf;

    amrex::Real dz_r_dT = V47_pf * dp_pf_dT + p_pf * dV47_pf_dT;
    amrex::Real dz_p_dT = dCr48_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_V47_to_He4_Ti44_derived_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // V47 + p --> He4 + Ti44

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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


    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real V47_pf, dV47_pf_dT;
    // interpolating V47 partition function
    get_partition_function(V47, tfactors, V47_pf, dV47_pf_dT);

    amrex::Real Ti44_pf, dTi44_pf_dT;
    // interpolating Ti44 partition function
    get_partition_function(Ti44, tfactors, Ti44_pf, dTi44_pf_dT);

    amrex::Real z_r = He4_pf * Ti44_pf;
    amrex::Real z_p = p_pf * V47_pf;

    amrex::Real dz_r_dT = Ti44_pf * dHe4_pf_dT + He4_pf * dTi44_pf_dT;
    amrex::Real dz_p_dT = V47_pf * dp_pf_dT + p_pf * dV47_pf_dT;

    amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

    drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    rate *= z_r/z_p;

}


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_S32_He4_to_Ar36_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_ag = rate_eval.screened_rates(k_He4_S32_to_Ar36_removed);
    amrex::Real r_ap = rate_eval.screened_rates(k_He4_S32_to_p_Cl35_derived_removed);
    amrex::Real r_pg = rate_eval.screened_rates(k_p_Cl35_to_Ar36_removed);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_Cl35_to_He4_S32_removed);
    amrex::Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ag + r_ap * r_pg * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ag = rate_eval.dscreened_rates_dT(k_He4_S32_to_Ar36_removed);
        amrex::Real drdT_ap = rate_eval.dscreened_rates_dT(k_He4_S32_to_p_Cl35_derived_removed);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_Cl35_to_Ar36_removed);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_Cl35_to_He4_S32_removed);
        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ar36_to_S32_He4_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_ga = rate_eval.screened_rates(k_Ar36_to_He4_S32_derived_removed);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_Cl35_to_He4_S32_removed);
    amrex::Real r_gp = rate_eval.screened_rates(k_Ar36_to_p_Cl35_derived_removed);
    amrex::Real r_pg = rate_eval.screened_rates(k_p_Cl35_to_Ar36_removed);
    amrex::Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ga + r_gp * r_pa * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ga = rate_eval.dscreened_rates_dT(k_Ar36_to_He4_S32_derived_removed);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_Cl35_to_He4_S32_removed);
        amrex::Real drdT_gp = rate_eval.dscreened_rates_dT(k_Ar36_to_p_Cl35_derived_removed);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_Cl35_to_Ar36_removed);
        drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ar36_He4_to_Ca40_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_ag = rate_eval.screened_rates(k_He4_Ar36_to_Ca40_removed);
    amrex::Real r_ap = rate_eval.screened_rates(k_He4_Ar36_to_p_K39_derived_removed);
    amrex::Real r_pg = rate_eval.screened_rates(k_p_K39_to_Ca40_removed);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_K39_to_He4_Ar36_removed);
    amrex::Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ag + r_ap * r_pg * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ag = rate_eval.dscreened_rates_dT(k_He4_Ar36_to_Ca40_removed);
        amrex::Real drdT_ap = rate_eval.dscreened_rates_dT(k_He4_Ar36_to_p_K39_derived_removed);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_K39_to_Ca40_removed);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_K39_to_He4_Ar36_removed);
        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ca40_to_Ar36_He4_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_ga = rate_eval.screened_rates(k_Ca40_to_He4_Ar36_derived_removed);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_K39_to_He4_Ar36_removed);
    amrex::Real r_gp = rate_eval.screened_rates(k_Ca40_to_p_K39_derived_removed);
    amrex::Real r_pg = rate_eval.screened_rates(k_p_K39_to_Ca40_removed);
    amrex::Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ga + r_gp * r_pa * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ga = rate_eval.dscreened_rates_dT(k_Ca40_to_He4_Ar36_derived_removed);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_K39_to_He4_Ar36_removed);
        amrex::Real drdT_gp = rate_eval.dscreened_rates_dT(k_Ca40_to_p_K39_derived_removed);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_K39_to_Ca40_removed);
        drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ca40_He4_to_Ti44_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_ag = rate_eval.screened_rates(k_He4_Ca40_to_Ti44_removed);
    amrex::Real r_ap = rate_eval.screened_rates(k_He4_Ca40_to_p_Sc43_derived_removed);
    amrex::Real r_pg = rate_eval.screened_rates(k_p_Sc43_to_Ti44_removed);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_Sc43_to_He4_Ca40_removed);
    amrex::Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ag + r_ap * r_pg * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ag = rate_eval.dscreened_rates_dT(k_He4_Ca40_to_Ti44_removed);
        amrex::Real drdT_ap = rate_eval.dscreened_rates_dT(k_He4_Ca40_to_p_Sc43_derived_removed);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_Sc43_to_Ti44_removed);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_Sc43_to_He4_Ca40_removed);
        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ti44_to_Ca40_He4_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_ga = rate_eval.screened_rates(k_Ti44_to_He4_Ca40_derived_removed);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_Sc43_to_He4_Ca40_removed);
    amrex::Real r_gp = rate_eval.screened_rates(k_Ti44_to_p_Sc43_derived_removed);
    amrex::Real r_pg = rate_eval.screened_rates(k_p_Sc43_to_Ti44_removed);
    amrex::Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ga + r_gp * r_pa * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ga = rate_eval.dscreened_rates_dT(k_Ti44_to_He4_Ca40_derived_removed);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_Sc43_to_He4_Ca40_removed);
        amrex::Real drdT_gp = rate_eval.dscreened_rates_dT(k_Ti44_to_p_Sc43_derived_removed);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_Sc43_to_Ti44_removed);
        drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ti44_He4_to_Cr48_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_ag = rate_eval.screened_rates(k_He4_Ti44_to_Cr48_removed);
    amrex::Real r_ap = rate_eval.screened_rates(k_He4_Ti44_to_p_V47_removed);
    amrex::Real r_pg = rate_eval.screened_rates(k_p_V47_to_Cr48_removed);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_V47_to_He4_Ti44_derived_removed);
    amrex::Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ag + r_ap * r_pg * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ag = rate_eval.dscreened_rates_dT(k_He4_Ti44_to_Cr48_removed);
        amrex::Real drdT_ap = rate_eval.dscreened_rates_dT(k_He4_Ti44_to_p_V47_removed);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_V47_to_Cr48_removed);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_V47_to_He4_Ti44_derived_removed);
        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Cr48_to_Ti44_He4_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_ga = rate_eval.screened_rates(k_Cr48_to_He4_Ti44_derived_removed);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_V47_to_He4_Ti44_derived_removed);
    amrex::Real r_gp = rate_eval.screened_rates(k_Cr48_to_p_V47_derived_removed);
    amrex::Real r_pg = rate_eval.screened_rates(k_p_V47_to_Cr48_removed);
    amrex::Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ga + r_gp * r_pa * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ga = rate_eval.dscreened_rates_dT(k_Cr48_to_He4_Ti44_derived_removed);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_V47_to_He4_Ti44_derived_removed);
        amrex::Real drdT_gp = rate_eval.dscreened_rates_dT(k_Cr48_to_p_V47_derived_removed);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_V47_to_Cr48_removed);
        drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa);
    }
}


template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
fill_reaclib_rates(const tf_t& tfactors, T& rate_eval)
{

    amrex::Real rate;
    amrex::Real drate_dT;

    rate_p_C12_to_N13<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_C12_to_N13) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_C12_to_N13) = drate_dT;

    }
    rate_He4_C12_to_O16<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_C12_to_O16) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_C12_to_O16) = drate_dT;

    }
    rate_He4_N14_to_F18<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_N14_to_F18) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_N14_to_F18) = drate_dT;

    }
    rate_He4_O16_to_Ne20<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_O16_to_Ne20) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_O16_to_Ne20) = drate_dT;

    }
    rate_He4_F18_to_Na22<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_F18_to_Na22) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_F18_to_Na22) = drate_dT;

    }
    rate_He4_Ne20_to_Mg24<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Ne20_to_Mg24) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Ne20_to_Mg24) = drate_dT;

    }
    rate_p_Ne21_to_Na22<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Ne21_to_Na22) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Ne21_to_Na22) = drate_dT;

    }
    rate_p_Na23_to_Mg24<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Na23_to_Mg24) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Na23_to_Mg24) = drate_dT;

    }
    rate_He4_Mg24_to_Si28<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Mg24_to_Si28) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Mg24_to_Si28) = drate_dT;

    }
    rate_p_Al27_to_Si28<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Al27_to_Si28) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Al27_to_Si28) = drate_dT;

    }
    rate_He4_Si28_to_S32<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Si28_to_S32) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Si28_to_S32) = drate_dT;

    }
    rate_p_P31_to_S32<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_P31_to_S32) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_P31_to_S32) = drate_dT;

    }
    rate_He4_Cr48_to_Fe52<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Cr48_to_Fe52) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Cr48_to_Fe52) = drate_dT;

    }
    rate_p_Mn51_to_Fe52<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Mn51_to_Fe52) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Mn51_to_Fe52) = drate_dT;

    }
    rate_He4_Mn51_to_Co55<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Mn51_to_Co55) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Mn51_to_Co55) = drate_dT;

    }
    rate_He4_Fe52_to_Ni56<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Fe52_to_Ni56) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Fe52_to_Ni56) = drate_dT;

    }
    rate_p_Co55_to_Ni56<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Co55_to_Ni56) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Co55_to_Ni56) = drate_dT;

    }
    rate_He4_Co55_to_Cu59<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Co55_to_Cu59) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Co55_to_Cu59) = drate_dT;

    }
    rate_He4_Ni56_to_Zn60<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Ni56_to_Zn60) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Ni56_to_Zn60) = drate_dT;

    }
    rate_p_Cu59_to_Zn60<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Cu59_to_Zn60) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Cu59_to_Zn60) = drate_dT;

    }
    rate_C12_C12_to_p_Na23<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_C12_C12_to_p_Na23) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_C12_C12_to_p_Na23) = drate_dT;

    }
    rate_C12_C12_to_He4_Ne20<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_C12_C12_to_He4_Ne20) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_C12_C12_to_He4_Ne20) = drate_dT;

    }
    rate_He4_N13_to_p_O16<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_N13_to_p_O16) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_N13_to_p_O16) = drate_dT;

    }
    rate_C12_O16_to_p_Al27<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_C12_O16_to_p_Al27) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_C12_O16_to_p_Al27) = drate_dT;

    }
    rate_C12_O16_to_He4_Mg24<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_C12_O16_to_He4_Mg24) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_C12_O16_to_He4_Mg24) = drate_dT;

    }
    rate_O16_O16_to_p_P31<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_O16_O16_to_p_P31) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_O16_O16_to_p_P31) = drate_dT;

    }
    rate_O16_O16_to_He4_Si28<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_O16_O16_to_He4_Si28) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_O16_O16_to_He4_Si28) = drate_dT;

    }
    rate_He4_F18_to_p_Ne21<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_F18_to_p_Ne21) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_F18_to_p_Ne21) = drate_dT;

    }
    rate_p_Na23_to_He4_Ne20<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Na23_to_He4_Ne20) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Na23_to_He4_Ne20) = drate_dT;

    }
    rate_p_Al27_to_He4_Mg24<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Al27_to_He4_Mg24) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Al27_to_He4_Mg24) = drate_dT;

    }
    rate_p_P31_to_He4_Si28<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_P31_to_He4_Si28) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_P31_to_He4_Si28) = drate_dT;

    }
    rate_He4_Cr48_to_p_Mn51<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Cr48_to_p_Mn51) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Cr48_to_p_Mn51) = drate_dT;

    }
    rate_He4_Fe52_to_p_Co55<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Fe52_to_p_Co55) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Fe52_to_p_Co55) = drate_dT;

    }
    rate_p_Cu59_to_He4_Ni56<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Cu59_to_He4_Ni56) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Cu59_to_He4_Ni56) = drate_dT;

    }
    rate_He4_He4_He4_to_C12<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_He4_He4_to_C12) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_He4_He4_to_C12) = drate_dT;

    }
    rate_C12_C12_to_Mg24_modified<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_C12_C12_to_Mg24_modified) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_C12_C12_to_Mg24_modified) = drate_dT;

    }
    rate_O16_O16_to_S32_modified<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_O16_O16_to_S32_modified) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_O16_O16_to_S32_modified) = drate_dT;

    }
    rate_C12_O16_to_Si28_modified<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_C12_O16_to_Si28_modified) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_C12_O16_to_Si28_modified) = drate_dT;

    }
    rate_p_Mn55_to_Fe56<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Mn55_to_Fe56) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Mn55_to_Fe56) = drate_dT;

    }
    rate_n_Fe52_to_Fe53<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Fe52_to_Fe53) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Fe52_to_Fe53) = drate_dT;

    }
    rate_n_Fe53_to_Fe54<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Fe53_to_Fe54) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Fe53_to_Fe54) = drate_dT;

    }
    rate_He4_Fe53_to_Ni57<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Fe53_to_Ni57) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Fe53_to_Ni57) = drate_dT;

    }
    rate_n_Fe54_to_Fe55<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Fe54_to_Fe55) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Fe54_to_Fe55) = drate_dT;

    }
    rate_p_Fe54_to_Co55<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Fe54_to_Co55) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Fe54_to_Co55) = drate_dT;

    }
    rate_He4_Fe54_to_Ni58<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Fe54_to_Ni58) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Fe54_to_Ni58) = drate_dT;

    }
    rate_n_Fe55_to_Fe56<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Fe55_to_Fe56) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Fe55_to_Fe56) = drate_dT;

    }
    rate_p_Fe55_to_Co56<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Fe55_to_Co56) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Fe55_to_Co56) = drate_dT;

    }
    rate_p_Fe56_to_Co57<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Fe56_to_Co57) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Fe56_to_Co57) = drate_dT;

    }
    rate_n_Co55_to_Co56<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Co55_to_Co56) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Co55_to_Co56) = drate_dT;

    }
    rate_n_Co56_to_Co57<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Co56_to_Co57) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Co56_to_Co57) = drate_dT;

    }
    rate_p_Co56_to_Ni57<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Co56_to_Ni57) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Co56_to_Ni57) = drate_dT;

    }
    rate_p_Co57_to_Ni58<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Co57_to_Ni58) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Co57_to_Ni58) = drate_dT;

    }
    rate_n_Ni56_to_Ni57<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Ni56_to_Ni57) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Ni56_to_Ni57) = drate_dT;

    }
    rate_n_Ni57_to_Ni58<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Ni57_to_Ni58) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Ni57_to_Ni58) = drate_dT;

    }
    rate_p_Ni58_to_Cu59<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Ni58_to_Cu59) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Ni58_to_Cu59) = drate_dT;

    }
    rate_He4_Mn51_to_p_Fe54<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Mn51_to_p_Fe54) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Mn51_to_p_Fe54) = drate_dT;

    }
    rate_He4_Fe53_to_p_Co56<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Fe53_to_p_Co56) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Fe53_to_p_Co56) = drate_dT;

    }
    rate_n_Fe55_to_p_Mn55<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Fe55_to_p_Mn55) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Fe55_to_p_Mn55) = drate_dT;

    }
    rate_n_Co55_to_p_Fe55<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Co55_to_p_Fe55) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Co55_to_p_Fe55) = drate_dT;

    }
    rate_He4_Co55_to_p_Ni58<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Co55_to_p_Ni58) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Co55_to_p_Ni58) = drate_dT;

    }
    rate_n_Co56_to_p_Fe56<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Co56_to_p_Fe56) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Co56_to_p_Fe56) = drate_dT;

    }
    rate_p_Co57_to_He4_Fe54<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Co57_to_He4_Fe54) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Co57_to_He4_Fe54) = drate_dT;

    }
    rate_n_Ni56_to_p_Co56<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Ni56_to_p_Co56) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Ni56_to_p_Co56) = drate_dT;

    }
    rate_n_Ni56_to_He4_Fe53<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Ni56_to_He4_Fe53) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Ni56_to_He4_Fe53) = drate_dT;

    }
    rate_n_Ni57_to_p_Co57<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Ni57_to_p_Co57) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Ni57_to_p_Co57) = drate_dT;

    }
    rate_n_Ni57_to_He4_Fe54<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Ni57_to_He4_Fe54) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Ni57_to_He4_Fe54) = drate_dT;

    }
    rate_n_Ni58_to_He4_Fe55<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Ni58_to_He4_Fe55) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Ni58_to_He4_Fe55) = drate_dT;

    }
    rate_n_Cu59_to_He4_Co56<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Cu59_to_He4_Co56) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Cu59_to_He4_Co56) = drate_dT;

    }
    rate_n_Zn60_to_He4_Ni57<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Zn60_to_He4_Ni57) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Zn60_to_He4_Ni57) = drate_dT;

    }
    rate_He4_S32_to_Ar36_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_S32_to_Ar36_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_S32_to_Ar36_removed) = drate_dT;

    }
    rate_p_Cl35_to_Ar36_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Cl35_to_Ar36_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Cl35_to_Ar36_removed) = drate_dT;

    }
    rate_p_Cl35_to_He4_S32_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Cl35_to_He4_S32_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Cl35_to_He4_S32_removed) = drate_dT;

    }
    rate_He4_Ar36_to_Ca40_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Ar36_to_Ca40_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Ar36_to_Ca40_removed) = drate_dT;

    }
    rate_p_K39_to_Ca40_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_K39_to_Ca40_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_K39_to_Ca40_removed) = drate_dT;

    }
    rate_p_K39_to_He4_Ar36_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_K39_to_He4_Ar36_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_K39_to_He4_Ar36_removed) = drate_dT;

    }
    rate_He4_Ca40_to_Ti44_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Ca40_to_Ti44_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Ca40_to_Ti44_removed) = drate_dT;

    }
    rate_p_Sc43_to_Ti44_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Sc43_to_Ti44_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Sc43_to_Ti44_removed) = drate_dT;

    }
    rate_p_Sc43_to_He4_Ca40_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Sc43_to_He4_Ca40_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Sc43_to_He4_Ca40_removed) = drate_dT;

    }
    rate_He4_Ti44_to_Cr48_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Ti44_to_Cr48_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Ti44_to_Cr48_removed) = drate_dT;

    }
    rate_He4_Ti44_to_p_V47_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Ti44_to_p_V47_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Ti44_to_p_V47_removed) = drate_dT;

    }
    rate_p_V47_to_Cr48_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_V47_to_Cr48_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_V47_to_Cr48_removed) = drate_dT;

    }
    rate_N13_to_p_C12_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_N13_to_p_C12_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_N13_to_p_C12_derived) = drate_dT;

    }
    rate_O16_to_He4_C12_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_O16_to_He4_C12_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_O16_to_He4_C12_derived) = drate_dT;

    }
    rate_F18_to_He4_N14_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_F18_to_He4_N14_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_F18_to_He4_N14_derived) = drate_dT;

    }
    rate_Ne20_to_He4_O16_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Ne20_to_He4_O16_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ne20_to_He4_O16_derived) = drate_dT;

    }
    rate_Na22_to_p_Ne21_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Na22_to_p_Ne21_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Na22_to_p_Ne21_derived) = drate_dT;

    }
    rate_Na22_to_He4_F18_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Na22_to_He4_F18_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Na22_to_He4_F18_derived) = drate_dT;

    }
    rate_Mg24_to_p_Na23_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Mg24_to_p_Na23_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Mg24_to_p_Na23_derived) = drate_dT;

    }
    rate_Mg24_to_He4_Ne20_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Mg24_to_He4_Ne20_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Mg24_to_He4_Ne20_derived) = drate_dT;

    }
    rate_Si28_to_p_Al27_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Si28_to_p_Al27_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Si28_to_p_Al27_derived) = drate_dT;

    }
    rate_Si28_to_He4_Mg24_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Si28_to_He4_Mg24_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Si28_to_He4_Mg24_derived) = drate_dT;

    }
    rate_S32_to_p_P31_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_S32_to_p_P31_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_S32_to_p_P31_derived) = drate_dT;

    }
    rate_S32_to_He4_Si28_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_S32_to_He4_Si28_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_S32_to_He4_Si28_derived) = drate_dT;

    }
    rate_Fe52_to_p_Mn51_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Fe52_to_p_Mn51_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe52_to_p_Mn51_derived) = drate_dT;

    }
    rate_Fe52_to_He4_Cr48_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Fe52_to_He4_Cr48_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe52_to_He4_Cr48_derived) = drate_dT;

    }
    rate_Co55_to_He4_Mn51_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Co55_to_He4_Mn51_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Co55_to_He4_Mn51_derived) = drate_dT;

    }
    rate_Ni56_to_p_Co55_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Ni56_to_p_Co55_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni56_to_p_Co55_derived) = drate_dT;

    }
    rate_Ni56_to_He4_Fe52_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Ni56_to_He4_Fe52_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni56_to_He4_Fe52_derived) = drate_dT;

    }
    rate_Cu59_to_He4_Co55_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Cu59_to_He4_Co55_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Cu59_to_He4_Co55_derived) = drate_dT;

    }
    rate_Zn60_to_p_Cu59_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Zn60_to_p_Cu59_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Zn60_to_p_Cu59_derived) = drate_dT;

    }
    rate_Zn60_to_He4_Ni56_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Zn60_to_He4_Ni56_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Zn60_to_He4_Ni56_derived) = drate_dT;

    }
    rate_C12_to_He4_He4_He4_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_C12_to_He4_He4_He4_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_C12_to_He4_He4_He4_derived) = drate_dT;

    }
    rate_p_O16_to_He4_N13_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_O16_to_He4_N13_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_O16_to_He4_N13_derived) = drate_dT;

    }
    rate_He4_Ne20_to_p_Na23_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Ne20_to_p_Na23_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Ne20_to_p_Na23_derived) = drate_dT;

    }
    rate_p_Ne21_to_He4_F18_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Ne21_to_He4_F18_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Ne21_to_He4_F18_derived) = drate_dT;

    }
    rate_He4_Mg24_to_p_Al27_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Mg24_to_p_Al27_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Mg24_to_p_Al27_derived) = drate_dT;

    }
    rate_He4_Si28_to_p_P31_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Si28_to_p_P31_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Si28_to_p_P31_derived) = drate_dT;

    }
    rate_p_Mn51_to_He4_Cr48_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Mn51_to_He4_Cr48_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Mn51_to_He4_Cr48_derived) = drate_dT;

    }
    rate_p_Co55_to_He4_Fe52_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Co55_to_He4_Fe52_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Co55_to_He4_Fe52_derived) = drate_dT;

    }
    rate_He4_Ni56_to_p_Cu59_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Ni56_to_p_Cu59_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Ni56_to_p_Cu59_derived) = drate_dT;

    }
    rate_Fe53_to_n_Fe52_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Fe53_to_n_Fe52_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe53_to_n_Fe52_derived) = drate_dT;

    }
    rate_Fe54_to_n_Fe53_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Fe54_to_n_Fe53_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe54_to_n_Fe53_derived) = drate_dT;

    }
    rate_Fe55_to_n_Fe54_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Fe55_to_n_Fe54_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe55_to_n_Fe54_derived) = drate_dT;

    }
    rate_Fe56_to_n_Fe55_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Fe56_to_n_Fe55_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe56_to_n_Fe55_derived) = drate_dT;

    }
    rate_Fe56_to_p_Mn55_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Fe56_to_p_Mn55_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe56_to_p_Mn55_derived) = drate_dT;

    }
    rate_Co55_to_p_Fe54_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Co55_to_p_Fe54_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Co55_to_p_Fe54_derived) = drate_dT;

    }
    rate_Co56_to_n_Co55_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Co56_to_n_Co55_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Co56_to_n_Co55_derived) = drate_dT;

    }
    rate_Co56_to_p_Fe55_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Co56_to_p_Fe55_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Co56_to_p_Fe55_derived) = drate_dT;

    }
    rate_Co57_to_n_Co56_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Co57_to_n_Co56_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Co57_to_n_Co56_derived) = drate_dT;

    }
    rate_Co57_to_p_Fe56_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Co57_to_p_Fe56_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Co57_to_p_Fe56_derived) = drate_dT;

    }
    rate_Ni57_to_n_Ni56_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Ni57_to_n_Ni56_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni57_to_n_Ni56_derived) = drate_dT;

    }
    rate_Ni57_to_p_Co56_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Ni57_to_p_Co56_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni57_to_p_Co56_derived) = drate_dT;

    }
    rate_Ni57_to_He4_Fe53_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Ni57_to_He4_Fe53_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni57_to_He4_Fe53_derived) = drate_dT;

    }
    rate_Ni58_to_n_Ni57_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Ni58_to_n_Ni57_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni58_to_n_Ni57_derived) = drate_dT;

    }
    rate_Ni58_to_p_Co57_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Ni58_to_p_Co57_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni58_to_p_Co57_derived) = drate_dT;

    }
    rate_Ni58_to_He4_Fe54_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Ni58_to_He4_Fe54_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni58_to_He4_Fe54_derived) = drate_dT;

    }
    rate_Cu59_to_p_Ni58_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Cu59_to_p_Ni58_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Cu59_to_p_Ni58_derived) = drate_dT;

    }
    rate_p_Mn55_to_n_Fe55_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Mn55_to_n_Fe55_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Mn55_to_n_Fe55_derived) = drate_dT;

    }
    rate_He4_Fe53_to_n_Ni56_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Fe53_to_n_Ni56_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Fe53_to_n_Ni56_derived) = drate_dT;

    }
    rate_p_Fe54_to_He4_Mn51_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Fe54_to_He4_Mn51_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Fe54_to_He4_Mn51_derived) = drate_dT;

    }
    rate_He4_Fe54_to_n_Ni57_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Fe54_to_n_Ni57_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Fe54_to_n_Ni57_derived) = drate_dT;

    }
    rate_He4_Fe54_to_p_Co57_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Fe54_to_p_Co57_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Fe54_to_p_Co57_derived) = drate_dT;

    }
    rate_p_Fe55_to_n_Co55_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Fe55_to_n_Co55_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Fe55_to_n_Co55_derived) = drate_dT;

    }
    rate_He4_Fe55_to_n_Ni58_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Fe55_to_n_Ni58_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Fe55_to_n_Ni58_derived) = drate_dT;

    }
    rate_p_Fe56_to_n_Co56_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Fe56_to_n_Co56_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Fe56_to_n_Co56_derived) = drate_dT;

    }
    rate_p_Co56_to_n_Ni56_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Co56_to_n_Ni56_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Co56_to_n_Ni56_derived) = drate_dT;

    }
    rate_p_Co56_to_He4_Fe53_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Co56_to_He4_Fe53_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Co56_to_He4_Fe53_derived) = drate_dT;

    }
    rate_He4_Co56_to_n_Cu59_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Co56_to_n_Cu59_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Co56_to_n_Cu59_derived) = drate_dT;

    }
    rate_p_Co57_to_n_Ni57_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Co57_to_n_Ni57_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Co57_to_n_Ni57_derived) = drate_dT;

    }
    rate_He4_Ni57_to_n_Zn60_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Ni57_to_n_Zn60_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Ni57_to_n_Zn60_derived) = drate_dT;

    }
    rate_p_Ni58_to_He4_Co55_derived<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Ni58_to_He4_Co55_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Ni58_to_He4_Co55_derived) = drate_dT;

    }
    rate_He4_S32_to_p_Cl35_derived_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_S32_to_p_Cl35_derived_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_S32_to_p_Cl35_derived_removed) = drate_dT;

    }
    rate_Ar36_to_He4_S32_derived_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Ar36_to_He4_S32_derived_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ar36_to_He4_S32_derived_removed) = drate_dT;

    }
    rate_Ar36_to_p_Cl35_derived_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Ar36_to_p_Cl35_derived_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ar36_to_p_Cl35_derived_removed) = drate_dT;

    }
    rate_He4_Ar36_to_p_K39_derived_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Ar36_to_p_K39_derived_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Ar36_to_p_K39_derived_removed) = drate_dT;

    }
    rate_Ca40_to_He4_Ar36_derived_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Ca40_to_He4_Ar36_derived_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ca40_to_He4_Ar36_derived_removed) = drate_dT;

    }
    rate_Ca40_to_p_K39_derived_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Ca40_to_p_K39_derived_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ca40_to_p_K39_derived_removed) = drate_dT;

    }
    rate_He4_Ca40_to_p_Sc43_derived_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Ca40_to_p_Sc43_derived_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Ca40_to_p_Sc43_derived_removed) = drate_dT;

    }
    rate_Ti44_to_He4_Ca40_derived_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Ti44_to_He4_Ca40_derived_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ti44_to_He4_Ca40_derived_removed) = drate_dT;

    }
    rate_Ti44_to_p_Sc43_derived_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Ti44_to_p_Sc43_derived_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ti44_to_p_Sc43_derived_removed) = drate_dT;

    }
    rate_Cr48_to_He4_Ti44_derived_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Cr48_to_He4_Ti44_derived_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Cr48_to_He4_Ti44_derived_removed) = drate_dT;

    }
    rate_Cr48_to_p_V47_derived_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_Cr48_to_p_V47_derived_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Cr48_to_p_V47_derived_removed) = drate_dT;

    }
    rate_p_V47_to_He4_Ti44_derived_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_V47_to_He4_Ti44_derived_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_V47_to_He4_Ti44_derived_removed) = drate_dT;

    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
fill_approx_rates([[maybe_unused]] const tf_t& tfactors, [[maybe_unused]] T& rate_eval)
{

    [[maybe_unused]] amrex::Real rate{};
    [[maybe_unused]] amrex::Real drate_dT{};

    rate_S32_He4_to_Ar36_approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_S32_He4_to_Ar36_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_S32_He4_to_Ar36_approx) = drate_dT;

    }
    rate_Ar36_to_S32_He4_approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_Ar36_to_S32_He4_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ar36_to_S32_He4_approx) = drate_dT;

    }
    rate_Ar36_He4_to_Ca40_approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_Ar36_He4_to_Ca40_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ar36_He4_to_Ca40_approx) = drate_dT;

    }
    rate_Ca40_to_Ar36_He4_approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_Ca40_to_Ar36_He4_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ca40_to_Ar36_He4_approx) = drate_dT;

    }
    rate_Ca40_He4_to_Ti44_approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_Ca40_He4_to_Ti44_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ca40_He4_to_Ti44_approx) = drate_dT;

    }
    rate_Ti44_to_Ca40_He4_approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_Ti44_to_Ca40_He4_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ti44_to_Ca40_He4_approx) = drate_dT;

    }
    rate_Ti44_He4_to_Cr48_approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_Ti44_He4_to_Cr48_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ti44_He4_to_Cr48_approx) = drate_dT;

    }
    rate_Cr48_to_Ti44_He4_approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_Cr48_to_Ti44_He4_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Cr48_to_Ti44_He4_approx) = drate_dT;

    }

}

#endif
