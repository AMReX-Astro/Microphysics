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
void rate_p_Na23_to_Mg24_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Na23 + p --> Mg24

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Mg24_to_Si28_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // Mg24 + He4 --> Si28

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Al27_to_Si28_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Si28_to_S32_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_P31_to_S32_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Cr48_to_Fe52_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Mn51_to_Fe52_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Mn51_to_Co55_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Fe52_to_Ni56_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Co55_to_Ni56_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_C12_C12_to_p_Na23_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
void rate_C12_O16_to_p_Al27_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
void rate_O16_O16_to_p_P31_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_O16_O16_to_He4_Si28_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Na23_to_He4_Ne20_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Al27_to_He4_Mg24_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_P31_to_He4_Si28_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Cr48_to_p_Mn51_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Fe52_to_p_Co55_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_C12_C12_to_n_Mg23_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // C12 + C12 --> n + Mg23

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_O16_O16_to_n_S31_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // O16 + O16 --> n + S31

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_C12_O16_to_n_Si27_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // O16 + C12 --> n + Si27

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Fe54_to_Co55_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Fe54_to_Ni58_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Fe56_to_Co57_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_n_Co55_to_Co56_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_n_Co56_to_Co57_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Co57_to_Ni58_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Mn51_to_p_Fe54_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Co55_to_p_Ni58_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_n_Co56_to_p_Fe56_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Co57_to_He4_Fe54_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_n_Ni56_to_p_Co56_reaclib(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_n_Fe52_to_Fe53_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_n_Fe53_to_Fe54_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_n_Fe54_to_Fe55_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_n_Fe55_to_Fe56_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_n_Ni56_to_Ni57_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_n_Ni57_to_Ni58_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

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
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_C12_C12_to_Mg24_modified(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // C12 + C12 --> Mg24 (calls the underlying rate)

    rate_C12_C12_to_n_Mg23_removed<do_T_derivatives>(tfactors, rate, drate_dT);
}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_O16_O16_to_S32_modified(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // O16 + O16 --> S32 (calls the underlying rate)

    rate_O16_O16_to_n_S31_removed<do_T_derivatives>(tfactors, rate, drate_dT);
}

template <int do_T_derivatives>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_C12_O16_to_Si28_modified(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT) {

    // O16 + C12 --> Si28 (calls the underlying rate)

    rate_C12_O16_to_n_Si27_removed<do_T_derivatives>(tfactors, rate, drate_dT);
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
void rate_Mg24_to_p_Na23_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mg24 --> p + Na23

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0023893015506829445;
    amrex::Real Q_kBT = -11.692695564000061 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_Na23_to_Mg24_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_Na23_to_Mg24_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real Na23_pf, dNa23_pf_dT;
    // interpolating Na23 partition function
    get_partition_function_cached(Na23, tfactors, pf_cache, Na23_pf, dNa23_pf_dT);

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Mg24_pf, dMg24_pf_dT;
    // interpolating Mg24 partition function
    get_partition_function_cached(Mg24, tfactors, pf_cache, Mg24_pf, dMg24_pf_dT);

    amrex::Real z_r = p_pf * Na23_pf;
    amrex::Real z_p = Mg24_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Na23_pf * dp_pf_dT + p_pf * dNa23_pf_dT;
        amrex::Real dz_p_dT = dMg24_pf_dT;

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
void rate_Si28_to_p_Al27_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Si28 --> p + Al27

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0036170500002771144;
    amrex::Real Q_kBT = -11.584908164001932 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_Al27_to_Si28_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_Al27_to_Si28_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Al27_pf, dAl27_pf_dT;
    // interpolating Al27 partition function
    get_partition_function_cached(Al27, tfactors, pf_cache, Al27_pf, dAl27_pf_dT);

    amrex::Real Si28_pf, dSi28_pf_dT;
    // interpolating Si28 partition function
    get_partition_function_cached(Si28, tfactors, pf_cache, Si28_pf, dSi28_pf_dT);

    amrex::Real z_r = p_pf * Al27_pf;
    amrex::Real z_p = Si28_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Al27_pf * dp_pf_dT + p_pf * dAl27_pf_dT;
        amrex::Real dz_p_dT = dSi28_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Si28_to_He4_Mg24_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Si28 --> He4 + Mg24

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0019853803931734024;
    amrex::Real Q_kBT = -9.98413496999995 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_Mg24_to_Si28_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_Mg24_to_Si28_reaclib);
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

    amrex::Real Si28_pf, dSi28_pf_dT;
    // interpolating Si28 partition function
    get_partition_function_cached(Si28, tfactors, pf_cache, Si28_pf, dSi28_pf_dT);

    amrex::Real z_r = He4_pf * Mg24_pf;
    amrex::Real z_p = Si28_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Mg24_pf * dHe4_pf_dT + He4_pf * dMg24_pf_dT;
        amrex::Real dz_p_dT = dSi28_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_S32_to_p_P31_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // S32 --> p + P31

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0012137199057359403;
    amrex::Real Q_kBT = -8.863963764002593 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_P31_to_S32_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_P31_to_S32_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real S32_pf, dS32_pf_dT;
    // interpolating S32 partition function
    get_partition_function_cached(S32, tfactors, pf_cache, S32_pf, dS32_pf_dT);

    amrex::Real P31_pf, dP31_pf_dT;
    // interpolating P31 partition function
    get_partition_function_cached(P31, tfactors, pf_cache, P31_pf, dP31_pf_dT);

    amrex::Real z_r = p_pf * P31_pf;
    amrex::Real z_p = S32_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = P31_pf * dp_pf_dT + p_pf * dP31_pf_dT;
        amrex::Real dz_p_dT = dS32_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_S32_to_He4_Si28_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // S32 --> He4 + Si28

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0020469666932336;
    amrex::Real Q_kBT = -6.9476558700022695 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_Si28_to_S32_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_Si28_to_S32_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real S32_pf, dS32_pf_dT;
    // interpolating S32 partition function
    get_partition_function_cached(S32, tfactors, pf_cache, S32_pf, dS32_pf_dT);

    amrex::Real Si28_pf, dSi28_pf_dT;
    // interpolating Si28 partition function
    get_partition_function_cached(Si28, tfactors, pf_cache, Si28_pf, dSi28_pf_dT);

    amrex::Real z_r = He4_pf * Si28_pf;
    amrex::Real z_p = S32_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Si28_pf * dHe4_pf_dT + He4_pf * dSi28_pf_dT;
        amrex::Real dz_p_dT = dS32_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Fe52_to_p_Mn51_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe52 --> p + Mn51

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0037087238653468386;
    amrex::Real Q_kBT = -7.377871064004735 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_Mn51_to_Fe52_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_Mn51_to_Fe52_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real Mn51_pf, dMn51_pf_dT;
    // interpolating Mn51 partition function
    get_partition_function_cached(Mn51, tfactors, pf_cache, Mn51_pf, dMn51_pf_dT);

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Fe52_pf, dFe52_pf_dT;
    // interpolating Fe52 partition function
    get_partition_function_cached(Fe52, tfactors, pf_cache, Fe52_pf, dFe52_pf_dT);

    amrex::Real z_r = p_pf * Mn51_pf;
    amrex::Real z_p = Fe52_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Mn51_pf * dp_pf_dT + p_pf * dMn51_pf_dT;
        amrex::Real dz_p_dT = dFe52_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Fe52_to_He4_Cr48_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe52 --> He4 + Cr48

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0022179215444048704;
    amrex::Real Q_kBT = -7.9360158700054875 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_Cr48_to_Fe52_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_Cr48_to_Fe52_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real Cr48_pf, dCr48_pf_dT;
    // interpolating Cr48 partition function
    get_partition_function_cached(Cr48, tfactors, pf_cache, Cr48_pf, dCr48_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Fe52_pf, dFe52_pf_dT;
    // interpolating Fe52 partition function
    get_partition_function_cached(Fe52, tfactors, pf_cache, Fe52_pf, dFe52_pf_dT);

    amrex::Real z_r = He4_pf * Cr48_pf;
    amrex::Real z_p = Fe52_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Cr48_pf * dHe4_pf_dT + He4_pf * dCr48_pf_dT;
        amrex::Real dz_p_dT = dFe52_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Co55_to_He4_Mn51_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co55 --> He4 + Mn51

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0016747914576913715;
    amrex::Real Q_kBT = -8.211715870003866 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_Mn51_to_Co55_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_Mn51_to_Co55_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real Mn51_pf, dMn51_pf_dT;
    // interpolating Mn51 partition function
    get_partition_function_cached(Mn51, tfactors, pf_cache, Mn51_pf, dMn51_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Co55_pf, dCo55_pf_dT;
    // interpolating Co55 partition function
    get_partition_function_cached(Co55, tfactors, pf_cache, Co55_pf, dCo55_pf_dT);

    amrex::Real z_r = He4_pf * Mn51_pf;
    amrex::Real z_p = Co55_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Mn51_pf * dHe4_pf_dT + He4_pf * dMn51_pf_dT;
        amrex::Real dz_p_dT = dCo55_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ni56_to_p_Co55_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni56 --> p + Co55

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.004955314145592046;
    amrex::Real Q_kBT = -7.166571063993911 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_Co55_to_Ni56_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_Co55_to_Ni56_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real Ni56_pf, dNi56_pf_dT;
    // interpolating Ni56 partition function
    get_partition_function_cached(Ni56, tfactors, pf_cache, Ni56_pf, dNi56_pf_dT);

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Co55_pf, dCo55_pf_dT;
    // interpolating Co55 partition function
    get_partition_function_cached(Co55, tfactors, pf_cache, Co55_pf, dCo55_pf_dT);

    amrex::Real z_r = p_pf * Co55_pf;
    amrex::Real z_p = Ni56_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Co55_pf * dp_pf_dT + p_pf * dCo55_pf_dT;
        amrex::Real dz_p_dT = dNi56_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ni56_to_He4_Fe52_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni56 --> He4 + Fe52

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0022377286911973004;
    amrex::Real Q_kBT = -8.000415869992421 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_Fe52_to_Ni56_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_Fe52_to_Ni56_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real Ni56_pf, dNi56_pf_dT;
    // interpolating Ni56 partition function
    get_partition_function_cached(Ni56, tfactors, pf_cache, Ni56_pf, dNi56_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Fe52_pf, dFe52_pf_dT;
    // interpolating Fe52 partition function
    get_partition_function_cached(Fe52, tfactors, pf_cache, Fe52_pf, dFe52_pf_dT);

    amrex::Real z_r = He4_pf * Fe52_pf;
    amrex::Real z_p = Ni56_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Fe52_pf * dHe4_pf_dT + He4_pf * dFe52_pf_dT;
        amrex::Real dz_p_dT = dNi56_pf_dT;

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
void rate_He4_Ne20_to_p_Na23_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ne20 + He4 --> p + Na23

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 1.2552512120460946;
    amrex::Real Q_kBT = -2.3761338939988508 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_Na23_to_He4_Ne20_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_Na23_to_He4_Ne20_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real Na23_pf, dNa23_pf_dT;
    // interpolating Na23 partition function
    get_partition_function_cached(Na23, tfactors, pf_cache, Na23_pf, dNa23_pf_dT);

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
    get_partition_function_cached(Ne20, tfactors, pf_cache, Ne20_pf, dNe20_pf_dT);

    amrex::Real z_r = p_pf * Na23_pf;
    amrex::Real z_p = He4_pf * Ne20_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Na23_pf * dp_pf_dT + p_pf * dNa23_pf_dT;
        amrex::Real dz_p_dT = Ne20_pf * dHe4_pf_dT + He4_pf * dNe20_pf_dT;

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
void rate_p_Na23_to_C12_C12_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Na23 + p --> C12 + C12

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.9612850916175811;
    amrex::Real Q_kBT = -2.2408824359984667 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_C12_C12_to_p_Na23_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_C12_C12_to_p_Na23_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real C12_pf, dC12_pf_dT;
    // setting C12 partition function to 1.0 by default, independent of T
    C12_pf = 1.0_rt;
    dC12_pf_dT = 0.0_rt;

    amrex::Real Na23_pf, dNa23_pf_dT;
    // interpolating Na23 partition function
    get_partition_function_cached(Na23, tfactors, pf_cache, Na23_pf, dNa23_pf_dT);

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real z_r = C12_pf * C12_pf;
    amrex::Real z_p = p_pf * Na23_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = dC12_pf_dT + dC12_pf_dT;
        amrex::Real dz_p_dT = Na23_pf * dp_pf_dT + p_pf * dNa23_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Mg24_to_p_Al27_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mg24 + He4 --> p + Al27

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 1.8218423092693468;
    amrex::Real Q_kBT = -1.6007731940019818 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_Al27_to_He4_Mg24_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_Al27_to_He4_Mg24_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Al27_pf, dAl27_pf_dT;
    // interpolating Al27 partition function
    get_partition_function_cached(Al27, tfactors, pf_cache, Al27_pf, dAl27_pf_dT);

    amrex::Real Mg24_pf, dMg24_pf_dT;
    // interpolating Mg24 partition function
    get_partition_function_cached(Mg24, tfactors, pf_cache, Mg24_pf, dMg24_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real z_r = p_pf * Al27_pf;
    amrex::Real z_p = He4_pf * Mg24_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Al27_pf * dp_pf_dT + p_pf * dAl27_pf_dT;
        amrex::Real dz_p_dT = Mg24_pf * dHe4_pf_dT + He4_pf * dMg24_pf_dT;

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

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Al27_to_C12_O16_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Al27 + p --> C12 + O16

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 1.5511707142760711;
    amrex::Real Q_kBT = -5.170886835997453 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_C12_O16_to_p_Al27_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_C12_O16_to_p_Al27_reaclib);
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

    amrex::Real Al27_pf, dAl27_pf_dT;
    // interpolating Al27 partition function
    get_partition_function_cached(Al27, tfactors, pf_cache, Al27_pf, dAl27_pf_dT);

    amrex::Real O16_pf, dO16_pf_dT;
    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors, pf_cache, O16_pf, dO16_pf_dT);

    amrex::Real z_r = C12_pf * O16_pf;
    amrex::Real z_p = p_pf * Al27_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = O16_pf * dC12_pf_dT + C12_pf * dO16_pf_dT;
        amrex::Real dz_p_dT = Al27_pf * dp_pf_dT + p_pf * dAl27_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Si28_to_p_P31_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Si28 + He4 --> p + P31

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.5929358351300884;
    amrex::Real Q_kBT = -1.9163078940003233 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_P31_to_He4_Si28_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_P31_to_He4_Si28_reaclib);
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

    amrex::Real P31_pf, dP31_pf_dT;
    // interpolating P31 partition function
    get_partition_function_cached(P31, tfactors, pf_cache, P31_pf, dP31_pf_dT);

    amrex::Real Si28_pf, dSi28_pf_dT;
    // interpolating Si28 partition function
    get_partition_function_cached(Si28, tfactors, pf_cache, Si28_pf, dSi28_pf_dT);

    amrex::Real z_r = p_pf * P31_pf;
    amrex::Real z_p = He4_pf * Si28_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = P31_pf * dp_pf_dT + p_pf * dP31_pf_dT;
        amrex::Real dz_p_dT = Si28_pf * dHe4_pf_dT + He4_pf * dSi28_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Si28_to_O16_O16_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Si28 + He4 --> O16 + O16

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 1.7258431621042238;
    amrex::Real Q_kBT = -9.593877029998845 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_O16_O16_to_He4_Si28_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_O16_O16_to_He4_Si28_reaclib);
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

    amrex::Real Si28_pf, dSi28_pf_dT;
    // interpolating Si28 partition function
    get_partition_function_cached(Si28, tfactors, pf_cache, Si28_pf, dSi28_pf_dT);

    amrex::Real z_r = O16_pf * O16_pf;
    amrex::Real z_p = He4_pf * Si28_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = dO16_pf_dT + dO16_pf_dT;
        amrex::Real dz_p_dT = Si28_pf * dHe4_pf_dT + He4_pf * dSi28_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_P31_to_O16_O16_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // P31 + p --> O16 + O16

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 2.910674410032206;
    amrex::Real Q_kBT = -7.677569135998522 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_O16_O16_to_p_P31_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_O16_O16_to_p_P31_reaclib);
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

    amrex::Real P31_pf, dP31_pf_dT;
    // interpolating P31 partition function
    get_partition_function_cached(P31, tfactors, pf_cache, P31_pf, dP31_pf_dT);

    amrex::Real z_r = O16_pf * O16_pf;
    amrex::Real z_p = p_pf * P31_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = dO16_pf_dT + dO16_pf_dT;
        amrex::Real dz_p_dT = P31_pf * dp_pf_dT + p_pf * dP31_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Mn51_to_He4_Cr48_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mn51 + p --> He4 + Cr48

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.5980282234351387;
    amrex::Real Q_kBT = -0.5581448060026141 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_Cr48_to_p_Mn51_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_Cr48_to_p_Mn51_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Cr48_pf, dCr48_pf_dT;
    // interpolating Cr48 partition function
    get_partition_function_cached(Cr48, tfactors, pf_cache, Cr48_pf, dCr48_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Mn51_pf, dMn51_pf_dT;
    // interpolating Mn51 partition function
    get_partition_function_cached(Mn51, tfactors, pf_cache, Mn51_pf, dMn51_pf_dT);

    amrex::Real z_r = He4_pf * Cr48_pf;
    amrex::Real z_p = p_pf * Mn51_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Cr48_pf * dHe4_pf_dT + He4_pf * dCr48_pf_dT;
        amrex::Real dz_p_dT = Mn51_pf * dp_pf_dT + p_pf * dMn51_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Co55_to_He4_Fe52_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co55 + p --> He4 + Fe52

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.45158160016713594;
    amrex::Real Q_kBT = -0.8338448059985103 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_Fe52_to_p_Co55_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_Fe52_to_p_Co55_reaclib);
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

    amrex::Real Co55_pf, dCo55_pf_dT;
    // interpolating Co55 partition function
    get_partition_function_cached(Co55, tfactors, pf_cache, Co55_pf, dCo55_pf_dT);

    amrex::Real Fe52_pf, dFe52_pf_dT;
    // interpolating Fe52 partition function
    get_partition_function_cached(Fe52, tfactors, pf_cache, Fe52_pf, dFe52_pf_dT);

    amrex::Real z_r = He4_pf * Fe52_pf;
    amrex::Real z_p = p_pf * Co55_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Fe52_pf * dHe4_pf_dT + He4_pf * dFe52_pf_dT;
        amrex::Real dz_p_dT = Co55_pf * dp_pf_dT + p_pf * dCo55_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Mg24_to_C12_C12_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mg24 --> C12 + C12

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.002296799960050283;
    amrex::Real Q_kBT = -13.933577999999976 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_C12_C12_to_Mg24_modified);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_C12_C12_to_Mg24_modified);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real C12_pf, dC12_pf_dT;
    // setting C12 partition function to 1.0 by default, independent of T
    C12_pf = 1.0_rt;
    dC12_pf_dT = 0.0_rt;

    amrex::Real Mg24_pf, dMg24_pf_dT;
    // interpolating Mg24 partition function
    get_partition_function_cached(Mg24, tfactors, pf_cache, Mg24_pf, dMg24_pf_dT);

    amrex::Real z_r = C12_pf * C12_pf;
    amrex::Real z_p = Mg24_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = dC12_pf_dT + dC12_pf_dT;
        amrex::Real dz_p_dT = dMg24_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_S32_to_O16_O16_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // S32 --> O16 + O16

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0035327434705723034;
    amrex::Real Q_kBT = -16.541532900001116 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_O16_O16_to_S32_modified);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_O16_O16_to_S32_modified);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real S32_pf, dS32_pf_dT;
    // interpolating S32 partition function
    get_partition_function_cached(S32, tfactors, pf_cache, S32_pf, dS32_pf_dT);

    amrex::Real O16_pf, dO16_pf_dT;
    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors, pf_cache, O16_pf, dO16_pf_dT);

    amrex::Real z_r = O16_pf * O16_pf;
    amrex::Real z_p = S32_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = dO16_pf_dT + dO16_pf_dT;
        amrex::Real dz_p_dT = dS32_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Si28_to_C12_O16_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Si28 --> C12 + O16

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.005610662032502117;
    amrex::Real Q_kBT = -16.755794999999384 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_C12_O16_to_Si28_modified);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_C12_O16_to_Si28_modified);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real C12_pf, dC12_pf_dT;
    // setting C12 partition function to 1.0 by default, independent of T
    C12_pf = 1.0_rt;
    dC12_pf_dT = 0.0_rt;

    amrex::Real O16_pf, dO16_pf_dT;
    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors, pf_cache, O16_pf, dO16_pf_dT);

    amrex::Real Si28_pf, dSi28_pf_dT;
    // interpolating Si28 partition function
    get_partition_function_cached(Si28, tfactors, pf_cache, Si28_pf, dSi28_pf_dT);

    amrex::Real z_r = C12_pf * O16_pf;
    amrex::Real z_p = Si28_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = O16_pf * dC12_pf_dT + C12_pf * dO16_pf_dT;
        amrex::Real dz_p_dT = dSi28_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Co55_to_p_Fe54_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co55 --> p + Fe54

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 7.738014712867733e-05;
    amrex::Real Q_kBT = -5.0643710640004365 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_Fe54_to_Co55_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_Fe54_to_Co55_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Co55_pf, dCo55_pf_dT;
    // interpolating Co55 partition function
    get_partition_function_cached(Co55, tfactors, pf_cache, Co55_pf, dCo55_pf_dT);

    amrex::Real Fe54_pf, dFe54_pf_dT;
    // interpolating Fe54 partition function
    get_partition_function_cached(Fe54, tfactors, pf_cache, Fe54_pf, dFe54_pf_dT);

    amrex::Real z_r = p_pf * Fe54_pf;
    amrex::Real z_p = Co55_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Fe54_pf * dp_pf_dT + p_pf * dFe54_pf_dT;
        amrex::Real dz_p_dT = dCo55_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Co56_to_n_Co55_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co56 --> n + Co55

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0005517947438389569;
    amrex::Real Q_kBT = -10.081818100002103 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_n_Co55_to_Co56_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_n_Co55_to_Co56_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real Co55_pf, dCo55_pf_dT;
    // interpolating Co55 partition function
    get_partition_function_cached(Co55, tfactors, pf_cache, Co55_pf, dCo55_pf_dT);

    amrex::Real Co56_pf, dCo56_pf_dT;
    // interpolating Co56 partition function
    get_partition_function_cached(Co56, tfactors, pf_cache, Co56_pf, dCo56_pf_dT);

    amrex::Real z_r = n_pf * Co55_pf;
    amrex::Real z_p = Co56_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Co55_pf * dn_pf_dT + n_pf * dCo55_pf_dT;
        amrex::Real dz_p_dT = dCo56_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Co57_to_n_Co56_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co57 --> n + Co56

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0006987408286997239;
    amrex::Real Q_kBT = -11.376518100002114 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_n_Co56_to_Co57_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_n_Co56_to_Co57_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real Co56_pf, dCo56_pf_dT;
    // interpolating Co56 partition function
    get_partition_function_cached(Co56, tfactors, pf_cache, Co56_pf, dCo56_pf_dT);

    amrex::Real Co57_pf, dCo57_pf_dT;
    // interpolating Co57 partition function
    get_partition_function_cached(Co57, tfactors, pf_cache, Co57_pf, dCo57_pf_dT);

    amrex::Real z_r = n_pf * Co56_pf;
    amrex::Real z_p = Co57_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Co56_pf * dn_pf_dT + n_pf * dCo56_pf_dT;
        amrex::Real dz_p_dT = dCo57_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Co57_to_p_Fe56_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co57 --> p + Fe56

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 7.745938038729025e-05;
    amrex::Real Q_kBT = -6.027511063998631 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_Fe56_to_Co57_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_Fe56_to_Co57_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Fe56_pf, dFe56_pf_dT;
    // interpolating Fe56 partition function
    get_partition_function_cached(Fe56, tfactors, pf_cache, Fe56_pf, dFe56_pf_dT);

    amrex::Real Co57_pf, dCo57_pf_dT;
    // interpolating Co57 partition function
    get_partition_function_cached(Co57, tfactors, pf_cache, Co57_pf, dCo57_pf_dT);

    amrex::Real z_r = p_pf * Fe56_pf;
    amrex::Real z_p = Co57_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Fe56_pf * dp_pf_dT + p_pf * dFe56_pf_dT;
        amrex::Real dz_p_dT = dCo57_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ni58_to_p_Co57_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni58 --> p + Co57

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.004960201648323489;
    amrex::Real Q_kBT = -8.172171063995753 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_Co57_to_Ni58_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_Co57_to_Ni58_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Ni58_pf, dNi58_pf_dT;
    // interpolating Ni58 partition function
    get_partition_function_cached(Ni58, tfactors, pf_cache, Ni58_pf, dNi58_pf_dT);

    amrex::Real Co57_pf, dCo57_pf_dT;
    // interpolating Co57 partition function
    get_partition_function_cached(Co57, tfactors, pf_cache, Co57_pf, dCo57_pf_dT);

    amrex::Real z_r = p_pf * Co57_pf;
    amrex::Real z_p = Ni58_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Co57_pf * dp_pf_dT + p_pf * dCo57_pf_dT;
        amrex::Real dz_p_dT = dNi58_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ni58_to_He4_Fe54_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni58 --> He4 + Fe54

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.002246417756289001;
    amrex::Real Q_kBT = -6.39921586999708 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_Fe54_to_Ni58_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_Fe54_to_Ni58_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Ni58_pf, dNi58_pf_dT;
    // interpolating Ni58 partition function
    get_partition_function_cached(Ni58, tfactors, pf_cache, Ni58_pf, dNi58_pf_dT);

    amrex::Real Fe54_pf, dFe54_pf_dT;
    // interpolating Fe54 partition function
    get_partition_function_cached(Fe54, tfactors, pf_cache, Fe54_pf, dFe54_pf_dT);

    amrex::Real z_r = He4_pf * Fe54_pf;
    amrex::Real z_p = Ni58_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Fe54_pf * dHe4_pf_dT + He4_pf * dFe54_pf_dT;
        amrex::Real dz_p_dT = dNi58_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Fe54_to_He4_Mn51_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe54 + p --> He4 + Mn51

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 21.64368406933008;
    amrex::Real Q_kBT = -3.1473448060034297 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_Mn51_to_p_Fe54_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_Mn51_to_p_Fe54_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Mn51_pf, dMn51_pf_dT;
    // interpolating Mn51 partition function
    get_partition_function_cached(Mn51, tfactors, pf_cache, Mn51_pf, dMn51_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Fe54_pf, dFe54_pf_dT;
    // interpolating Fe54 partition function
    get_partition_function_cached(Fe54, tfactors, pf_cache, Fe54_pf, dFe54_pf_dT);

    amrex::Real z_r = He4_pf * Mn51_pf;
    amrex::Real z_p = p_pf * Fe54_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Mn51_pf * dHe4_pf_dT + He4_pf * dMn51_pf_dT;
        amrex::Real dz_p_dT = Fe54_pf * dp_pf_dT + p_pf * dFe54_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Fe54_to_p_Co57_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe54 + He4 --> p + Co57

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 2.208049519924361;
    amrex::Real Q_kBT = -1.772955193998673 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_Co57_to_He4_Fe54_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_Co57_to_He4_Fe54_reaclib);
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

    amrex::Real Co57_pf, dCo57_pf_dT;
    // interpolating Co57 partition function
    get_partition_function_cached(Co57, tfactors, pf_cache, Co57_pf, dCo57_pf_dT);

    amrex::Real Fe54_pf, dFe54_pf_dT;
    // interpolating Fe54 partition function
    get_partition_function_cached(Fe54, tfactors, pf_cache, Fe54_pf, dFe54_pf_dT);

    amrex::Real z_r = p_pf * Co57_pf;
    amrex::Real z_p = He4_pf * Fe54_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Co57_pf * dp_pf_dT + p_pf * dCo57_pf_dT;
        amrex::Real dz_p_dT = Fe54_pf * dHe4_pf_dT + He4_pf * dFe54_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Fe56_to_n_Co56_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe56 + p --> n + Co56

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 9.020738678854384;
    amrex::Real Q_kBT = -5.349007036003482 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_n_Co56_to_p_Fe56_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_n_Co56_to_p_Fe56_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
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
    get_partition_function_cached(Fe56, tfactors, pf_cache, Fe56_pf, dFe56_pf_dT);

    amrex::Real Co56_pf, dCo56_pf_dT;
    // interpolating Co56 partition function
    get_partition_function_cached(Co56, tfactors, pf_cache, Co56_pf, dCo56_pf_dT);

    amrex::Real z_r = n_pf * Co56_pf;
    amrex::Real z_p = p_pf * Fe56_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Co56_pf * dn_pf_dT + n_pf * dCo56_pf_dT;
        amrex::Real dz_p_dT = Fe56_pf * dp_pf_dT + p_pf * dFe56_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Co56_to_n_Ni56_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co56 + p --> n + Ni56

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.11135413974304752;
    amrex::Real Q_kBT = -2.915247036008191 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_n_Ni56_to_p_Co56_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_n_Ni56_to_p_Co56_reaclib);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real Ni56_pf, dNi56_pf_dT;
    // interpolating Ni56 partition function
    get_partition_function_cached(Ni56, tfactors, pf_cache, Ni56_pf, dNi56_pf_dT);

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Co56_pf, dCo56_pf_dT;
    // interpolating Co56 partition function
    get_partition_function_cached(Co56, tfactors, pf_cache, Co56_pf, dCo56_pf_dT);

    amrex::Real z_r = n_pf * Ni56_pf;
    amrex::Real z_p = p_pf * Co56_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Ni56_pf * dn_pf_dT + n_pf * dNi56_pf_dT;
        amrex::Real dz_p_dT = Co56_pf * dp_pf_dT + p_pf * dCo56_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_Ni58_to_He4_Co55_derived(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni58 + p --> He4 + Co55

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 29.030931571548685;
    amrex::Real Q_kBT = -1.3348448059966436 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_Co55_to_p_Ni58_reaclib);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_Co55_to_p_Ni58_reaclib);
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

    amrex::Real Co55_pf, dCo55_pf_dT;
    // interpolating Co55 partition function
    get_partition_function_cached(Co55, tfactors, pf_cache, Co55_pf, dCo55_pf_dT);

    amrex::Real Ni58_pf, dNi58_pf_dT;
    // interpolating Ni58 partition function
    get_partition_function_cached(Ni58, tfactors, pf_cache, Ni58_pf, dNi58_pf_dT);

    amrex::Real z_r = He4_pf * Co55_pf;
    amrex::Real z_p = p_pf * Ni58_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Co55_pf * dHe4_pf_dT + He4_pf * dCo55_pf_dT;
        amrex::Real dz_p_dT = Ni58_pf * dp_pf_dT + p_pf * dNi58_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_S32_to_p_Cl35_removed(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // S32 + He4 --> p + Cl35

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 1.164434155428163;
    amrex::Real Q_kBT = -1.8660622939993667 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_Cl35_to_He4_S32_removed);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_Cl35_to_He4_S32_removed);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Cl35_pf, dCl35_pf_dT;
    // interpolating Cl35 partition function
    get_partition_function_cached(Cl35, tfactors, pf_cache, Cl35_pf, dCl35_pf_dT);

    amrex::Real S32_pf, dS32_pf_dT;
    // interpolating S32 partition function
    get_partition_function_cached(S32, tfactors, pf_cache, S32_pf, dS32_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real z_r = p_pf * Cl35_pf;
    amrex::Real z_p = He4_pf * S32_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Cl35_pf * dp_pf_dT + p_pf * dCl35_pf_dT;
        amrex::Real dz_p_dT = S32_pf * dHe4_pf_dT + He4_pf * dS32_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ar36_to_He4_S32_removed(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ar36 --> He4 + S32

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0020957921472510906;
    amrex::Real Q_kBT = -6.640920770003172 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_S32_to_Ar36_removed);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_S32_to_Ar36_removed);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real S32_pf, dS32_pf_dT;
    // interpolating S32 partition function
    get_partition_function_cached(S32, tfactors, pf_cache, S32_pf, dS32_pf_dT);

    amrex::Real Ar36_pf, dAr36_pf_dT;
    // interpolating Ar36 partition function
    get_partition_function_cached(Ar36, tfactors, pf_cache, Ar36_pf, dAr36_pf_dT);

    amrex::Real z_r = He4_pf * S32_pf;
    amrex::Real z_p = Ar36_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = S32_pf * dHe4_pf_dT + He4_pf * dS32_pf_dT;
        amrex::Real dz_p_dT = dAr36_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ar36_to_p_Cl35_removed(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ar36 --> p + Cl35

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0024404119589373;
    amrex::Real Q_kBT = -8.506983063999229 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_Cl35_to_Ar36_removed);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_Cl35_to_Ar36_removed);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Cl35_pf, dCl35_pf_dT;
    // interpolating Cl35 partition function
    get_partition_function_cached(Cl35, tfactors, pf_cache, Cl35_pf, dCl35_pf_dT);

    amrex::Real Ar36_pf, dAr36_pf_dT;
    // interpolating Ar36 partition function
    get_partition_function_cached(Ar36, tfactors, pf_cache, Ar36_pf, dAr36_pf_dT);

    amrex::Real z_r = p_pf * Cl35_pf;
    amrex::Real z_p = Ar36_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Cl35_pf * dp_pf_dT + p_pf * dCl35_pf_dT;
        amrex::Real dz_p_dT = dAr36_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Ar36_to_p_K39_removed(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ar36 + He4 --> p + K39

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 1.1478047817875796;
    amrex::Real Q_kBT = -1.288402194006114 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_K39_to_He4_Ar36_removed);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_K39_to_He4_Ar36_removed);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real K39_pf, dK39_pf_dT;
    // interpolating K39 partition function
    get_partition_function_cached(K39, tfactors, pf_cache, K39_pf, dK39_pf_dT);

    amrex::Real Ar36_pf, dAr36_pf_dT;
    // interpolating Ar36 partition function
    get_partition_function_cached(Ar36, tfactors, pf_cache, Ar36_pf, dAr36_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real z_r = p_pf * K39_pf;
    amrex::Real z_p = He4_pf * Ar36_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = K39_pf * dp_pf_dT + p_pf * dK39_pf_dT;
        amrex::Real dz_p_dT = Ar36_pf * dHe4_pf_dT + He4_pf * dAr36_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ca40_to_He4_Ar36_removed(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ca40 --> He4 + Ar36

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0021352372744340166;
    amrex::Real Q_kBT = -7.039775869997547 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_Ar36_to_Ca40_removed);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_Ar36_to_Ca40_removed);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Ca40_pf, dCa40_pf_dT;
    // interpolating Ca40 partition function
    get_partition_function_cached(Ca40, tfactors, pf_cache, Ca40_pf, dCa40_pf_dT);

    amrex::Real Ar36_pf, dAr36_pf_dT;
    // interpolating Ar36 partition function
    get_partition_function_cached(Ar36, tfactors, pf_cache, Ar36_pf, dAr36_pf_dT);

    amrex::Real z_r = He4_pf * Ar36_pf;
    amrex::Real z_p = Ca40_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Ar36_pf * dHe4_pf_dT + He4_pf * dAr36_pf_dT;
        amrex::Real dz_p_dT = dCa40_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ca40_to_p_K39_removed(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ca40 --> p + K39

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.002450835553846443;
    amrex::Real Q_kBT = -8.32817806400366 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_K39_to_Ca40_removed);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_K39_to_Ca40_removed);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Ca40_pf, dCa40_pf_dT;
    // interpolating Ca40 partition function
    get_partition_function_cached(Ca40, tfactors, pf_cache, Ca40_pf, dCa40_pf_dT);

    amrex::Real K39_pf, dK39_pf_dT;
    // interpolating K39 partition function
    get_partition_function_cached(K39, tfactors, pf_cache, K39_pf, dK39_pf_dT);

    amrex::Real z_r = p_pf * K39_pf;
    amrex::Real z_p = Ca40_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = K39_pf * dp_pf_dT + p_pf * dK39_pf_dT;
        amrex::Real dz_p_dT = dCa40_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_He4_Ca40_to_p_Sc43_removed(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ca40 + He4 --> p + Sc43

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 2.269518190046766;
    amrex::Real Q_kBT = -3.522357194003101 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_Sc43_to_He4_Ca40_removed);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_Sc43_to_He4_Ca40_removed);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Ca40_pf, dCa40_pf_dT;
    // interpolating Ca40 partition function
    get_partition_function_cached(Ca40, tfactors, pf_cache, Ca40_pf, dCa40_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Sc43_pf, dSc43_pf_dT;
    // interpolating Sc43 partition function
    get_partition_function_cached(Sc43, tfactors, pf_cache, Sc43_pf, dSc43_pf_dT);

    amrex::Real z_r = p_pf * Sc43_pf;
    amrex::Real z_p = He4_pf * Ca40_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Sc43_pf * dp_pf_dT + p_pf * dSc43_pf_dT;
        amrex::Real dz_p_dT = Ca40_pf * dHe4_pf_dT + He4_pf * dCa40_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ti44_to_He4_Ca40_removed(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ti44 --> He4 + Ca40

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0021673845668759716;
    amrex::Real Q_kBT = -5.127113869992141 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_Ca40_to_Ti44_removed);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_Ca40_to_Ti44_removed);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Ca40_pf, dCa40_pf_dT;
    // interpolating Ca40 partition function
    get_partition_function_cached(Ca40, tfactors, pf_cache, Ca40_pf, dCa40_pf_dT);

    amrex::Real Ti44_pf, dTi44_pf_dT;
    // interpolating Ti44 partition function
    get_partition_function_cached(Ti44, tfactors, pf_cache, Ti44_pf, dTi44_pf_dT);

    amrex::Real z_r = He4_pf * Ca40_pf;
    amrex::Real z_p = Ti44_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Ca40_pf * dHe4_pf_dT + He4_pf * dCa40_pf_dT;
        amrex::Real dz_p_dT = dTi44_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ti44_to_p_Sc43_removed(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ti44 --> p + Sc43

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.004918918699351648;
    amrex::Real Q_kBT = -8.649471063995243 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_Sc43_to_Ti44_removed);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_Sc43_to_Ti44_removed);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real Ti44_pf, dTi44_pf_dT;
    // interpolating Ti44 partition function
    get_partition_function_cached(Ti44, tfactors, pf_cache, Ti44_pf, dTi44_pf_dT);

    amrex::Real Sc43_pf, dSc43_pf_dT;
    // interpolating Sc43 partition function
    get_partition_function_cached(Sc43, tfactors, pf_cache, Sc43_pf, dSc43_pf_dT);

    amrex::Real z_r = p_pf * Sc43_pf;
    amrex::Real z_p = Ti44_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Sc43_pf * dp_pf_dT + p_pf * dSc43_pf_dT;
        amrex::Real dz_p_dT = dTi44_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Cr48_to_He4_Ti44_removed(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Cr48 --> He4 + Ti44

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0021948664061815595;
    amrex::Real Q_kBT = -7.6973158699968405 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_Ti44_to_Cr48_removed);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_Ti44_to_Cr48_removed);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real Cr48_pf, dCr48_pf_dT;
    // interpolating Cr48 partition function
    get_partition_function_cached(Cr48, tfactors, pf_cache, Cr48_pf, dCr48_pf_dT);

    amrex::Real He4_pf, dHe4_pf_dT;
    // setting He4 partition function to 1.0 by default, independent of T
    He4_pf = 1.0_rt;
    dHe4_pf_dT = 0.0_rt;

    amrex::Real Ti44_pf, dTi44_pf_dT;
    // interpolating Ti44 partition function
    get_partition_function_cached(Ti44, tfactors, pf_cache, Ti44_pf, dTi44_pf_dT);

    amrex::Real z_r = He4_pf * Ti44_pf;
    amrex::Real z_p = Cr48_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Ti44_pf * dHe4_pf_dT + He4_pf * dTi44_pf_dT;
        amrex::Real dz_p_dT = dCr48_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Cr48_to_p_V47_removed(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Cr48 --> p + V47

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.0024665232387831902;
    amrex::Real Q_kBT = -8.102901063992563 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_p_V47_to_Cr48_removed);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_p_V47_to_Cr48_removed);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real Cr48_pf, dCr48_pf_dT;
    // interpolating Cr48 partition function
    get_partition_function_cached(Cr48, tfactors, pf_cache, Cr48_pf, dCr48_pf_dT);

    amrex::Real p_pf, dp_pf_dT;
    // setting p partition function to 1.0 by default, independent of T
    p_pf = 1.0_rt;
    dp_pf_dT = 0.0_rt;

    amrex::Real V47_pf, dV47_pf_dT;
    // interpolating V47 partition function
    get_partition_function_cached(V47, tfactors, pf_cache, V47_pf, dV47_pf_dT);

    amrex::Real z_r = p_pf * V47_pf;
    amrex::Real z_p = Cr48_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = V47_pf * dp_pf_dT + p_pf * dV47_pf_dT;
        amrex::Real dz_p_dT = dCr48_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_p_V47_to_He4_Ti44_removed(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // V47 + p --> He4 + Ti44

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.8898624475414845;
    amrex::Real Q_kBT = 0.4055851939957223 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_He4_Ti44_to_p_V47_removed);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (0.0 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_He4_Ti44_to_p_V47_removed);
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

    amrex::Real V47_pf, dV47_pf_dT;
    // interpolating V47 partition function
    get_partition_function_cached(V47, tfactors, pf_cache, V47_pf, dV47_pf_dT);

    amrex::Real Ti44_pf, dTi44_pf_dT;
    // interpolating Ti44 partition function
    get_partition_function_cached(Ti44, tfactors, pf_cache, Ti44_pf, dTi44_pf_dT);

    amrex::Real z_r = He4_pf * Ti44_pf;
    amrex::Real z_p = p_pf * V47_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Ti44_pf * dHe4_pf_dT + He4_pf * dTi44_pf_dT;
        amrex::Real dz_p_dT = V47_pf * dp_pf_dT + p_pf * dV47_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Fe54_to_n_Fe53_removed(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe54 --> n + Fe53

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.004961965295535373;
    amrex::Real Q_kBT = -13.378418100002976 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_n_Fe53_to_Fe54_removed);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_n_Fe53_to_Fe54_removed);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real Fe53_pf, dFe53_pf_dT;
    // interpolating Fe53 partition function
    get_partition_function_cached(Fe53, tfactors, pf_cache, Fe53_pf, dFe53_pf_dT);

    amrex::Real Fe54_pf, dFe54_pf_dT;
    // interpolating Fe54 partition function
    get_partition_function_cached(Fe54, tfactors, pf_cache, Fe54_pf, dFe54_pf_dT);

    amrex::Real z_r = n_pf * Fe53_pf;
    amrex::Real z_p = Fe54_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Fe53_pf * dn_pf_dT + n_pf * dFe53_pf_dT;
        amrex::Real dz_p_dT = dFe54_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Fe53_to_n_Fe52_removed(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe53 --> n + Fe52

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 7.747917802697315e-05;
    amrex::Real Q_kBT = -10.686718099998421 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_n_Fe52_to_Fe53_removed);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_n_Fe52_to_Fe53_removed);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real Fe52_pf, dFe52_pf_dT;
    // interpolating Fe52 partition function
    get_partition_function_cached(Fe52, tfactors, pf_cache, Fe52_pf, dFe52_pf_dT);

    amrex::Real Fe53_pf, dFe53_pf_dT;
    // interpolating Fe53 partition function
    get_partition_function_cached(Fe53, tfactors, pf_cache, Fe53_pf, dFe53_pf_dT);

    amrex::Real z_r = n_pf * Fe52_pf;
    amrex::Real z_p = Fe53_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Fe52_pf * dn_pf_dT + n_pf * dFe52_pf_dT;
        amrex::Real dz_p_dT = dFe53_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Fe56_to_n_Fe55_removed(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe56 --> n + Fe55

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.002483201755819107;
    amrex::Real Q_kBT = -11.197078100005754 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_n_Fe55_to_Fe56_removed);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_n_Fe55_to_Fe56_removed);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real Fe56_pf, dFe56_pf_dT;
    // interpolating Fe56 partition function
    get_partition_function_cached(Fe56, tfactors, pf_cache, Fe56_pf, dFe56_pf_dT);

    amrex::Real Fe55_pf, dFe55_pf_dT;
    // interpolating Fe55 partition function
    get_partition_function_cached(Fe55, tfactors, pf_cache, Fe55_pf, dFe55_pf_dT);

    amrex::Real z_r = n_pf * Fe55_pf;
    amrex::Real z_p = Fe56_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Fe55_pf * dn_pf_dT + n_pf * dFe55_pf_dT;
        amrex::Real dz_p_dT = dFe56_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Fe55_to_n_Fe54_removed(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe55 --> n + Fe54

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.00015510907420772846;
    amrex::Real Q_kBT = -9.298118100000266 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_n_Fe54_to_Fe55_removed);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_n_Fe54_to_Fe55_removed);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real Fe55_pf, dFe55_pf_dT;
    // interpolating Fe55 partition function
    get_partition_function_cached(Fe55, tfactors, pf_cache, Fe55_pf, dFe55_pf_dT);

    amrex::Real Fe54_pf, dFe54_pf_dT;
    // interpolating Fe54 partition function
    get_partition_function_cached(Fe54, tfactors, pf_cache, Fe54_pf, dFe54_pf_dT);

    amrex::Real z_r = n_pf * Fe54_pf;
    amrex::Real z_p = Fe55_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Fe54_pf * dn_pf_dT + n_pf * dFe54_pf_dT;
        amrex::Real dz_p_dT = dFe55_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ni58_to_n_Ni57_removed(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni58 --> n + Ni57

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.002485653242244135;
    amrex::Real Q_kBT = -12.216218100000214 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_n_Ni57_to_Ni58_removed);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_n_Ni57_to_Ni58_removed);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real Ni57_pf, dNi57_pf_dT;
    // interpolating Ni57 partition function
    get_partition_function_cached(Ni57, tfactors, pf_cache, Ni57_pf, dNi57_pf_dT);

    amrex::Real Ni58_pf, dNi58_pf_dT;
    // interpolating Ni58 partition function
    get_partition_function_cached(Ni58, tfactors, pf_cache, Ni58_pf, dNi58_pf_dT);

    amrex::Real z_r = n_pf * Ni57_pf;
    amrex::Real z_p = Ni58_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Ni57_pf * dn_pf_dT + n_pf * dNi57_pf_dT;
        amrex::Real dz_p_dT = dNi58_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ni57_to_n_Ni56_removed(const T& rate_eval, const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni57 --> n + Ni56

    // Evaluate the equilibrium ratio without partition function
    amrex::Real ratio = 0.00015526775629461862;
    amrex::Real Q_kBT = -10.247718100005844 / (C::k_MeV * tfactors.T9 * 1.0e9_rt);
    ratio *= std::exp(Q_kBT);
    ratio *= std::sqrt(amrex::Math::powi<3>(tfactors.T9 * 1.0e9_rt));

    // Apply the ratio without partition function
    // Note that screening is not yet applied to the inverse rate
    rate = rate_eval.screened_rates(k_n_Ni56_to_Ni57_removed);
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dratio_dT = ratio * tfactors.T9i * 1.0e-9_rt * (1.5 - Q_kBT);
        drate_dT = rate_eval.dscreened_rates_dT(k_n_Ni56_to_Ni57_removed);
        drate_dT = drate_dT * ratio + rate * dratio_dT;
    }
    rate *= ratio;

    // Now apply partition function effects
    amrex::Real n_pf, dn_pf_dT;
    // setting n partition function to 1.0 by default, independent of T
    n_pf = 1.0_rt;
    dn_pf_dT = 0.0_rt;

    amrex::Real Ni56_pf, dNi56_pf_dT;
    // interpolating Ni56 partition function
    get_partition_function_cached(Ni56, tfactors, pf_cache, Ni56_pf, dNi56_pf_dT);

    amrex::Real Ni57_pf, dNi57_pf_dT;
    // interpolating Ni57 partition function
    get_partition_function_cached(Ni57, tfactors, pf_cache, Ni57_pf, dNi57_pf_dT);

    amrex::Real z_r = n_pf * Ni56_pf;
    amrex::Real z_p = Ni57_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dz_r_dT = Ni56_pf * dn_pf_dT + n_pf * dNi56_pf_dT;
        amrex::Real dz_p_dT = dNi57_pf_dT;

        amrex::Real dzterm_dT = (z_p * dz_r_dT - z_r * dz_p_dT) / (z_p * z_p);

        drate_dT = dzterm_dT * rate + drate_dT * (z_r / z_p);
    }
    rate *= z_r / z_p;

}


template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_S32_He4_to_Ar36_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_ag = rate_eval.screened_rates(k_He4_S32_to_Ar36_removed);
    amrex::Real r_ap = rate_eval.screened_rates(k_He4_S32_to_p_Cl35_removed);
    amrex::Real r_pg = rate_eval.screened_rates(k_p_Cl35_to_Ar36_removed);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_Cl35_to_He4_S32_removed);
    amrex::Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ag + r_ap * r_pg * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ag = rate_eval.dscreened_rates_dT(k_He4_S32_to_Ar36_removed);
        amrex::Real drdT_ap = rate_eval.dscreened_rates_dT(k_He4_S32_to_p_Cl35_removed);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_Cl35_to_Ar36_removed);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_Cl35_to_He4_S32_removed);
        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ar36_to_S32_He4_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_ga = rate_eval.screened_rates(k_Ar36_to_He4_S32_removed);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_Cl35_to_He4_S32_removed);
    amrex::Real r_gp = rate_eval.screened_rates(k_Ar36_to_p_Cl35_removed);
    amrex::Real r_pg = rate_eval.screened_rates(k_p_Cl35_to_Ar36_removed);
    amrex::Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ga + r_gp * r_pa * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ga = rate_eval.dscreened_rates_dT(k_Ar36_to_He4_S32_removed);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_Cl35_to_He4_S32_removed);
        amrex::Real drdT_gp = rate_eval.dscreened_rates_dT(k_Ar36_to_p_Cl35_removed);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_Cl35_to_Ar36_removed);
        drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ar36_He4_to_Ca40_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_ag = rate_eval.screened_rates(k_He4_Ar36_to_Ca40_removed);
    amrex::Real r_ap = rate_eval.screened_rates(k_He4_Ar36_to_p_K39_removed);
    amrex::Real r_pg = rate_eval.screened_rates(k_p_K39_to_Ca40_removed);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_K39_to_He4_Ar36_removed);
    amrex::Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ag + r_ap * r_pg * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ag = rate_eval.dscreened_rates_dT(k_He4_Ar36_to_Ca40_removed);
        amrex::Real drdT_ap = rate_eval.dscreened_rates_dT(k_He4_Ar36_to_p_K39_removed);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_K39_to_Ca40_removed);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_K39_to_He4_Ar36_removed);
        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ca40_to_Ar36_He4_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_ga = rate_eval.screened_rates(k_Ca40_to_He4_Ar36_removed);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_K39_to_He4_Ar36_removed);
    amrex::Real r_gp = rate_eval.screened_rates(k_Ca40_to_p_K39_removed);
    amrex::Real r_pg = rate_eval.screened_rates(k_p_K39_to_Ca40_removed);
    amrex::Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ga + r_gp * r_pa * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ga = rate_eval.dscreened_rates_dT(k_Ca40_to_He4_Ar36_removed);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_K39_to_He4_Ar36_removed);
        amrex::Real drdT_gp = rate_eval.dscreened_rates_dT(k_Ca40_to_p_K39_removed);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_K39_to_Ca40_removed);
        drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ca40_He4_to_Ti44_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_ag = rate_eval.screened_rates(k_He4_Ca40_to_Ti44_removed);
    amrex::Real r_ap = rate_eval.screened_rates(k_He4_Ca40_to_p_Sc43_removed);
    amrex::Real r_pg = rate_eval.screened_rates(k_p_Sc43_to_Ti44_removed);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_Sc43_to_He4_Ca40_removed);
    amrex::Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ag + r_ap * r_pg * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ag = rate_eval.dscreened_rates_dT(k_He4_Ca40_to_Ti44_removed);
        amrex::Real drdT_ap = rate_eval.dscreened_rates_dT(k_He4_Ca40_to_p_Sc43_removed);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_Sc43_to_Ti44_removed);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_Sc43_to_He4_Ca40_removed);
        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ti44_to_Ca40_He4_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_ga = rate_eval.screened_rates(k_Ti44_to_He4_Ca40_removed);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_Sc43_to_He4_Ca40_removed);
    amrex::Real r_gp = rate_eval.screened_rates(k_Ti44_to_p_Sc43_removed);
    amrex::Real r_pg = rate_eval.screened_rates(k_p_Sc43_to_Ti44_removed);
    amrex::Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ga + r_gp * r_pa * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ga = rate_eval.dscreened_rates_dT(k_Ti44_to_He4_Ca40_removed);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_Sc43_to_He4_Ca40_removed);
        amrex::Real drdT_gp = rate_eval.dscreened_rates_dT(k_Ti44_to_p_Sc43_removed);
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
    amrex::Real r_pa = rate_eval.screened_rates(k_p_V47_to_He4_Ti44_removed);
    amrex::Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ag + r_ap * r_pg * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ag = rate_eval.dscreened_rates_dT(k_He4_Ti44_to_Cr48_removed);
        amrex::Real drdT_ap = rate_eval.dscreened_rates_dT(k_He4_Ti44_to_p_V47_removed);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_V47_to_Cr48_removed);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_V47_to_He4_Ti44_removed);
        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Cr48_to_Ti44_He4_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_ga = rate_eval.screened_rates(k_Cr48_to_He4_Ti44_removed);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_V47_to_He4_Ti44_removed);
    amrex::Real r_gp = rate_eval.screened_rates(k_Cr48_to_p_V47_removed);
    amrex::Real r_pg = rate_eval.screened_rates(k_p_V47_to_Cr48_removed);
    amrex::Real dd = 1.0_rt / (r_pg + r_pa);
    rate = r_ga + r_gp * r_pa * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ga = rate_eval.dscreened_rates_dT(k_Cr48_to_He4_Ti44_removed);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_V47_to_He4_Ti44_removed);
        amrex::Real drdT_gp = rate_eval.dscreened_rates_dT(k_Cr48_to_p_V47_removed);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_V47_to_Cr48_removed);
        drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Fe52_n_n_to_Fe54_approx(const T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real Yn = Y(N);
    amrex::Real r1_ng = rate_eval.screened_rates(k_n_Fe52_to_Fe53_removed);
    amrex::Real r2_ng = rate_eval.screened_rates(k_n_Fe53_to_Fe54_removed);
    amrex::Real r1_gn = rate_eval.screened_rates(k_Fe53_to_n_Fe52_removed);
    amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
    rate = r1_ng * r2_ng * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_n_Fe52_to_Fe53_removed);
        amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_n_Fe53_to_Fe54_removed);
        amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Fe53_to_n_Fe52_removed);
        drate_dT = dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Fe54_to_Fe52_n_n_approx(const T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real Yn = Y(N);
    amrex::Real r1_gn = rate_eval.screened_rates(k_Fe53_to_n_Fe52_removed);
    amrex::Real r2_gn = rate_eval.screened_rates(k_Fe54_to_n_Fe53_removed);
    amrex::Real r2_ng = rate_eval.screened_rates(k_n_Fe53_to_Fe54_removed);
    amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
    rate = r1_gn * r2_gn * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Fe53_to_n_Fe52_removed);
        amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Fe54_to_n_Fe53_removed);
        amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_n_Fe53_to_Fe54_removed);
        drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Fe54_n_n_to_Fe56_approx(const T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real Yn = Y(N);
    amrex::Real r1_ng = rate_eval.screened_rates(k_n_Fe54_to_Fe55_removed);
    amrex::Real r2_ng = rate_eval.screened_rates(k_n_Fe55_to_Fe56_removed);
    amrex::Real r1_gn = rate_eval.screened_rates(k_Fe55_to_n_Fe54_removed);
    amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
    rate = r1_ng * r2_ng * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_n_Fe54_to_Fe55_removed);
        amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_n_Fe55_to_Fe56_removed);
        amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Fe55_to_n_Fe54_removed);
        drate_dT = dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Fe56_to_Fe54_n_n_approx(const T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real Yn = Y(N);
    amrex::Real r1_gn = rate_eval.screened_rates(k_Fe55_to_n_Fe54_removed);
    amrex::Real r2_gn = rate_eval.screened_rates(k_Fe56_to_n_Fe55_removed);
    amrex::Real r2_ng = rate_eval.screened_rates(k_n_Fe55_to_Fe56_removed);
    amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
    rate = r1_gn * r2_gn * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Fe55_to_n_Fe54_removed);
        amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Fe56_to_n_Fe55_removed);
        amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_n_Fe55_to_Fe56_removed);
        drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ni56_n_n_to_Ni58_approx(const T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real Yn = Y(N);
    amrex::Real r1_ng = rate_eval.screened_rates(k_n_Ni56_to_Ni57_removed);
    amrex::Real r2_ng = rate_eval.screened_rates(k_n_Ni57_to_Ni58_removed);
    amrex::Real r1_gn = rate_eval.screened_rates(k_Ni57_to_n_Ni56_removed);
    amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
    rate = r1_ng * r2_ng * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_n_Ni56_to_Ni57_removed);
        amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_n_Ni57_to_Ni58_removed);
        amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ni57_to_n_Ni56_removed);
        drate_dT = dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ni58_to_Ni56_n_n_approx(const T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real Yn = Y(N);
    amrex::Real r1_gn = rate_eval.screened_rates(k_Ni57_to_n_Ni56_removed);
    amrex::Real r2_gn = rate_eval.screened_rates(k_Ni58_to_n_Ni57_removed);
    amrex::Real r2_ng = rate_eval.screened_rates(k_n_Ni57_to_Ni58_removed);
    amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
    rate = r1_gn * r2_gn * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ni57_to_n_Ni56_removed);
        amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Ni58_to_n_Ni57_removed);
        amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_n_Ni57_to_Ni58_removed);
        drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
    }
}


template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
fill_reaclib_rates(const tf_t& tfactors, T& rate_eval)
{

    amrex::Real rate;
    amrex::Real drate_dT;

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
    rate_He4_O16_to_Ne20_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_O16_to_Ne20_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_O16_to_Ne20_reaclib) = drate_dT;

    }
    rate_He4_Ne20_to_Mg24_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Ne20_to_Mg24_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Ne20_to_Mg24_reaclib) = drate_dT;

    }
    rate_p_Na23_to_Mg24_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Na23_to_Mg24_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Na23_to_Mg24_reaclib) = drate_dT;

    }
    rate_He4_Mg24_to_Si28_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Mg24_to_Si28_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Mg24_to_Si28_reaclib) = drate_dT;

    }
    rate_p_Al27_to_Si28_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Al27_to_Si28_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Al27_to_Si28_reaclib) = drate_dT;

    }
    rate_He4_Si28_to_S32_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Si28_to_S32_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Si28_to_S32_reaclib) = drate_dT;

    }
    rate_p_P31_to_S32_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_P31_to_S32_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_P31_to_S32_reaclib) = drate_dT;

    }
    rate_He4_Cr48_to_Fe52_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Cr48_to_Fe52_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Cr48_to_Fe52_reaclib) = drate_dT;

    }
    rate_p_Mn51_to_Fe52_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Mn51_to_Fe52_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Mn51_to_Fe52_reaclib) = drate_dT;

    }
    rate_He4_Mn51_to_Co55_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Mn51_to_Co55_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Mn51_to_Co55_reaclib) = drate_dT;

    }
    rate_He4_Fe52_to_Ni56_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Fe52_to_Ni56_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Fe52_to_Ni56_reaclib) = drate_dT;

    }
    rate_p_Co55_to_Ni56_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Co55_to_Ni56_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Co55_to_Ni56_reaclib) = drate_dT;

    }
    rate_C12_C12_to_p_Na23_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_C12_C12_to_p_Na23_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_C12_C12_to_p_Na23_reaclib) = drate_dT;

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
    rate_C12_O16_to_p_Al27_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_C12_O16_to_p_Al27_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_C12_O16_to_p_Al27_reaclib) = drate_dT;

    }
    rate_C12_O16_to_He4_Mg24_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_C12_O16_to_He4_Mg24_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_C12_O16_to_He4_Mg24_reaclib) = drate_dT;

    }
    rate_O16_O16_to_p_P31_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_O16_O16_to_p_P31_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_O16_O16_to_p_P31_reaclib) = drate_dT;

    }
    rate_O16_O16_to_He4_Si28_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_O16_O16_to_He4_Si28_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_O16_O16_to_He4_Si28_reaclib) = drate_dT;

    }
    rate_p_Na23_to_He4_Ne20_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Na23_to_He4_Ne20_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Na23_to_He4_Ne20_reaclib) = drate_dT;

    }
    rate_p_Al27_to_He4_Mg24_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Al27_to_He4_Mg24_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Al27_to_He4_Mg24_reaclib) = drate_dT;

    }
    rate_p_P31_to_He4_Si28_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_P31_to_He4_Si28_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_P31_to_He4_Si28_reaclib) = drate_dT;

    }
    rate_He4_Cr48_to_p_Mn51_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Cr48_to_p_Mn51_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Cr48_to_p_Mn51_reaclib) = drate_dT;

    }
    rate_He4_Fe52_to_p_Co55_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Fe52_to_p_Co55_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Fe52_to_p_Co55_reaclib) = drate_dT;

    }
    rate_He4_He4_He4_to_C12_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_He4_He4_to_C12_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_He4_He4_to_C12_reaclib) = drate_dT;

    }
    rate_C12_C12_to_n_Mg23_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_C12_C12_to_n_Mg23_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_C12_C12_to_n_Mg23_removed) = drate_dT;

    }
    rate_O16_O16_to_n_S31_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_O16_O16_to_n_S31_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_O16_O16_to_n_S31_removed) = drate_dT;

    }
    rate_C12_O16_to_n_Si27_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_C12_O16_to_n_Si27_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_C12_O16_to_n_Si27_removed) = drate_dT;

    }
    rate_p_Fe54_to_Co55_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Fe54_to_Co55_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Fe54_to_Co55_reaclib) = drate_dT;

    }
    rate_He4_Fe54_to_Ni58_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Fe54_to_Ni58_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Fe54_to_Ni58_reaclib) = drate_dT;

    }
    rate_p_Fe56_to_Co57_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Fe56_to_Co57_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Fe56_to_Co57_reaclib) = drate_dT;

    }
    rate_n_Co55_to_Co56_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Co55_to_Co56_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Co55_to_Co56_reaclib) = drate_dT;

    }
    rate_n_Co56_to_Co57_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Co56_to_Co57_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Co56_to_Co57_reaclib) = drate_dT;

    }
    rate_p_Co57_to_Ni58_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Co57_to_Ni58_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Co57_to_Ni58_reaclib) = drate_dT;

    }
    rate_He4_Mn51_to_p_Fe54_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Mn51_to_p_Fe54_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Mn51_to_p_Fe54_reaclib) = drate_dT;

    }
    rate_He4_Co55_to_p_Ni58_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_He4_Co55_to_p_Ni58_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Co55_to_p_Ni58_reaclib) = drate_dT;

    }
    rate_n_Co56_to_p_Fe56_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Co56_to_p_Fe56_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Co56_to_p_Fe56_reaclib) = drate_dT;

    }
    rate_p_Co57_to_He4_Fe54_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_p_Co57_to_He4_Fe54_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Co57_to_He4_Fe54_reaclib) = drate_dT;

    }
    rate_n_Ni56_to_p_Co56_reaclib<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Ni56_to_p_Co56_reaclib) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Ni56_to_p_Co56_reaclib) = drate_dT;

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
    rate_n_Fe52_to_Fe53_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Fe52_to_Fe53_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Fe52_to_Fe53_removed) = drate_dT;

    }
    rate_n_Fe53_to_Fe54_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Fe53_to_Fe54_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Fe53_to_Fe54_removed) = drate_dT;

    }
    rate_n_Fe54_to_Fe55_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Fe54_to_Fe55_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Fe54_to_Fe55_removed) = drate_dT;

    }
    rate_n_Fe55_to_Fe56_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Fe55_to_Fe56_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Fe55_to_Fe56_removed) = drate_dT;

    }
    rate_n_Ni56_to_Ni57_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Ni56_to_Ni57_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Ni56_to_Ni57_removed) = drate_dT;

    }
    rate_n_Ni57_to_Ni58_removed<do_T_derivatives>(tfactors, rate, drate_dT);
    rate_eval.screened_rates(k_n_Ni57_to_Ni58_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_n_Ni57_to_Ni58_removed) = drate_dT;

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
    rate_O16_to_He4_C12_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_O16_to_He4_C12_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_O16_to_He4_C12_derived) = drate_dT;

    }
    rate_Ne20_to_He4_O16_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Ne20_to_He4_O16_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ne20_to_He4_O16_derived) = drate_dT;

    }
    rate_Mg24_to_p_Na23_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Mg24_to_p_Na23_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Mg24_to_p_Na23_derived) = drate_dT;

    }
    rate_Mg24_to_He4_Ne20_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Mg24_to_He4_Ne20_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Mg24_to_He4_Ne20_derived) = drate_dT;

    }
    rate_Si28_to_p_Al27_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Si28_to_p_Al27_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Si28_to_p_Al27_derived) = drate_dT;

    }
    rate_Si28_to_He4_Mg24_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Si28_to_He4_Mg24_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Si28_to_He4_Mg24_derived) = drate_dT;

    }
    rate_S32_to_p_P31_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_S32_to_p_P31_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_S32_to_p_P31_derived) = drate_dT;

    }
    rate_S32_to_He4_Si28_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_S32_to_He4_Si28_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_S32_to_He4_Si28_derived) = drate_dT;

    }
    rate_Fe52_to_p_Mn51_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Fe52_to_p_Mn51_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe52_to_p_Mn51_derived) = drate_dT;

    }
    rate_Fe52_to_He4_Cr48_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Fe52_to_He4_Cr48_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe52_to_He4_Cr48_derived) = drate_dT;

    }
    rate_Co55_to_He4_Mn51_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Co55_to_He4_Mn51_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Co55_to_He4_Mn51_derived) = drate_dT;

    }
    rate_Ni56_to_p_Co55_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Ni56_to_p_Co55_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni56_to_p_Co55_derived) = drate_dT;

    }
    rate_Ni56_to_He4_Fe52_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Ni56_to_He4_Fe52_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni56_to_He4_Fe52_derived) = drate_dT;

    }
    rate_C12_to_He4_He4_He4_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_C12_to_He4_He4_He4_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_C12_to_He4_He4_He4_derived) = drate_dT;

    }
    rate_p_O16_to_He4_N13_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_p_O16_to_He4_N13_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_O16_to_He4_N13_derived) = drate_dT;

    }
    rate_He4_Ne20_to_p_Na23_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_He4_Ne20_to_p_Na23_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Ne20_to_p_Na23_derived) = drate_dT;

    }
    rate_He4_Ne20_to_C12_C12_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_He4_Ne20_to_C12_C12_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Ne20_to_C12_C12_derived) = drate_dT;

    }
    rate_p_Na23_to_C12_C12_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_p_Na23_to_C12_C12_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Na23_to_C12_C12_derived) = drate_dT;

    }
    rate_He4_Mg24_to_p_Al27_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_He4_Mg24_to_p_Al27_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Mg24_to_p_Al27_derived) = drate_dT;

    }
    rate_He4_Mg24_to_C12_O16_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_He4_Mg24_to_C12_O16_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Mg24_to_C12_O16_derived) = drate_dT;

    }
    rate_p_Al27_to_C12_O16_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_p_Al27_to_C12_O16_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Al27_to_C12_O16_derived) = drate_dT;

    }
    rate_He4_Si28_to_p_P31_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_He4_Si28_to_p_P31_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Si28_to_p_P31_derived) = drate_dT;

    }
    rate_He4_Si28_to_O16_O16_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_He4_Si28_to_O16_O16_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Si28_to_O16_O16_derived) = drate_dT;

    }
    rate_p_P31_to_O16_O16_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_p_P31_to_O16_O16_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_P31_to_O16_O16_derived) = drate_dT;

    }
    rate_p_Mn51_to_He4_Cr48_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_p_Mn51_to_He4_Cr48_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Mn51_to_He4_Cr48_derived) = drate_dT;

    }
    rate_p_Co55_to_He4_Fe52_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_p_Co55_to_He4_Fe52_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Co55_to_He4_Fe52_derived) = drate_dT;

    }
    rate_Mg24_to_C12_C12_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Mg24_to_C12_C12_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Mg24_to_C12_C12_derived) = drate_dT;

    }
    rate_S32_to_O16_O16_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_S32_to_O16_O16_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_S32_to_O16_O16_derived) = drate_dT;

    }
    rate_Si28_to_C12_O16_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Si28_to_C12_O16_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Si28_to_C12_O16_derived) = drate_dT;

    }
    rate_Co55_to_p_Fe54_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Co55_to_p_Fe54_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Co55_to_p_Fe54_derived) = drate_dT;

    }
    rate_Co56_to_n_Co55_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Co56_to_n_Co55_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Co56_to_n_Co55_derived) = drate_dT;

    }
    rate_Co57_to_n_Co56_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Co57_to_n_Co56_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Co57_to_n_Co56_derived) = drate_dT;

    }
    rate_Co57_to_p_Fe56_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Co57_to_p_Fe56_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Co57_to_p_Fe56_derived) = drate_dT;

    }
    rate_Ni58_to_p_Co57_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Ni58_to_p_Co57_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni58_to_p_Co57_derived) = drate_dT;

    }
    rate_Ni58_to_He4_Fe54_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Ni58_to_He4_Fe54_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni58_to_He4_Fe54_derived) = drate_dT;

    }
    rate_p_Fe54_to_He4_Mn51_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_p_Fe54_to_He4_Mn51_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Fe54_to_He4_Mn51_derived) = drate_dT;

    }
    rate_He4_Fe54_to_p_Co57_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_He4_Fe54_to_p_Co57_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Fe54_to_p_Co57_derived) = drate_dT;

    }
    rate_p_Fe56_to_n_Co56_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_p_Fe56_to_n_Co56_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Fe56_to_n_Co56_derived) = drate_dT;

    }
    rate_p_Co56_to_n_Ni56_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_p_Co56_to_n_Ni56_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Co56_to_n_Ni56_derived) = drate_dT;

    }
    rate_p_Ni58_to_He4_Co55_derived<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_p_Ni58_to_He4_Co55_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Ni58_to_He4_Co55_derived) = drate_dT;

    }
    rate_He4_S32_to_p_Cl35_removed<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_He4_S32_to_p_Cl35_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_S32_to_p_Cl35_removed) = drate_dT;

    }
    rate_Ar36_to_He4_S32_removed<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Ar36_to_He4_S32_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ar36_to_He4_S32_removed) = drate_dT;

    }
    rate_Ar36_to_p_Cl35_removed<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Ar36_to_p_Cl35_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ar36_to_p_Cl35_removed) = drate_dT;

    }
    rate_He4_Ar36_to_p_K39_removed<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_He4_Ar36_to_p_K39_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Ar36_to_p_K39_removed) = drate_dT;

    }
    rate_Ca40_to_He4_Ar36_removed<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Ca40_to_He4_Ar36_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ca40_to_He4_Ar36_removed) = drate_dT;

    }
    rate_Ca40_to_p_K39_removed<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Ca40_to_p_K39_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ca40_to_p_K39_removed) = drate_dT;

    }
    rate_He4_Ca40_to_p_Sc43_removed<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_He4_Ca40_to_p_Sc43_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Ca40_to_p_Sc43_removed) = drate_dT;

    }
    rate_Ti44_to_He4_Ca40_removed<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Ti44_to_He4_Ca40_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ti44_to_He4_Ca40_removed) = drate_dT;

    }
    rate_Ti44_to_p_Sc43_removed<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Ti44_to_p_Sc43_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ti44_to_p_Sc43_removed) = drate_dT;

    }
    rate_Cr48_to_He4_Ti44_removed<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Cr48_to_He4_Ti44_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Cr48_to_He4_Ti44_removed) = drate_dT;

    }
    rate_Cr48_to_p_V47_removed<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Cr48_to_p_V47_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Cr48_to_p_V47_removed) = drate_dT;

    }
    rate_p_V47_to_He4_Ti44_removed<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_p_V47_to_He4_Ti44_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_V47_to_He4_Ti44_removed) = drate_dT;

    }
    rate_Fe54_to_n_Fe53_removed<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Fe54_to_n_Fe53_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe54_to_n_Fe53_removed) = drate_dT;

    }
    rate_Fe53_to_n_Fe52_removed<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Fe53_to_n_Fe52_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe53_to_n_Fe52_removed) = drate_dT;

    }
    rate_Fe56_to_n_Fe55_removed<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Fe56_to_n_Fe55_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe56_to_n_Fe55_removed) = drate_dT;

    }
    rate_Fe55_to_n_Fe54_removed<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Fe55_to_n_Fe54_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe55_to_n_Fe54_removed) = drate_dT;

    }
    rate_Ni58_to_n_Ni57_removed<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Ni58_to_n_Ni57_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni58_to_n_Ni57_removed) = drate_dT;

    }
    rate_Ni57_to_n_Ni56_removed<T>(rate_eval, tfactors, rate, drate_dT, pf_cache);
    rate_eval.screened_rates(k_Ni57_to_n_Ni56_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni57_to_n_Ni56_removed) = drate_dT;

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
    rate_Fe52_n_n_to_Fe54_approx<T>(rate_eval, rho, Y, rate, drate_dT);
    rate_eval.screened_rates(k_Fe52_n_n_to_Fe54_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe52_n_n_to_Fe54_approx) = drate_dT;

    }
    rate_Fe54_to_Fe52_n_n_approx<T>(rate_eval, rho, Y, rate, drate_dT);
    rate_eval.screened_rates(k_Fe54_to_Fe52_n_n_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe54_to_Fe52_n_n_approx) = drate_dT;

    }
    rate_Fe54_n_n_to_Fe56_approx<T>(rate_eval, rho, Y, rate, drate_dT);
    rate_eval.screened_rates(k_Fe54_n_n_to_Fe56_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe54_n_n_to_Fe56_approx) = drate_dT;

    }
    rate_Fe56_to_Fe54_n_n_approx<T>(rate_eval, rho, Y, rate, drate_dT);
    rate_eval.screened_rates(k_Fe56_to_Fe54_n_n_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe56_to_Fe54_n_n_approx) = drate_dT;

    }
    rate_Ni56_n_n_to_Ni58_approx<T>(rate_eval, rho, Y, rate, drate_dT);
    rate_eval.screened_rates(k_Ni56_n_n_to_Ni58_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni56_n_n_to_Ni58_approx) = drate_dT;

    }
    rate_Ni58_to_Ni56_n_n_approx<T>(rate_eval, rho, Y, rate, drate_dT);
    rate_eval.screened_rates(k_Ni58_to_Ni56_n_n_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni58_to_Ni56_n_n_approx) = drate_dT;

    }

}

#endif
