#ifndef REACLIB_RATES_H
#define REACLIB_RATES_H

#include <AMReX.H>
#include <AMReX_Print.H>

#include <tfactors.H>
#include <actual_network.H>
#include <rate_type.H>

using namespace Rates;
using namespace Species;

template <int do_T_derivatives>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_C14_to_N14_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // C14 --> N14

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // wc12w
    ln_set_rate =  -26.2827;
    amrex::ignore_unused(tfactors);

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 = 0.0;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_N13_to_C13_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // N13 --> C13

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // wc12w
    ln_set_rate =  -6.7601;
    amrex::ignore_unused(tfactors);

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 = 0.0;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_C12_to_N13_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // C12 + p --> N13

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ls09n
    ln_set_rate =  17.1482 + -13.692 * tfactors.T913i + -0.230881 * tfactors.T913
                         + 4.44362 * tfactors.T9 + -3.15898 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -13.692 * tfactors.T943i + (1.0/3.0) * -0.230881 * tfactors.T923i
                                  + 4.44362 + (5.0/3.0) * -3.15898 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  3.77849 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -5.10735 * tfactors.T943i + (1.0/3.0) * -2.24111 * tfactors.T923i
                                  + 0.148883 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_C12_to_O16_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // C12 + He4 --> O16

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // nac2 
    ln_set_rate =  254.634 + -1.84097 * tfactors.T9i + 103.411 * tfactors.T913i + -420.567 * tfactors.T913
                         + 64.0874 * tfactors.T9 + -12.4624 * tfactors.T953 + 137.303 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.84097 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 103.411 * tfactors.T943i + (1.0/3.0) * -420.567 * tfactors.T923i
                                  + 64.0874 + (5.0/3.0) * -12.4624 * tfactors.T923 + 137.303 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.39254 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 58.9128 * tfactors.T943i + (1.0/3.0) * -148.273 * tfactors.T923i
                                  + 9.08324 + (5.0/3.0) * -0.541041 * tfactors.T923 + 70.3554 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_C13_to_N14_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // C13 + p --> N14

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // nacrn
    ln_set_rate =  18.5155 + -13.72 * tfactors.T913i + -0.450018 * tfactors.T913
                         + 3.70823 * tfactors.T9 + -1.70545 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -13.72 * tfactors.T943i + (1.0/3.0) * -0.450018 * tfactors.T923i
                                  + 3.70823 + (5.0/3.0) * -1.70545 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  5.78147 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.196703 * tfactors.T923i
                                  + 0.142126 + (5.0/3.0) * -0.0238912 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  13.5543 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_C14_to_N15_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // C14 + p --> N15

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10r
    ln_set_rate =  17.1024 + -1.52341 * tfactors.T9i + -10.658 * tfactors.T913i + 1.73644 * tfactors.T913
                         + -0.350498 * tfactors.T9 + 0.0279902 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.52341 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -10.658 * tfactors.T943i + (1.0/3.0) * 1.73644 * tfactors.T923i
                                  + -0.350498 + (5.0/3.0) * 0.0279902 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

    // il10n
    ln_set_rate =  20.119 + -13.9619 * tfactors.T913i + -4.34315 * tfactors.T913
                         + 6.64922 * tfactors.T9 + -3.22592 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -13.9619 * tfactors.T943i + (1.0/3.0) * -4.34315 * tfactors.T923i
                                  + 6.64922 + (5.0/3.0) * -3.22592 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_C14_to_O18_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // C14 + He4 --> O18

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10r
    ln_set_rate =  11.8309 + -10.3983 * tfactors.T9i + -3.83188 * tfactors.T913
                         + 1.64358 * tfactors.T9 + -0.177785 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  10.3983 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -3.83188 * tfactors.T923i
                                  + 1.64358 + (5.0/3.0) * -0.177785 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

    // il10r
    ln_set_rate =  -23.805 + -2.06876 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  2.06876 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

    // il10n
    ln_set_rate =  18.4877 + -31.7222 * tfactors.T913i + 11.3923 * tfactors.T913
                         + -9.92249 * tfactors.T9 + -2.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -31.7222 * tfactors.T943i + (1.0/3.0) * 11.3923 * tfactors.T923i
                                  + -9.92249 + (5.0/3.0) * -2.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_N14_to_F18_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // N14 + He4 --> F18

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10r
    ln_set_rate =  13.8995 + -10.9656 * tfactors.T9i + -5.6227 * tfactors.T913i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  10.9656 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -5.6227 * tfactors.T943i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  5.16034 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -36.2504 * tfactors.T943i
                                  + (5.0/3.0) * -5.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_N15_to_O16_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // N15 + p --> O16

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // li10r
    ln_set_rate =  14.5444 + -10.2295 * tfactors.T9i
                         + 0.0459037 * tfactors.T9 + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  10.2295 * tfactors.T9i * tfactors.T9i
                                  + 0.0459037 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  2.92315 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -15.24 * tfactors.T943i + (1.0/3.0) * 0.334926 * tfactors.T923i
                                  + 4.59088 + (5.0/3.0) * -4.78468 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_N15_to_F19_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // N15 + He4 --> F19

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10r
    ln_set_rate =  -28.7989 + -4.19986 * tfactors.T9i + 35.4292 * tfactors.T913
                         + -5.5767 * tfactors.T9 + 0.441293 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  4.19986 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 35.4292 * tfactors.T923i
                                  + -5.5767 + (5.0/3.0) * 0.441293 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  6.98462 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  4.17795 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -36.2324 * tfactors.T943i
                                  + (5.0/3.0) * -2.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_O16_to_Ne20_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // O16 + He4 --> Ne20

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // co10r
    ln_set_rate =  3.88571 + -10.3585 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  10.3585 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -39.7262 * tfactors.T943i + (1.0/3.0) * -0.210799 * tfactors.T923i
                                  + 0.442879 + (5.0/3.0) * -0.0797753 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  12.7643 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -3.65925 * tfactors.T923i
                                  + 0.714224 + (5.0/3.0) * -0.00107508 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_O17_to_F18_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // O17 + p --> F18

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10r
    ln_set_rate =  9.39048 + -6.22828 * tfactors.T9i + 2.31435 * tfactors.T913
                         + -0.302835 * tfactors.T9 + 0.020133 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  6.22828 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 2.31435 * tfactors.T923i
                                  + -0.302835 + (5.0/3.0) * 0.020133 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.746296 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -16.4035 * tfactors.T943i + (1.0/3.0) * 4.31885 * tfactors.T923i
                                  + -0.709921 + (5.0/3.0) * -2.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_O17_to_Ne21_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // O17 + He4 --> Ne21

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // be13r
    ln_set_rate =  -25.0898 + -5.50926 * tfactors.T9i + 123.363 * tfactors.T913i + -87.4351 * tfactors.T913
                         + -3.40974e-06 * tfactors.T9 + -57.0469 * tfactors.T953 + 82.2218 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  5.50926 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 123.363 * tfactors.T943i + (1.0/3.0) * -87.4351 * tfactors.T923i
                                  + -3.40974e-06 + (5.0/3.0) * -57.0469 * tfactors.T923 + 82.2218 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

    // be13r
    ln_set_rate =  -117.134 + -13.6759 * tfactors.T9i + 3.31162e-08 * tfactors.T913i + 130.258 * tfactors.T913
                         + -7.92551e-05 * tfactors.T9 + -4.13772 * tfactors.T953 + -42.7753 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  13.6759 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 3.31162e-08 * tfactors.T943i + (1.0/3.0) * 130.258 * tfactors.T923i
                                  + -7.92551e-05 + (5.0/3.0) * -4.13772 * tfactors.T923 + -42.7753 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

    // be13r
    ln_set_rate =  2.14 + -5.99952 * tfactors.T9i + 2.87641 * tfactors.T913i + -3.54489 * tfactors.T913
                         + -2.11222e-08 * tfactors.T9 + -3.90649e-09 * tfactors.T953 + 4.75778 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  5.99952 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 2.87641 * tfactors.T943i + (1.0/3.0) * -3.54489 * tfactors.T923i
                                  + -2.11222e-08 + (5.0/3.0) * -3.90649e-09 * tfactors.T923 + 4.75778 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_O18_to_F19_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // O18 + p --> F19

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10n
    ln_set_rate =  19.917 + -16.7246 * tfactors.T913i
                         + -3.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -16.7246 * tfactors.T943i
                                  + (5.0/3.0) * -3.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  6.7253 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 3.99059 * tfactors.T923i
                                  + -0.593127 + (5.0/3.0) * 0.0877534 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.65681 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.244743 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_O18_to_Ne22_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // O18 + He4 --> Ne22

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10r
    ln_set_rate =  -31.9126 + -2.00306 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  2.00306 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

    // il10r
    ln_set_rate =  -81.3036 + -0.676112 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.676112 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

    // il10r
    ln_set_rate =  14.9748 + -31.0468 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  31.0468 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

    // il10r
    ln_set_rate =  82.2053 + -1.58534 * tfactors.T9i + -44.3823 * tfactors.T913i + -46.6617 * tfactors.T913
                         + 7.88059 * tfactors.T9 + -0.590829 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.58534 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -44.3823 * tfactors.T943i + (1.0/3.0) * -46.6617 * tfactors.T923i
                                  + 7.88059 + (5.0/3.0) * -0.590829 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_F18_to_Na22_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // F18 + He4 --> Na22

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // rpsmr
    ln_set_rate =  35.3786 + -1.82957 * tfactors.T9i + 18.8956 * tfactors.T913i + -65.6134 * tfactors.T913
                         + 1.71114 * tfactors.T9 + -0.0260999 * tfactors.T953 + 37.8396 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.82957 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 18.8956 * tfactors.T943i + (1.0/3.0) * -65.6134 * tfactors.T923i
                                  + 1.71114 + (5.0/3.0) * -0.0260999 * tfactors.T923 + 37.8396 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_F19_to_Ne20_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // F19 + p --> Ne20

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // nacrr
    ln_set_rate =  -5.63093 + -7.74414 * tfactors.T9i + 31.6442 * tfactors.T913i + -58.6563 * tfactors.T913
                         + 67.7365 * tfactors.T9 + -22.9721 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  7.74414 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 31.6442 * tfactors.T943i + (1.0/3.0) * -58.6563 * tfactors.T923i
                                  + 67.7365 + (5.0/3.0) * -22.9721 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.71383 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -11.3832 * tfactors.T943i + (1.0/3.0) * 5.47872 * tfactors.T923i
                                  + -1.07203 + (5.0/3.0) * 0.11196 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -18.116 * tfactors.T943i + (1.0/3.0) * -1.4622 * tfactors.T923i
                                  + 6.95113 + (5.0/3.0) * -2.90366 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_F19_to_Na23_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // F19 + He4 --> Na23

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // rpsmr
    ln_set_rate =  52.7856 + -2.11408 * tfactors.T9i + 39.7219 * tfactors.T913i + -100.401 * tfactors.T913
                         + 3.15808 * tfactors.T9 + -0.0629822 * tfactors.T953 + 54.4823 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  2.11408 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 39.7219 * tfactors.T943i + (1.0/3.0) * -100.401 * tfactors.T923i
                                  + 3.15808 + (5.0/3.0) * -0.0629822 * tfactors.T923 + 54.4823 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Ne20_to_Mg24_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ne20 + He4 --> Mg24

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10r
    ln_set_rate =  -8.79827 + -12.7809 * tfactors.T9i + 16.9229 * tfactors.T913
                         + -2.57325 * tfactors.T9 + 0.208997 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  12.7809 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 16.9229 * tfactors.T923i
                                  + -2.57325 + (5.0/3.0) * 0.208997 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  9.22026 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  2.50605 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -46.2525 * tfactors.T943i + (1.0/3.0) * 5.58901 * tfactors.T923i
                                  + 7.61843 + (5.0/3.0) * -3.683 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Ne21_to_Na22_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ne21 + p --> Na22

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10n
    ln_set_rate =  19.0696 + -19.2096 * tfactors.T913i
                         + -1.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -19.2096 * tfactors.T943i
                                  + (5.0/3.0) * -1.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

    // il10r
    ln_set_rate =  -39.4862 + -4.21385 * tfactors.T9i + 21.1176 * tfactors.T913i + 34.0411 * tfactors.T913
                         + -4.45593 * tfactors.T9 + 0.328613 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  4.21385 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 21.1176 * tfactors.T943i + (1.0/3.0) * 34.0411 * tfactors.T923i
                                  + -4.45593 + (5.0/3.0) * 0.328613 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

    // il10r
    ln_set_rate =  1.75704 + -1.39957 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.39957 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

    // il10r
    ln_set_rate =  -47.6554 + -0.19618 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.19618 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Ne21_to_Mg25_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ne21 + He4 --> Mg25

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // cf88r
    ln_set_rate =  26.2429 + -22.049 * tfactors.T9i + -29.4583 * tfactors.T913
                         + 14.6328 * tfactors.T9 + -3.47392 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  22.049 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -29.4583 * tfactors.T923i
                                  + 14.6328 + (5.0/3.0) * -3.47392 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

    // cf88n
    ln_set_rate =  36.6927 + -46.89 * tfactors.T913i + -0.72642 * tfactors.T913
                         + -0.76406 * tfactors.T9 + 0.0797483 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -46.89 * tfactors.T943i + (1.0/3.0) * -0.72642 * tfactors.T923i
                                  + -0.76406 + (5.0/3.0) * 0.0797483 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Ne22_to_Na23_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ne22 + p --> Na23

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ke17r
    ln_set_rate =  -0.599331 + -1.72007 * tfactors.T9i
                         + -0.310765 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.72007 * tfactors.T9i * tfactors.T9i
                                  + -0.310765 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

    // ke17r
    ln_set_rate =  -21.4345 + -0.410962 * tfactors.T9i
                         + -1.49019 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.410962 * tfactors.T9i * tfactors.T9i
                                  + -1.49019 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

    // ke17r
    ln_set_rate =  -4.09035 + -0.799756 * tfactors.T9i
                         + 3.23558 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.799756 * tfactors.T9i * tfactors.T9i
                                  + 3.23558 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

    // ke17r
    ln_set_rate =  14.0373 + -8.72377 * tfactors.T9i
                         + -0.767467 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  8.72377 * tfactors.T9i * tfactors.T9i
                                  + -0.767467 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

    // ke17r
    ln_set_rate =  11.6408 + -4.59936 * tfactors.T9i
                         + 0.156226 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  4.59936 * tfactors.T9i * tfactors.T9i
                                  + 0.156226 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

    // ke17r
    ln_set_rate =  -4.00597 + -2.6179 * tfactors.T9i
                         + -4.29964 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  2.6179 * tfactors.T9i * tfactors.T9i
                                  + -4.29964 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Na22_to_Al26_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Na22 + He4 --> Al26

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  36.3797 + -50.0924 * tfactors.T913i + -0.390826 * tfactors.T913
                         + -0.99531 * tfactors.T9 + 0.101354 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -50.0924 * tfactors.T943i + (1.0/3.0) * -0.390826 * tfactors.T923i
                                  + -0.99531 + (5.0/3.0) * 0.101354 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Na23_to_Mg24_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Na23 + p --> Mg24

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10r
    ln_set_rate =  9.0594 + -3.28029 * tfactors.T9i + -0.360588 * tfactors.T913
                         + 1.4187 * tfactors.T9 + -0.184061 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  3.28029 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.360588 * tfactors.T923i
                                  + 1.4187 + (5.0/3.0) * -0.184061 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.61219 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -20.6428 * tfactors.T943i + (1.0/3.0) * 1.52954 * tfactors.T923i
                                  + 2.7487 + (5.0/3.0) * -1.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Na23_to_Al27_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Na23 + He4 --> Al27

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  44.7724 + -50.2042 * tfactors.T913i + -1.64239 * tfactors.T913
                         + -1.59995 * tfactors.T9 + 0.184933 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -50.2042 * tfactors.T943i + (1.0/3.0) * -1.64239 * tfactors.T923i
                                  + -1.59995 + (5.0/3.0) * 0.184933 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Mg24_to_Si28_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Mg24 + He4 --> Si28

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // st08r
    ln_set_rate =  8.03977 + -15.629 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  15.629 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  12.8332 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 21.3721 * tfactors.T943i + (1.0/3.0) * 37.7649 * tfactors.T923i
                                  + -4.10635 + (5.0/3.0) * 0.249618 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Mg25_to_Al26_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Mg25 + p --> Al26

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10r
    ln_set_rate =  4.21826 + -0.71983 * tfactors.T9i + -88.9297 * tfactors.T913
                         + 302.948 * tfactors.T9 + -346.461 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.71983 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -88.9297 * tfactors.T923i
                                  + 302.948 + (5.0/3.0) * -346.461 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

    // il10r
    ln_set_rate =  -37.1963 + -0.429366 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.429366 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

    // il10r
    ln_set_rate =  2.22778 + -3.22353 * tfactors.T9i + 8.46334 * tfactors.T913
                         + -0.907024 * tfactors.T9 + 0.0642981 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  3.22353 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 8.46334 * tfactors.T923i
                                  + -0.907024 + (5.0/3.0) * 0.0642981 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Mg25_to_Si29_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Mg25 + He4 --> Si29

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // cf88n
    ln_set_rate =  40.3715 + -53.41 * tfactors.T913i + -1.83266 * tfactors.T913
                         + -0.573073 * tfactors.T9 + 0.0462678 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -53.41 * tfactors.T943i + (1.0/3.0) * -1.83266 * tfactors.T923i
                                  + -0.573073 + (5.0/3.0) * 0.0462678 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Al26_to_P30_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Al26 + He4 --> P30

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  42.9778 + -56.4422 * tfactors.T913i + -2.44848 * tfactors.T913
                         + -1.17578 * tfactors.T9 + 0.150757 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -56.4422 * tfactors.T943i + (1.0/3.0) * -2.44848 * tfactors.T923i
                                  + -1.17578 + (5.0/3.0) * 0.150757 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Al27_to_Si28_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Al27 + p --> Si28

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10r
    ln_set_rate =  -13.6664 + -1.90396 * tfactors.T9i + 23.8634 * tfactors.T913
                         + -3.70135 * tfactors.T9 + 0.28964 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.90396 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 23.8634 * tfactors.T923i
                                  + -3.70135 + (5.0/3.0) * 0.28964 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.387313 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -26.8327 * tfactors.T943i + (1.0/3.0) * -116.137 * tfactors.T923i
                                  + 0.00950567 + (5.0/3.0) * 0.00999755 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -23.2205 * tfactors.T943i
                                  + (5.0/3.0) * -2.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Al27_to_P31_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Al27 + He4 --> P31

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  47.2333 + -56.5351 * tfactors.T913i + -0.896208 * tfactors.T913
                         + -1.72024 * tfactors.T9 + 0.185409 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -56.5351 * tfactors.T943i + (1.0/3.0) * -0.896208 * tfactors.T923i
                                  + -1.72024 + (5.0/3.0) * 0.185409 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Si28_to_S32_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Si28 + He4 --> S32

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  47.9212 + -59.4896 * tfactors.T913i + 4.47205 * tfactors.T913
                         + -4.78989 * tfactors.T9 + 0.557201 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -59.4896 * tfactors.T943i + (1.0/3.0) * 4.47205 * tfactors.T923i
                                  + -4.78989 + (5.0/3.0) * 0.557201 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Si29_to_P30_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Si29 + p --> P30

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10r
    ln_set_rate =  -1.24791 + -3.33929 * tfactors.T9i + 14.0921 * tfactors.T913
                         + -3.92096 * tfactors.T9 + 0.447706 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  3.33929 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 14.0921 * tfactors.T923i
                                  + -3.92096 + (5.0/3.0) * 0.447706 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

    // il10r
    ln_set_rate =  -13.4701 + -1.25026 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.25026 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

    // il10n
    ln_set_rate =  16.5182 + -23.9101 * tfactors.T913i + 10.7796 * tfactors.T913
                         + -3.04181 * tfactors.T9 + 0.274565 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -23.9101 * tfactors.T943i + (1.0/3.0) * 10.7796 * tfactors.T923i
                                  + -3.04181 + (5.0/3.0) * 0.274565 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Si29_to_S33_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Si29 + He4 --> S33

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  49.5657 + -59.5755 * tfactors.T913i + 1.06274 * tfactors.T913
                         + -3.07529 * tfactors.T9 + 0.372011 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -59.5755 * tfactors.T943i + (1.0/3.0) * 1.06274 * tfactors.T923i
                                  + -3.07529 + (5.0/3.0) * 0.372011 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_P30_to_Cl34_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // P30 + He4 --> Cl34

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  45.3321 + -62.4643 * tfactors.T913i + -3.19028 * tfactors.T913
                         + -0.832633 * tfactors.T9 + 0.0987525 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -62.4643 * tfactors.T943i + (1.0/3.0) * -3.19028 * tfactors.T923i
                                  + -0.832633 + (5.0/3.0) * 0.0987525 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_P31_to_S32_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // P31 + p --> S32

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10r
    ln_set_rate =  0.821556 + -3.77704 * tfactors.T9i + 8.09341 * tfactors.T913
                         + -0.615971 * tfactors.T9 + 0.031159 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  3.77704 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 8.09341 * tfactors.T923i
                                  + -0.615971 + (5.0/3.0) * 0.031159 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  2.25958 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -25.3278 * tfactors.T943i + (1.0/3.0) * 6.4931 * tfactors.T923i
                                  + -9.27513 + (5.0/3.0) * -0.610439 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_P31_to_Cl35_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // P31 + He4 --> Cl35

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  50.451 + -62.5433 * tfactors.T913i + -2.95026 * tfactors.T913
                         + -0.89652 * tfactors.T9 + 0.0774126 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -62.5433 * tfactors.T943i + (1.0/3.0) * -2.95026 * tfactors.T923i
                                  + -0.89652 + (5.0/3.0) * 0.0774126 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_S32_to_Ar36_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // S32 + He4 --> Ar36

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  48.901 + -65.3709 * tfactors.T913i + 5.68294 * tfactors.T913
                         + -5.00388 * tfactors.T9 + 0.571407 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -65.3709 * tfactors.T943i + (1.0/3.0) * 5.68294 * tfactors.T923i
                                  + -5.00388 + (5.0/3.0) * 0.571407 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_S33_to_Cl34_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // S33 + p --> Cl34

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  36.4908 + -26.777 * tfactors.T913i + -5.96882 * tfactors.T913
                         + -1.0706 * tfactors.T9 + 0.19692 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -26.777 * tfactors.T943i + (1.0/3.0) * -5.96882 * tfactors.T923i
                                  + -1.0706 + (5.0/3.0) * 0.19692 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_S33_to_Ar37_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // S33 + He4 --> Ar37

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  49.9315 + -65.4446 * tfactors.T913i + 3.59607 * tfactors.T913
                         + -3.40501 * tfactors.T9 + 0.363961 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -65.4446 * tfactors.T943i + (1.0/3.0) * 3.59607 * tfactors.T923i
                                  + -3.40501 + (5.0/3.0) * 0.363961 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Cl34_to_K38_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Cl34 + He4 --> K38

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  47.9589 + -68.2165 * tfactors.T913i + -2.18792 * tfactors.T913
                         + -1.25322 * tfactors.T9 + 0.140901 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -68.2165 * tfactors.T943i + (1.0/3.0) * -2.18792 * tfactors.T923i
                                  + -1.25322 + (5.0/3.0) * 0.140901 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Cl35_to_Ar36_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Cl35 + p --> Ar36

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10r
    ln_set_rate =  -42.5249 + -0.564651 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.564651 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -27.8971 * tfactors.T943i + (1.0/3.0) * -16.2304 * tfactors.T923i
                                  + 35.255 + (5.0/3.0) * -25.8411 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  3.65092 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 18.0179 * tfactors.T923i
                                  + -2.86304 + (5.0/3.0) * 0.250854 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  2.00996 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Cl35_to_K39_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Cl35 + He4 --> K39

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  52.718 + -68.2848 * tfactors.T913i + 0.0178545 * tfactors.T913
                         + -2.06783 * tfactors.T9 + 0.199659 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -68.2848 * tfactors.T943i + (1.0/3.0) * 0.0178545 * tfactors.T923i
                                  + -2.06783 + (5.0/3.0) * 0.199659 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Ar36_to_Ca40_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ar36 + He4 --> Ca40

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  52.3486 + -71.0046 * tfactors.T913i + 4.0656 * tfactors.T913
                         + -5.26509 * tfactors.T9 + 0.683546 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -71.0046 * tfactors.T943i + (1.0/3.0) * 4.0656 * tfactors.T923i
                                  + -5.26509 + (5.0/3.0) * 0.683546 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Ar37_to_K38_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ar37 + p --> K38

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  37.7254 + -28.9954 * tfactors.T913i + -5.57966 * tfactors.T913
                         + -1.27516 * tfactors.T9 + 0.217568 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -28.9954 * tfactors.T943i + (1.0/3.0) * -5.57966 * tfactors.T923i
                                  + -1.27516 + (5.0/3.0) * 0.217568 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Ar37_to_Ca41_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ar37 + He4 --> Ca41

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  52.7522 + -71.0688 * tfactors.T913i + 2.71847 * tfactors.T913
                         + -3.17333 * tfactors.T9 + 0.335323 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -71.0688 * tfactors.T943i + (1.0/3.0) * 2.71847 * tfactors.T923i
                                  + -3.17333 + (5.0/3.0) * 0.335323 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_K38_to_Sc42_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // K38 + He4 --> Sc42

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  53.2593 + -73.7407 * tfactors.T913i + -7.08848 * tfactors.T913
                         + 0.322473 * tfactors.T9 + -0.0523557 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -73.7407 * tfactors.T943i + (1.0/3.0) * -7.08848 * tfactors.T923i
                                  + 0.322473 + (5.0/3.0) * -0.0523557 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_K39_to_Ca40_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // K39 + p --> Ca40

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // lo18r
    ln_set_rate =  2761.38 + -5.22234 * tfactors.T9i + 802.18 * tfactors.T913i + -4010.27 * tfactors.T913
                         + 1136.19 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  5.22234 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 802.18 * tfactors.T943i + (1.0/3.0) * -4010.27 * tfactors.T923i
                                  + 1136.19 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  12.5647 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 641.844 * tfactors.T943i + (1.0/3.0) * -1248.49 * tfactors.T923i
                                  + 564.926 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.66508 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 41.1723 * tfactors.T943i + (1.0/3.0) * -149.299 * tfactors.T923i
                                  + 10.5229 + (5.0/3.0) * -0.68208 * tfactors.T923 + 59.2367 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_K39_to_Sc43_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // K39 + He4 --> Sc43

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  54.1202 + -73.8006 * tfactors.T913i + 1.87885 * tfactors.T913
                         + -2.75862 * tfactors.T9 + 0.279678 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -73.8006 * tfactors.T943i + (1.0/3.0) * 1.87885 * tfactors.T923i
                                  + -2.75862 + (5.0/3.0) * 0.279678 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Ca40_to_Ti44_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ca40 + He4 --> Ti44

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // chw0 
    ln_set_rate =  53.75 + -76.4273 * tfactors.T913i + 3.87451 * tfactors.T913
                         + -3.61477 * tfactors.T9 + 0.367451 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -76.4273 * tfactors.T943i + (1.0/3.0) * 3.87451 * tfactors.T923i
                                  + -3.61477 + (5.0/3.0) * 0.367451 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Ca41_to_Sc42_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ca41 + p --> Sc42

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  39.315 + -31.1321 * tfactors.T913i + -6.33393 * tfactors.T913
                         + -1.39316 * tfactors.T9 + 0.255024 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -31.1321 * tfactors.T943i + (1.0/3.0) * -6.33393 * tfactors.T923i
                                  + -1.39316 + (5.0/3.0) * 0.255024 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Ca41_to_Ti45_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ca41 + He4 --> Ti45

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  54.6305 + -76.4839 * tfactors.T913i + 3.03748 * tfactors.T913
                         + -2.59814 * tfactors.T9 + 0.210582 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -76.4839 * tfactors.T943i + (1.0/3.0) * 3.03748 * tfactors.T923i
                                  + -2.59814 + (5.0/3.0) * 0.210582 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Sc42_to_V46_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Sc42 + He4 --> V46

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  51.7727 + -79.0687 * tfactors.T913i + -4.25493 * tfactors.T913
                         + 0.184703 * tfactors.T9 + -0.0785452 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -79.0687 * tfactors.T943i + (1.0/3.0) * -4.25493 * tfactors.T923i
                                  + 0.184703 + (5.0/3.0) * -0.0785452 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Sc43_to_Ti44_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Sc43 + p --> Ti44

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  36.8432 + -32.1734 * tfactors.T913i + -1.77078 * tfactors.T913
                         + -2.21706 * tfactors.T9 + 0.298499 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -32.1734 * tfactors.T943i + (1.0/3.0) * -1.77078 * tfactors.T923i
                                  + -2.21706 + (5.0/3.0) * 0.298499 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Sc43_to_V47_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Sc43 + He4 --> V47

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  59.0195 + -79.122 * tfactors.T913i + -7.07006 * tfactors.T913
                         + 0.424183 * tfactors.T9 + -0.0665231 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -79.122 * tfactors.T943i + (1.0/3.0) * -7.07006 * tfactors.T923i
                                  + 0.424183 + (5.0/3.0) * -0.0665231 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Ti44_to_Cr48_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ti44 + He4 --> Cr48

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  64.7958 + -81.667 * tfactors.T913i + -10.6333 * tfactors.T913
                         + -0.672613 * tfactors.T9 + 0.161209 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -81.667 * tfactors.T943i + (1.0/3.0) * -10.6333 * tfactors.T923i
                                  + -0.672613 + (5.0/3.0) * 0.161209 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Ti45_to_V46_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ti45 + p --> V46

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  37.4468 + -33.1981 * tfactors.T913i + -1.66837 * tfactors.T913
                         + -2.50285 * tfactors.T9 + 0.349152 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -33.1981 * tfactors.T943i + (1.0/3.0) * -1.66837 * tfactors.T923i
                                  + -2.50285 + (5.0/3.0) * 0.349152 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Ti45_to_Cr49_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ti45 + He4 --> Cr49

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  64.023 + -81.7175 * tfactors.T913i + -10.1755 * tfactors.T913
                         + 0.364167 * tfactors.T9 + -0.000317666 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -81.7175 * tfactors.T943i + (1.0/3.0) * -10.1755 * tfactors.T923i
                                  + 0.364167 + (5.0/3.0) * -0.000317666 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_V46_to_Mn50_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // V46 + He4 --> Mn50

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  54.311 + -84.2255 * tfactors.T913i + -4.85634 * tfactors.T913
                         + 0.0528515 * tfactors.T9 + -0.0425496 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -84.2255 * tfactors.T943i + (1.0/3.0) * -4.85634 * tfactors.T923i
                                  + 0.0528515 + (5.0/3.0) * -0.0425496 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_V47_to_Cr48_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // V47 + p --> Cr48

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // nfisn
    ln_set_rate =  42.6798 + -6.0593 * tfactors.T9i + -34.0548 * tfactors.T913i + -3.41973 * tfactors.T913
                         + 1.16501 * tfactors.T9 + -0.105543 * tfactors.T953 + -7.70886 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  6.0593 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -34.0548 * tfactors.T943i + (1.0/3.0) * -3.41973 * tfactors.T923i
                                  + 1.16501 + (5.0/3.0) * -0.105543 * tfactors.T923 + -7.70886 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  5.29491 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 317.171 * tfactors.T943i + (1.0/3.0) * -911.679 * tfactors.T923i
                                  + 94.4245 + (5.0/3.0) * -10.1973 * tfactors.T923 + 330.727 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  -0.246665 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -45.9868 * tfactors.T943i + (1.0/3.0) * 13.6822 * tfactors.T923i
                                  + -0.376902 + (5.0/3.0) * -0.0194875 * tfactors.T923 + -8.42325 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.514414 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -110.655 * tfactors.T943i + (1.0/3.0) * 83.0232 * tfactors.T923i
                                  + -19.7762 + (5.0/3.0) * 3.03961 * tfactors.T923 + -49.4742 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_V47_to_Mn51_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // V47 + He4 --> Mn51

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  56.8618 + -84.2732 * tfactors.T913i + -2.98061 * tfactors.T913
                         + -0.531361 * tfactors.T9 + 0.023612 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -84.2732 * tfactors.T943i + (1.0/3.0) * -2.98061 * tfactors.T923i
                                  + -0.531361 + (5.0/3.0) * 0.023612 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Cr48_to_Fe52_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Cr48 + He4 --> Fe52

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  65.1754 + -86.7459 * tfactors.T913i + -9.79373 * tfactors.T913
                         + -0.772169 * tfactors.T9 + 0.155883 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -86.7459 * tfactors.T943i + (1.0/3.0) * -9.79373 * tfactors.T923i
                                  + -0.772169 + (5.0/3.0) * 0.155883 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Cr49_to_Mn50_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Cr49 + p --> Mn50

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  36.3884 + -35.2018 * tfactors.T913i + 0.168579 * tfactors.T913
                         + -2.87983 * tfactors.T9 + 0.378768 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -35.2018 * tfactors.T943i + (1.0/3.0) * 0.168579 * tfactors.T923i
                                  + -2.87983 + (5.0/3.0) * 0.378768 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Mn51_to_Fe52_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Mn51 + p --> Fe52

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  36.2596 + -36.1825 * tfactors.T913i + 0.873042 * tfactors.T913
                         + -2.89731 * tfactors.T9 + 0.364394 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -36.1825 * tfactors.T943i + (1.0/3.0) * 0.873042 * tfactors.T923i
                                  + -2.89731 + (5.0/3.0) * 0.364394 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Mn51_to_Co55_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Mn51 + He4 --> Co55

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  65.9219 + -89.274 * tfactors.T913i + -10.4373 * tfactors.T913
                         + 1.00492 * tfactors.T9 + -0.125548 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -89.274 * tfactors.T943i + (1.0/3.0) * -10.4373 * tfactors.T923i
                                  + 1.00492 + (5.0/3.0) * -0.125548 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Fe52_to_Ni56_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe52 + He4 --> Ni56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  66.6417 + -91.6819 * tfactors.T913i + -9.51885 * tfactors.T913
                         + -0.533014 * tfactors.T9 + 0.0892607 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -91.6819 * tfactors.T943i + (1.0/3.0) * -9.51885 * tfactors.T923i
                                  + -0.533014 + (5.0/3.0) * 0.0892607 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Co55_to_Ni56_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co55 + p --> Ni56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  37.3736 + -38.1053 * tfactors.T913i + -0.210947 * tfactors.T913
                         + -2.68377 * tfactors.T9 + 0.355814 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -38.1053 * tfactors.T943i + (1.0/3.0) * -0.210947 * tfactors.T923i
                                  + -2.68377 + (5.0/3.0) * 0.355814 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_C12_C12_to_p_Na23_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // C12 + C12 --> p + Na23

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // cf88r
    ln_set_rate =  60.9649 + -84.165 * tfactors.T913i + -1.4191 * tfactors.T913
                         + -0.114619 * tfactors.T9 + -0.070307 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -84.165 * tfactors.T943i + (1.0/3.0) * -1.4191 * tfactors.T923i
                                  + -0.114619 + (5.0/3.0) * -0.070307 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_C12_C12_to_He4_Ne20_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // C12 + C12 --> He4 + Ne20

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // cf88r
    ln_set_rate =  61.2863 + -84.165 * tfactors.T913i + -1.56627 * tfactors.T913
                         + -0.0736084 * tfactors.T9 + -0.072797 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -84.165 * tfactors.T943i + (1.0/3.0) * -1.56627 * tfactors.T923i
                                  + -0.0736084 + (5.0/3.0) * -0.072797 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_N13_to_p_O16_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // N13 + He4 --> p + O16

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // cf88n
    ln_set_rate =  40.4644 + -35.829 * tfactors.T913i + -0.530275 * tfactors.T913
                         + -0.982462 * tfactors.T9 + 0.0808059 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -35.829 * tfactors.T943i + (1.0/3.0) * -0.530275 * tfactors.T923i
                                  + -0.982462 + (5.0/3.0) * 0.0808059 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_N15_to_He4_C12_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // N15 + p --> He4 + C12

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // nacrn
    ln_set_rate =  27.4764 + -15.253 * tfactors.T913i + 1.59318 * tfactors.T913
                         + 2.4479 * tfactors.T9 + -2.19708 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -15.253 * tfactors.T943i + (1.0/3.0) * 1.59318 * tfactors.T923i
                                  + 2.4479 + (5.0/3.0) * -2.19708 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.1638 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 22.7105 * tfactors.T923i
                                  + -2.90707 + (5.0/3.0) * 0.205754 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  7.406 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  2.02117 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 30.8497 * tfactors.T923i
                                  + -8.50433 + (5.0/3.0) * -1.54426 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_C12_O16_to_p_Al27_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // O16 + C12 --> p + Al27

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // cf88r
    ln_set_rate =  68.5253 + 0.205134 * tfactors.T9i + -119.242 * tfactors.T913i + 13.3667 * tfactors.T913
                         + 0.295425 * tfactors.T9 + -0.267288 * tfactors.T953 + -9.91729 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  -0.205134 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -119.242 * tfactors.T943i + (1.0/3.0) * 13.3667 * tfactors.T923i
                                  + 0.295425 + (5.0/3.0) * -0.267288 * tfactors.T923 + -9.91729 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_C12_O16_to_He4_Mg24_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // O16 + C12 --> He4 + Mg24

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // cf88r
    ln_set_rate =  48.5341 + 0.37204 * tfactors.T9i + -133.413 * tfactors.T913i + 50.1572 * tfactors.T913
                         + -3.15987 * tfactors.T9 + 0.0178251 * tfactors.T953 + -23.7027 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  -0.37204 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -133.413 * tfactors.T943i + (1.0/3.0) * 50.1572 * tfactors.T923i
                                  + -3.15987 + (5.0/3.0) * 0.0178251 * tfactors.T923 + -23.7027 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_O16_O16_to_p_P31_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // O16 + O16 --> p + P31

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // cf88r
    ln_set_rate =  85.2628 + 0.223453 * tfactors.T9i + -145.844 * tfactors.T913i + 8.72612 * tfactors.T913
                         + -0.554035 * tfactors.T9 + -0.137562 * tfactors.T953 + -6.88807 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  -0.223453 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -145.844 * tfactors.T943i + (1.0/3.0) * 8.72612 * tfactors.T923i
                                  + -0.554035 + (5.0/3.0) * -0.137562 * tfactors.T923 + -6.88807 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_O16_O16_to_He4_Si28_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // O16 + O16 --> He4 + Si28

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // cf88r
    ln_set_rate =  97.2435 + -0.268514 * tfactors.T9i + -119.324 * tfactors.T913i + -32.2497 * tfactors.T913
                         + 1.46214 * tfactors.T9 + -0.200893 * tfactors.T953 + 13.2148 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.268514 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -119.324 * tfactors.T943i + (1.0/3.0) * -32.2497 * tfactors.T923i
                                  + 1.46214 + (5.0/3.0) * -0.200893 * tfactors.T923 + 13.2148 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_O17_to_He4_N14_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // O17 + p --> He4 + N14

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10r
    ln_set_rate =  -7.20763 + -0.753395 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.753395 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -16.9078 * tfactors.T943i
                                  + (5.0/3.0) * -2.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  4.95865 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 5.10182 * tfactors.T923i
                                  + 0.379373 + (5.0/3.0) * -0.0672515 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  2.11477 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_O18_to_He4_N15_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // O18 + p --> He4 + N15

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10r
    ln_set_rate =  -27.9044 + -0.245884 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.245884 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -16.6979 * tfactors.T943i
                                  + (5.0/3.0) * -3.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  5.32335 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 11.6568 * tfactors.T923i
                                  + -2.16303 + (5.0/3.0) * 0.209965 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.663 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_F18_to_p_Ne21_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // F18 + He4 --> p + Ne21

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // rpsmr
    ln_set_rate =  49.7863 + -1.84559 * tfactors.T9i + 21.4461 * tfactors.T913i + -73.252 * tfactors.T913
                         + 2.42329 * tfactors.T9 + -0.077278 * tfactors.T953 + 40.7604 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.84559 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 21.4461 * tfactors.T943i + (1.0/3.0) * -73.252 * tfactors.T923i
                                  + 2.42329 + (5.0/3.0) * -0.077278 * tfactors.T923 + 40.7604 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_F19_to_He4_O16_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // F19 + p --> He4 + O16

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // nacr 
    ln_set_rate =  -52.7043 + -0.12765 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.12765 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -18.116 * tfactors.T943i
                                  + 1.86674 + (5.0/3.0) * -7.5666 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  3.286 * tfactors.T9i * tfactors.T9i
                                  + -0.21103 + 2.87702 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  3.75185 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  2.46828 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_F19_to_p_Ne22_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // F19 + He4 --> p + Ne22

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // da18r
    ln_set_rate =  52.9317 + -2.8444 * tfactors.T9i + -38.7722 * tfactors.T913i + -13.3654 * tfactors.T913
                         + 0.863648 * tfactors.T9 + -0.0451491 * tfactors.T953 + 1.33333 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  2.8444 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -38.7722 * tfactors.T943i + (1.0/3.0) * -13.3654 * tfactors.T923i
                                  + 0.863648 + (5.0/3.0) * -0.0451491 * tfactors.T923 + 1.33333 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

    // da18r
    ln_set_rate =  51.6709 + -45.7808 * tfactors.T9i + -34.5008 * tfactors.T913i + 56.9316 * tfactors.T913
                         + 2.09613 * tfactors.T9 + -32.496 * tfactors.T953 + 0.333333 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  45.7808 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -34.5008 * tfactors.T943i + (1.0/3.0) * 56.9316 * tfactors.T923i
                                  + 2.09613 + (5.0/3.0) * -32.496 * tfactors.T923 + 0.333333 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

    // da18r
    ln_set_rate =  29430.6 + -133.026 * tfactors.T9i + 12625.1 * tfactors.T913i + -49107.1 * tfactors.T913
                         + 9227.53 * tfactors.T9 + -2086.65 * tfactors.T953 + 14520.2 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  133.026 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 12625.1 * tfactors.T943i + (1.0/3.0) * -49107.1 * tfactors.T923i
                                  + 9227.53 + (5.0/3.0) * -2086.65 * tfactors.T923 + 14520.2 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_C12_Ne20_to_p_P31_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ne20 + C12 --> p + P31

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // rolfr
    ln_set_rate =  -268.136 + -38.7624 * tfactors.T9i + 361.154 * tfactors.T913i + -92.643 * tfactors.T913
                         + -9.98738 * tfactors.T9 + 0.892737 * tfactors.T953 + 161.042 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  38.7624 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 361.154 * tfactors.T943i + (1.0/3.0) * -92.643 * tfactors.T923i
                                  + -9.98738 + (5.0/3.0) * 0.892737 * tfactors.T923 + 161.042 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_C12_Ne20_to_He4_Si28_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ne20 + C12 --> He4 + Si28

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // rolfr
    ln_set_rate =  -308.905 + -47.2175 * tfactors.T9i + 514.197 * tfactors.T913i + -200.896 * tfactors.T913
                         + -6.42713 * tfactors.T9 + 0.758256 * tfactors.T953 + 236.359 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  47.2175 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 514.197 * tfactors.T943i + (1.0/3.0) * -200.896 * tfactors.T923i
                                  + -6.42713 + (5.0/3.0) * 0.758256 * tfactors.T923 + 236.359 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Na22_to_p_Mg25_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Na22 + He4 --> p + Mg25

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  44.973 + -50.0924 * tfactors.T913i + 0.807739 * tfactors.T913
                         + -0.956029 * tfactors.T9 + 0.0793321 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -50.0924 * tfactors.T943i + (1.0/3.0) * 0.807739 * tfactors.T923i
                                  + -0.956029 + (5.0/3.0) * 0.0793321 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Na23_to_He4_Ne20_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Na23 + p --> He4 + Ne20

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10r
    ln_set_rate =  -6.58736 + -2.31577 * tfactors.T9i + 19.7297 * tfactors.T913
                         + -2.20987 * tfactors.T9 + 0.153374 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  2.31577 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 19.7297 * tfactors.T923i
                                  + -2.20987 + (5.0/3.0) * 0.153374 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.86103 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -20.0024 * tfactors.T943i + (1.0/3.0) * 11.5988 * tfactors.T923i
                                  + -1.37398 + (5.0/3.0) * -1.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Al26_to_p_Si29_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Al26 + He4 --> p + Si29

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  47.7092 + -56.4422 * tfactors.T913i + 0.705353 * tfactors.T913
                         + -0.957427 * tfactors.T9 + 0.0756045 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -56.4422 * tfactors.T943i + (1.0/3.0) * 0.705353 * tfactors.T923i
                                  + -0.957427 + (5.0/3.0) * 0.0756045 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Al27_to_He4_Mg24_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Al27 + p --> He4 + Mg24

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10r
    ln_set_rate =  -7.02789 + -4.2425 * tfactors.T9i + 18.0416 * tfactors.T913
                         + -1.54137 * tfactors.T9 + 0.0847506 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  4.2425 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 18.0416 * tfactors.T923i
                                  + -1.54137 + (5.0/3.0) * 0.0847506 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.963012 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 5.18642 * tfactors.T943i + (1.0/3.0) * -34.7936 * tfactors.T923i
                                  + 168.225 + (5.0/3.0) * -115.825 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -26.4162 * tfactors.T943i
                                  + (5.0/3.0) * -2.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_P30_to_p_S33_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // P30 + He4 --> p + S33

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  49.813 + -62.4643 * tfactors.T913i + 0.492934 * tfactors.T913
                         + -0.841855 * tfactors.T9 + 0.059263 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -62.4643 * tfactors.T943i + (1.0/3.0) * 0.492934 * tfactors.T923i
                                  + -0.841855 + (5.0/3.0) * 0.059263 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_P31_to_He4_Si28_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // P31 + p --> He4 + Si28

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10r
    ln_set_rate =  -10.893 + -3.42575 * tfactors.T9i + 21.521 * tfactors.T913
                         + -1.90355 * tfactors.T9 + 0.092724 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  3.42575 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 21.521 * tfactors.T923i
                                  + -1.90355 + (5.0/3.0) * 0.092724 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.87716 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -31.932 * tfactors.T943i + (1.0/3.0) * -77.0334 * tfactors.T923i
                                  + -43.6847 + (5.0/3.0) * -4.28955 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Cl34_to_p_Ar37_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Cl34 + He4 --> p + Ar37

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  52.795 + -68.2165 * tfactors.T913i + 0.330057 * tfactors.T913
                         + -0.873334 * tfactors.T9 + 0.0592127 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -68.2165 * tfactors.T943i + (1.0/3.0) * 0.330057 * tfactors.T923i
                                  + -0.873334 + (5.0/3.0) * 0.0592127 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Cl35_to_He4_S32_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Cl35 + p --> He4 + S32

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // il10r
    ln_set_rate =  -1.01202 + -3.93495 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  3.93495 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.532931 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 25.5338 * tfactors.T923i
                                  + 6.45824 + (5.0/3.0) * -0.950294 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -30.9147 * tfactors.T943i + (1.0/3.0) * -1.2345 * tfactors.T923i
                                  + 22.5118 + (5.0/3.0) * -33.0589 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  6.00976 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 5.33756 * tfactors.T923i
                                  + 1.64418 + (5.0/3.0) * -0.246167 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_K38_to_p_Ca41_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // K38 + He4 --> p + Ca41

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  54.4204 + -73.7407 * tfactors.T913i + 0.600725 * tfactors.T913
                         + -0.957475 * tfactors.T9 + 0.0689286 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -73.7407 * tfactors.T943i + (1.0/3.0) * 0.600725 * tfactors.T923i
                                  + -0.957475 + (5.0/3.0) * 0.0689286 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_K39_to_He4_Ar36_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // K39 + p --> He4 + Ar36

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  20.5166 + -30.0732 * tfactors.T913i + 7.03263 * tfactors.T913
                         + -1.10085 * tfactors.T9 + 0.133768 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -30.0732 * tfactors.T943i + (1.0/3.0) * 7.03263 * tfactors.T923i
                                  + -1.10085 + (5.0/3.0) * 0.133768 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Sc42_to_p_Ti45_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Sc42 + He4 --> p + Ti45

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  57.6956 + -79.0687 * tfactors.T913i + -0.981235 * tfactors.T913
                         + -0.48951 * tfactors.T9 + 0.0114272 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -79.0687 * tfactors.T943i + (1.0/3.0) * -0.981235 * tfactors.T923i
                                  + -0.48951 + (5.0/3.0) * 0.0114272 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Sc43_to_He4_Ca40_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Sc43 + p --> He4 + Ca40

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  34.8559 + -32.1734 * tfactors.T913i + 0.0296879 * tfactors.T913
                         + -0.95232 * tfactors.T9 + 0.129022 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -32.1734 * tfactors.T943i + (1.0/3.0) * 0.0296879 * tfactors.T923i
                                  + -0.95232 + (5.0/3.0) * 0.129022 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Ti44_to_p_V47_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ti44 + He4 --> p + V47

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // chw0r
    ln_set_rate =  -76.5154 + -10.7931 * tfactors.T9i + 70.2835 * tfactors.T913
                         + -7.99061 * tfactors.T9 + 0.486213 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  10.7931 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 70.2835 * tfactors.T923i
                                  + -7.99061 + (5.0/3.0) * 0.486213 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_V46_to_p_Cr49_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // V46 + He4 --> p + Cr49

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  59.798 + -84.2255 * tfactors.T913i + -1.09776 * tfactors.T913
                         + -0.495216 * tfactors.T9 + 0.0124719 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -84.2255 * tfactors.T943i + (1.0/3.0) * -1.09776 * tfactors.T923i
                                  + -0.495216 + (5.0/3.0) * 0.0124719 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Cr48_to_p_Mn51_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Cr48 + He4 --> p + Mn51

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  59.2276 + -86.7459 * tfactors.T913i + 1.05653 * tfactors.T913
                         + -1.15757 * tfactors.T9 + 0.0877546 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -86.7459 * tfactors.T943i + (1.0/3.0) * 1.05653 * tfactors.T923i
                                  + -1.15757 + (5.0/3.0) * 0.0877546 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Cr49_to_p_Mn52_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Cr49 + He4 --> p + Mn52

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  52.5075 + -86.7913 * tfactors.T913i + 10.6892 * tfactors.T913
                         + -4.01586 * tfactors.T9 + 0.411275 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -86.7913 * tfactors.T943i + (1.0/3.0) * 10.6892 * tfactors.T923i
                                  + -4.01586 + (5.0/3.0) * 0.411275 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Fe52_to_p_Co55_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe52 + He4 --> p + Co55

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  62.2207 + -91.6819 * tfactors.T913i + -0.329235 * tfactors.T913
                         + -0.780924 * tfactors.T9 + 0.0425179 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -91.6819 * tfactors.T943i + (1.0/3.0) * -0.329235 * tfactors.T923i
                                  + -0.780924 + (5.0/3.0) * 0.0425179 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_He4_He4_to_C12_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // 3 He4 --> C12

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // fy05r
    ln_set_rate =  -11.7884 + -1.02446 * tfactors.T9i + -23.57 * tfactors.T913i + 20.4886 * tfactors.T913
                         + -12.9882 * tfactors.T9 + -20.0 * tfactors.T953 + -2.16667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  1.02446 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -23.57 * tfactors.T943i + (1.0/3.0) * 20.4886 * tfactors.T923i
                                  + -12.9882 + (5.0/3.0) * -20.0 * tfactors.T923 + -2.16667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -37.06 * tfactors.T943i + (1.0/3.0) * 29.3493 * tfactors.T923i
                                  + -115.507 + (5.0/3.0) * -10.0 * tfactors.T923 + -1.33333 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  4.12656 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -13.49 * tfactors.T943i + (1.0/3.0) * 21.4259 * tfactors.T923i
                                  + -1.34769 + (5.0/3.0) * 0.0879816 * tfactors.T923 + -13.1653 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_n_Mn50_to_Mn51_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Mn50 + n --> Mn51

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  15.4434 + -0.538879 * tfactors.T913
                         + 0.284528 * tfactors.T9 + -0.0459849 * tfactors.T953;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * -0.538879 * tfactors.T923i
                                  + 0.284528 + (5.0/3.0) * -0.0459849 * tfactors.T923;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_n_Mn51_to_Mn52_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Mn51 + n --> Mn52

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  17.2055 + -1.83611 * tfactors.T913
                         + 0.460384 * tfactors.T9 + -0.0584947 * tfactors.T953;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * -1.83611 * tfactors.T923i
                                  + 0.460384 + (5.0/3.0) * -0.0584947 * tfactors.T923;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Mn52_to_Fe53_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Mn52 + p --> Fe53

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  35.4319 + -36.187 * tfactors.T913i + 1.79171 * tfactors.T913
                         + -1.77786 * tfactors.T9 + 0.179936 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -36.187 * tfactors.T943i + (1.0/3.0) * 1.79171 * tfactors.T923i
                                  + -1.77786 + (5.0/3.0) * 0.179936 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Mn52_to_Co56_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Mn52 + He4 --> Co56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  63.5657 + -89.3157 * tfactors.T913i + -1.95982 * tfactors.T913
                         + -2.45671 * tfactors.T9 + 0.296095 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -89.3157 * tfactors.T943i + (1.0/3.0) * -1.95982 * tfactors.T923i
                                  + -2.45671 + (5.0/3.0) * 0.296095 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_n_Fe52_to_Fe53_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe52 + n --> Fe53

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  15.8885 + -0.344319 * tfactors.T913
                         + 0.178277 * tfactors.T9 + -0.0334326 * tfactors.T953;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * -0.344319 * tfactors.T923i
                                  + 0.178277 + (5.0/3.0) * -0.0334326 * tfactors.T923;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_n_Fe53_to_Fe54_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe53 + n --> Fe54

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  16.4534 + -1.10421 * tfactors.T913
                         + 0.379905 * tfactors.T9 + -0.0581878 * tfactors.T953;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * -1.10421 * tfactors.T923i
                                  + 0.379905 + (5.0/3.0) * -0.0581878 * tfactors.T923;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Fe53_to_Ni57_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe53 + He4 --> Ni57

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  69.0636 + -91.7231 * tfactors.T913i + -10.1816 * tfactors.T913
                         + -0.0406303 * tfactors.T9 + 0.0345056 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -91.7231 * tfactors.T943i + (1.0/3.0) * -10.1816 * tfactors.T923i
                                  + -0.0406303 + (5.0/3.0) * 0.0345056 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_n_Fe54_to_Fe55_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe54 + n --> Fe55

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ks03 
    ln_set_rate =  -0.80864 + 0.0591716 * tfactors.T9i + -8.66617 * tfactors.T913i + 26.4472 * tfactors.T913
                         + -1.9222 * tfactors.T9 + 0.0986404 * tfactors.T953 + -9.78317 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  -0.0591716 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -8.66617 * tfactors.T943i + (1.0/3.0) * 26.4472 * tfactors.T923i
                                  + -1.9222 + (5.0/3.0) * 0.0986404 * tfactors.T923 + -9.78317 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Fe54_to_Co55_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe54 + p --> Co55

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  36.2304 + -37.1544 * tfactors.T913i + 0.950364 * tfactors.T913
                         + -1.77529 * tfactors.T9 + 0.198562 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -37.1544 * tfactors.T943i + (1.0/3.0) * 0.950364 * tfactors.T923i
                                  + -1.77529 + (5.0/3.0) * 0.198562 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Fe54_to_Ni58_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe54 + He4 --> Ni58

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  60.2478 + -91.7628 * tfactors.T913i + 4.23027 * tfactors.T913
                         + -3.31305 * tfactors.T9 + 0.271293 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -91.7628 * tfactors.T943i + (1.0/3.0) * 4.23027 * tfactors.T923i
                                  + -3.31305 + (5.0/3.0) * 0.271293 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_n_Fe55_to_Fe56_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe55 + n --> Fe56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ks03 
    ln_set_rate =  21.7202 + -0.0955677 * tfactors.T9i + 8.06062 * tfactors.T913i + -14.4809 * tfactors.T913
                         + 0.94252 * tfactors.T9 + -0.0776007 * tfactors.T953 + 6.47093 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  0.0955677 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 8.06062 * tfactors.T943i + (1.0/3.0) * -14.4809 * tfactors.T923i
                                  + 0.94252 + (5.0/3.0) * -0.0776007 * tfactors.T923 + 6.47093 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Fe55_to_Co56_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe55 + p --> Co56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  35.7333 + -37.1585 * tfactors.T913i + 1.66198 * tfactors.T913
                         + -1.60842 * tfactors.T9 + 0.148916 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -37.1585 * tfactors.T943i + (1.0/3.0) * 1.66198 * tfactors.T923i
                                  + -1.60842 + (5.0/3.0) * 0.148916 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Fe56_to_Co57_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe56 + p --> Co57

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  36.0665 + -37.1625 * tfactors.T913i + 1.06776 * tfactors.T913
                         + -1.31689 * tfactors.T9 + 0.122089 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -37.1625 * tfactors.T943i + (1.0/3.0) * 1.06776 * tfactors.T923i
                                  + -1.31689 + (5.0/3.0) * 0.122089 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_n_Co55_to_Co56_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co55 + n --> Co56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  17.065 + -1.86357 * tfactors.T913
                         + 0.616591 * tfactors.T9 + -0.0839313 * tfactors.T953;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * -1.86357 * tfactors.T923i
                                  + 0.616591 + (5.0/3.0) * -0.0839313 * tfactors.T923;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Co55_to_Cu59_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co55 + He4 --> Cu59

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  65.2921 + -94.1404 * tfactors.T913i + -2.62786 * tfactors.T913
                         + -2.12066 * tfactors.T9 + 0.237999 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -94.1404 * tfactors.T943i + (1.0/3.0) * -2.62786 * tfactors.T923i
                                  + -2.12066 + (5.0/3.0) * 0.237999 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_n_Co56_to_Co57_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co56 + n --> Co57

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  17.3552 + -1.37855 * tfactors.T913
                         + 0.299896 * tfactors.T9 + -0.04382 * tfactors.T953;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * -1.37855 * tfactors.T923i
                                  + 0.299896 + (5.0/3.0) * -0.04382 * tfactors.T923;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Co56_to_Ni57_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co56 + p --> Ni57

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  35.5937 + -38.1094 * tfactors.T913i + 2.57091 * tfactors.T913
                         + -2.07795 * tfactors.T9 + 0.20757 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -38.1094 * tfactors.T943i + (1.0/3.0) * 2.57091 * tfactors.T923i
                                  + -2.07795 + (5.0/3.0) * 0.20757 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Co57_to_Ni58_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co57 + p --> Ni58

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  36.0159 + -38.1133 * tfactors.T913i + 1.77414 * tfactors.T913
                         + -1.48268 * tfactors.T9 + 0.121073 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -38.1133 * tfactors.T943i + (1.0/3.0) * 1.77414 * tfactors.T923i
                                  + -1.48268 + (5.0/3.0) * 0.121073 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_n_Ni56_to_Ni57_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni56 + n --> Ni57

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  16.0765 + -1.19665 * tfactors.T913
                         + 0.507179 * tfactors.T9 + -0.074604 * tfactors.T953;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * -1.19665 * tfactors.T923i
                                  + 0.507179 + (5.0/3.0) * -0.074604 * tfactors.T923;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Ni56_to_Zn60_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni56 + He4 --> Zn60

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  61.0733 + -96.4898 * tfactors.T913i + 6.47209 * tfactors.T913
                         + -5.2029 * tfactors.T9 + 0.533391 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -96.4898 * tfactors.T943i + (1.0/3.0) * 6.47209 * tfactors.T923i
                                  + -5.2029 + (5.0/3.0) * 0.533391 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_n_Ni57_to_Ni58_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni57 + n --> Ni58

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  17.2731 + -1.90814 * tfactors.T913
                         + 0.493188 * tfactors.T9 + -0.0684633 * tfactors.T953;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * -1.90814 * tfactors.T923i
                                  + 0.493188 + (5.0/3.0) * -0.0684633 * tfactors.T923;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Ni58_to_Cu59_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni58 + p --> Cu59

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // si13n
    ln_set_rate =  36.6662 + -39.0526 * tfactors.T913i + 1.0436 * tfactors.T913
                         + -2.10834 * tfactors.T9 + 0.239559 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -39.0526 * tfactors.T943i + (1.0/3.0) * 1.0436 * tfactors.T923i
                                  + -2.10834 + (5.0/3.0) * 0.239559 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Cu59_to_Zn60_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Cu59 + p --> Zn60

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  37.226 + -39.9808 * tfactors.T913i + 1.17311 * tfactors.T913
                         + -2.90486 * tfactors.T9 + 0.339644 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -39.9808 * tfactors.T943i + (1.0/3.0) * 1.17311 * tfactors.T923i
                                  + -2.90486 + (5.0/3.0) * 0.339644 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Mn50_to_p_Fe53_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Mn50 + He4 --> p + Fe53

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  62.6721 + -89.2309 * tfactors.T913i + -1.4758 * tfactors.T913
                         + -0.437831 * tfactors.T9 + -0.00354594 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -89.2309 * tfactors.T943i + (1.0/3.0) * -1.4758 * tfactors.T923i
                                  + -0.437831 + (5.0/3.0) * -0.00354594 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Mn51_to_p_Fe54_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Mn51 + He4 --> p + Fe54

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  62.2777 + -89.274 * tfactors.T913i + -0.862452 * tfactors.T913
                         + -0.635672 * tfactors.T9 + 0.0196464 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -89.274 * tfactors.T943i + (1.0/3.0) * -0.862452 * tfactors.T923i
                                  + -0.635672 + (5.0/3.0) * 0.0196464 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Mn52_to_p_Fe55_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Mn52 + He4 --> p + Fe55

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  56.6689 + -89.3157 * tfactors.T913i + 5.65639 * tfactors.T913
                         + -2.28952 * tfactors.T9 + 0.186017 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -89.3157 * tfactors.T943i + (1.0/3.0) * 5.65639 * tfactors.T923i
                                  + -2.28952 + (5.0/3.0) * 0.186017 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_n_Fe52_to_p_Mn52_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe52 + n --> p + Mn52

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  17.9485 + 0.569856 * tfactors.T913
                         + 0.390117 * tfactors.T9 + -0.0479043 * tfactors.T953;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * 0.569856 * tfactors.T923i
                                  + 0.390117 + (5.0/3.0) * -0.0479043 * tfactors.T923;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Fe53_to_p_Co56_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Fe53 + He4 --> p + Co56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  54.8614 + -91.7231 * tfactors.T913i + 10.052 * tfactors.T913
                         + -3.86332 * tfactors.T9 + 0.391105 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -91.7231 * tfactors.T943i + (1.0/3.0) * 10.052 * tfactors.T923i
                                  + -3.86332 + (5.0/3.0) * 0.391105 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_n_Co55_to_p_Fe55_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co55 + n --> p + Fe55

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  20.7398 + -1.62382 * tfactors.T913
                         + 0.58115 * tfactors.T9 + -0.0537057 * tfactors.T953;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * -1.62382 * tfactors.T923i
                                  + 0.58115 + (5.0/3.0) * -0.0537057 * tfactors.T923;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_n_Co55_to_He4_Mn52_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co55 + n --> He4 + Mn52

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  -12.5469 + 8.15876 * tfactors.T913
                         + 3.05978 * tfactors.T9 + -0.406535 * tfactors.T953;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * 8.15876 * tfactors.T923i
                                  + 3.05978 + (5.0/3.0) * -0.406535 * tfactors.T923;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Co55_to_p_Ni58_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co55 + He4 --> p + Ni58

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  60.2281 + -94.1404 * tfactors.T913i + 3.39179 * tfactors.T913
                         + -1.71062 * tfactors.T9 + 0.133003 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -94.1404 * tfactors.T943i + (1.0/3.0) * 3.39179 * tfactors.T923i
                                  + -1.71062 + (5.0/3.0) * 0.133003 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_n_Co56_to_p_Fe56_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co56 + n --> p + Fe56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  20.4539 + -1.13331 * tfactors.T913
                         + 0.347185 * tfactors.T9 + -0.0328879 * tfactors.T953;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * -1.13331 * tfactors.T923i
                                  + 0.347185 + (5.0/3.0) * -0.0328879 * tfactors.T923;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Co57_to_He4_Fe54_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Co57 + p --> He4 + Fe54

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  -2.1962 + -38.1133 * tfactors.T913i + 29.3541 * tfactors.T913
                         + -4.75966 * tfactors.T9 + 0.40418 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -38.1133 * tfactors.T943i + (1.0/3.0) * 29.3541 * tfactors.T923i
                                  + -4.75966 + (5.0/3.0) * 0.40418 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_n_Ni56_to_p_Co56_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni56 + n --> p + Co56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  15.5693 + 1.76846 * tfactors.T913
                         + 0.197992 * tfactors.T9 + -0.017494 * tfactors.T953;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * 1.76846 * tfactors.T923i
                                  + 0.197992 + (5.0/3.0) * -0.017494 * tfactors.T923;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_n_Ni56_to_He4_Fe53_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni56 + n --> He4 + Fe53

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  -7.84393 + 6.43259 * tfactors.T913
                         + 2.29591 * tfactors.T9 + -0.254724 * tfactors.T953;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * 6.43259 * tfactors.T923i
                                  + 2.29591 + (5.0/3.0) * -0.254724 * tfactors.T923;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_n_Ni57_to_p_Co57_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni57 + n --> p + Co57

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  19.9939 + -1.13177 * tfactors.T913
                         + 0.56875 * tfactors.T9 + -0.0579913 * tfactors.T953;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * -1.13177 * tfactors.T923i
                                  + 0.56875 + (5.0/3.0) * -0.0579913 * tfactors.T923;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_n_Ni57_to_He4_Fe54_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni57 + n --> He4 + Fe54

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  14.5327 + 0.48562 * tfactors.T913
                         + -0.0249302 * tfactors.T9 + 0.0509605 * tfactors.T953;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * 0.48562 * tfactors.T923i
                                  + -0.0249302 + (5.0/3.0) * 0.0509605 * tfactors.T923;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_n_Ni58_to_He4_Fe55_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Ni58 + n --> He4 + Fe55

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  -2.03901 + 3.44996 * tfactors.T913
                         + 2.98226 * tfactors.T9 + -0.387699 * tfactors.T953;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * 3.44996 * tfactors.T923i
                                  + 2.98226 + (5.0/3.0) * -0.387699 * tfactors.T923;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_n_Cu59_to_He4_Co56_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Cu59 + n --> He4 + Co56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  11.0108 + -1.00882 * tfactors.T913
                         + 1.88993 * tfactors.T9 + -0.202319 * tfactors.T953;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * -1.00882 * tfactors.T923i
                                  + 1.88993 + (5.0/3.0) * -0.202319 * tfactors.T923;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Cu59_to_He4_Ni56_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Cu59 + p --> He4 + Ni56

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  20.7688 + -39.9808 * tfactors.T913i + 13.6773 * tfactors.T913
                         + -3.76429 * tfactors.T9 + 0.438096 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + -(1.0/3.0) * -39.9808 * tfactors.T943i + (1.0/3.0) * 13.6773 * tfactors.T923i
                                  + -3.76429 + (5.0/3.0) * 0.438096 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
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
AMREX_GPU_DEVICE AMREX_INLINE
void rate_n_Zn60_to_He4_Ni57_reaclib(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

    // Zn60 + n --> He4 + Ni57

    rate = 0.0;
    drate_dT = 0.0;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ths8r
    ln_set_rate =  12.4241 + 5.43552 * tfactors.T913
                         + -1.25169 * tfactors.T9 + 0.148383 * tfactors.T953;

    ln_set_rate += log_scor;

    if constexpr (do_T_derivatives) {
        dln_set_rate_dT9 =  + (1.0/3.0) * 5.43552 * tfactors.T923i
                                  + -1.25169 + (5.0/3.0) * 0.148383 * tfactors.T923;

        dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (do_T_derivatives) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
    }

}


template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
fill_reaclib_rates([[maybe_unused]] const tf_t& tfactors,
                   [[maybe_unused]] T& rate_eval)
{

    [[maybe_unused]] amrex::Real rate;
    [[maybe_unused]] amrex::Real drate_dT;

    {
        // C14_to_N14_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_C14_to_N14_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_C14_to_N14_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_C14_to_N14_reaclib) = drate_dT;
        }
    }

    {
        // N13_to_C13_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_N13_to_C13_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_N13_to_C13_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_N13_to_C13_reaclib) = drate_dT;
        }
    }

    {
        // p_C12_to_N13_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_C12);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_C12);
        }
#endif
        rate_p_C12_to_N13_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_C12_to_N13_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_C12_to_N13_reaclib) = drate_dT;
        }
    }

    {
        // He4_C12_to_O16_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_C12);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_C12);
        }
#endif
        rate_He4_C12_to_O16_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_C12_to_O16_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_C12_to_O16_reaclib) = drate_dT;
        }
    }

    {
        // p_C13_to_N14_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_C13);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_C13);
        }
#endif
        rate_p_C13_to_N14_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_C13_to_N14_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_C13_to_N14_reaclib) = drate_dT;
        }
    }

    {
        // p_C14_to_N15_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_C14);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_C14);
        }
#endif
        rate_p_C14_to_N15_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_C14_to_N15_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_C14_to_N15_reaclib) = drate_dT;
        }
    }

    {
        // He4_C14_to_O18_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_C14);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_C14);
        }
#endif
        rate_He4_C14_to_O18_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_C14_to_O18_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_C14_to_O18_reaclib) = drate_dT;
        }
    }

    {
        // He4_N14_to_F18_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_N14);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_N14);
        }
#endif
        rate_He4_N14_to_F18_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_N14_to_F18_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_N14_to_F18_reaclib) = drate_dT;
        }
    }

    {
        // p_N15_to_O16_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_N15);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_N15);
        }
#endif
        rate_p_N15_to_O16_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_N15_to_O16_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_N15_to_O16_reaclib) = drate_dT;
        }
    }

    {
        // He4_N15_to_F19_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_N15);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_N15);
        }
#endif
        rate_He4_N15_to_F19_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_N15_to_F19_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_N15_to_F19_reaclib) = drate_dT;
        }
    }

    {
        // He4_O16_to_Ne20_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_O16);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_O16);
        }
#endif
        rate_He4_O16_to_Ne20_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_O16_to_Ne20_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_O16_to_Ne20_reaclib) = drate_dT;
        }
    }

    {
        // p_O17_to_F18_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_O17);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_O17);
        }
#endif
        rate_p_O17_to_F18_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_O17_to_F18_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_O17_to_F18_reaclib) = drate_dT;
        }
    }

    {
        // He4_O17_to_Ne21_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_O17);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_O17);
        }
#endif
        rate_He4_O17_to_Ne21_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_O17_to_Ne21_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_O17_to_Ne21_reaclib) = drate_dT;
        }
    }

    {
        // p_O18_to_F19_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_O18);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_O18);
        }
#endif
        rate_p_O18_to_F19_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_O18_to_F19_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_O18_to_F19_reaclib) = drate_dT;
        }
    }

    {
        // He4_O18_to_Ne22_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_O18);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_O18);
        }
#endif
        rate_He4_O18_to_Ne22_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_O18_to_Ne22_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_O18_to_Ne22_reaclib) = drate_dT;
        }
    }

    {
        // He4_F18_to_Na22_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_F18);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_F18);
        }
#endif
        rate_He4_F18_to_Na22_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_F18_to_Na22_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_F18_to_Na22_reaclib) = drate_dT;
        }
    }

    {
        // p_F19_to_Ne20_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_F19);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_F19);
        }
#endif
        rate_p_F19_to_Ne20_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_F19_to_Ne20_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_F19_to_Ne20_reaclib) = drate_dT;
        }
    }

    {
        // He4_F19_to_Na23_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_F19);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_F19);
        }
#endif
        rate_He4_F19_to_Na23_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_F19_to_Na23_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_F19_to_Na23_reaclib) = drate_dT;
        }
    }

    {
        // He4_Ne20_to_Mg24_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Ne20);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ne20);
        }
#endif
        rate_He4_Ne20_to_Mg24_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Ne20_to_Mg24_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Ne20_to_Mg24_reaclib) = drate_dT;
        }
    }

    {
        // p_Ne21_to_Na22_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Ne21);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ne21);
        }
#endif
        rate_p_Ne21_to_Na22_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Ne21_to_Na22_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Ne21_to_Na22_reaclib) = drate_dT;
        }
    }

    {
        // He4_Ne21_to_Mg25_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Ne21);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ne21);
        }
#endif
        rate_He4_Ne21_to_Mg25_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Ne21_to_Mg25_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Ne21_to_Mg25_reaclib) = drate_dT;
        }
    }

    {
        // p_Ne22_to_Na23_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Ne22);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ne22);
        }
#endif
        rate_p_Ne22_to_Na23_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Ne22_to_Na23_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Ne22_to_Na23_reaclib) = drate_dT;
        }
    }

    {
        // He4_Na22_to_Al26_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Na22);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Na22);
        }
#endif
        rate_He4_Na22_to_Al26_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Na22_to_Al26_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Na22_to_Al26_reaclib) = drate_dT;
        }
    }

    {
        // p_Na23_to_Mg24_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Na23);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Na23);
        }
#endif
        rate_p_Na23_to_Mg24_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Na23_to_Mg24_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Na23_to_Mg24_reaclib) = drate_dT;
        }
    }

    {
        // He4_Na23_to_Al27_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Na23);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Na23);
        }
#endif
        rate_He4_Na23_to_Al27_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Na23_to_Al27_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Na23_to_Al27_reaclib) = drate_dT;
        }
    }

    {
        // He4_Mg24_to_Si28_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Mg24);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mg24);
        }
#endif
        rate_He4_Mg24_to_Si28_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Mg24_to_Si28_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Mg24_to_Si28_reaclib) = drate_dT;
        }
    }

    {
        // p_Mg25_to_Al26_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Mg25);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mg25);
        }
#endif
        rate_p_Mg25_to_Al26_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Mg25_to_Al26_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Mg25_to_Al26_reaclib) = drate_dT;
        }
    }

    {
        // He4_Mg25_to_Si29_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Mg25);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mg25);
        }
#endif
        rate_He4_Mg25_to_Si29_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Mg25_to_Si29_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Mg25_to_Si29_reaclib) = drate_dT;
        }
    }

    {
        // He4_Al26_to_P30_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Al26);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Al26);
        }
#endif
        rate_He4_Al26_to_P30_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Al26_to_P30_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Al26_to_P30_reaclib) = drate_dT;
        }
    }

    {
        // p_Al27_to_Si28_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Al27);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Al27);
        }
#endif
        rate_p_Al27_to_Si28_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Al27_to_Si28_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Al27_to_Si28_reaclib) = drate_dT;
        }
    }

    {
        // He4_Al27_to_P31_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Al27);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Al27);
        }
#endif
        rate_He4_Al27_to_P31_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Al27_to_P31_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Al27_to_P31_reaclib) = drate_dT;
        }
    }

    {
        // He4_Si28_to_S32_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Si28);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si28);
        }
#endif
        rate_He4_Si28_to_S32_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Si28_to_S32_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Si28_to_S32_reaclib) = drate_dT;
        }
    }

    {
        // p_Si29_to_P30_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Si29);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Si29);
        }
#endif
        rate_p_Si29_to_P30_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Si29_to_P30_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Si29_to_P30_reaclib) = drate_dT;
        }
    }

    {
        // He4_Si29_to_S33_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Si29);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si29);
        }
#endif
        rate_He4_Si29_to_S33_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Si29_to_S33_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Si29_to_S33_reaclib) = drate_dT;
        }
    }

    {
        // He4_P30_to_Cl34_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_P30);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_P30);
        }
#endif
        rate_He4_P30_to_Cl34_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_P30_to_Cl34_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_P30_to_Cl34_reaclib) = drate_dT;
        }
    }

    {
        // p_P31_to_S32_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_P31);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_P31);
        }
#endif
        rate_p_P31_to_S32_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_P31_to_S32_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_P31_to_S32_reaclib) = drate_dT;
        }
    }

    {
        // He4_P31_to_Cl35_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_P31);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_P31);
        }
#endif
        rate_He4_P31_to_Cl35_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_P31_to_Cl35_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_P31_to_Cl35_reaclib) = drate_dT;
        }
    }

    {
        // He4_S32_to_Ar36_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_S32);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_S32);
        }
#endif
        rate_He4_S32_to_Ar36_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_S32_to_Ar36_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_S32_to_Ar36_reaclib) = drate_dT;
        }
    }

    {
        // p_S33_to_Cl34_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_S33);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_S33);
        }
#endif
        rate_p_S33_to_Cl34_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_S33_to_Cl34_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_S33_to_Cl34_reaclib) = drate_dT;
        }
    }

    {
        // He4_S33_to_Ar37_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_S33);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_S33);
        }
#endif
        rate_He4_S33_to_Ar37_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_S33_to_Ar37_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_S33_to_Ar37_reaclib) = drate_dT;
        }
    }

    {
        // He4_Cl34_to_K38_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Cl34);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cl34);
        }
#endif
        rate_He4_Cl34_to_K38_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Cl34_to_K38_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Cl34_to_K38_reaclib) = drate_dT;
        }
    }

    {
        // p_Cl35_to_Ar36_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Cl35);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cl35);
        }
#endif
        rate_p_Cl35_to_Ar36_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Cl35_to_Ar36_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Cl35_to_Ar36_reaclib) = drate_dT;
        }
    }

    {
        // He4_Cl35_to_K39_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Cl35);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cl35);
        }
#endif
        rate_He4_Cl35_to_K39_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Cl35_to_K39_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Cl35_to_K39_reaclib) = drate_dT;
        }
    }

    {
        // He4_Ar36_to_Ca40_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Ar36);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ar36);
        }
#endif
        rate_He4_Ar36_to_Ca40_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Ar36_to_Ca40_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Ar36_to_Ca40_reaclib) = drate_dT;
        }
    }

    {
        // p_Ar37_to_K38_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Ar37);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ar37);
        }
#endif
        rate_p_Ar37_to_K38_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Ar37_to_K38_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Ar37_to_K38_reaclib) = drate_dT;
        }
    }

    {
        // He4_Ar37_to_Ca41_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Ar37);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ar37);
        }
#endif
        rate_He4_Ar37_to_Ca41_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Ar37_to_Ca41_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Ar37_to_Ca41_reaclib) = drate_dT;
        }
    }

    {
        // He4_K38_to_Sc42_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_K38);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_K38);
        }
#endif
        rate_He4_K38_to_Sc42_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_K38_to_Sc42_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_K38_to_Sc42_reaclib) = drate_dT;
        }
    }

    {
        // p_K39_to_Ca40_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_K39);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_K39);
        }
#endif
        rate_p_K39_to_Ca40_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_K39_to_Ca40_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_K39_to_Ca40_reaclib) = drate_dT;
        }
    }

    {
        // He4_K39_to_Sc43_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_K39);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_K39);
        }
#endif
        rate_He4_K39_to_Sc43_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_K39_to_Sc43_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_K39_to_Sc43_reaclib) = drate_dT;
        }
    }

    {
        // He4_Ca40_to_Ti44_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Ca40);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca40);
        }
#endif
        rate_He4_Ca40_to_Ti44_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Ca40_to_Ti44_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Ca40_to_Ti44_reaclib) = drate_dT;
        }
    }

    {
        // p_Ca41_to_Sc42_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Ca41);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca41);
        }
#endif
        rate_p_Ca41_to_Sc42_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Ca41_to_Sc42_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Ca41_to_Sc42_reaclib) = drate_dT;
        }
    }

    {
        // He4_Ca41_to_Ti45_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Ca41);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca41);
        }
#endif
        rate_He4_Ca41_to_Ti45_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Ca41_to_Ti45_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Ca41_to_Ti45_reaclib) = drate_dT;
        }
    }

    {
        // He4_Sc42_to_V46_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Sc42);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc42);
        }
#endif
        rate_He4_Sc42_to_V46_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Sc42_to_V46_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Sc42_to_V46_reaclib) = drate_dT;
        }
    }

    {
        // p_Sc43_to_Ti44_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Sc43);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc43);
        }
#endif
        rate_p_Sc43_to_Ti44_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Sc43_to_Ti44_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Sc43_to_Ti44_reaclib) = drate_dT;
        }
    }

    {
        // He4_Sc43_to_V47_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Sc43);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc43);
        }
#endif
        rate_He4_Sc43_to_V47_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Sc43_to_V47_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Sc43_to_V47_reaclib) = drate_dT;
        }
    }

    {
        // He4_Ti44_to_Cr48_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Ti44);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti44);
        }
#endif
        rate_He4_Ti44_to_Cr48_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Ti44_to_Cr48_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Ti44_to_Cr48_reaclib) = drate_dT;
        }
    }

    {
        // p_Ti45_to_V46_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Ti45);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti45);
        }
#endif
        rate_p_Ti45_to_V46_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Ti45_to_V46_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Ti45_to_V46_reaclib) = drate_dT;
        }
    }

    {
        // He4_Ti45_to_Cr49_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Ti45);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti45);
        }
#endif
        rate_He4_Ti45_to_Cr49_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Ti45_to_Cr49_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Ti45_to_Cr49_reaclib) = drate_dT;
        }
    }

    {
        // He4_V46_to_Mn50_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_V46);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V46);
        }
#endif
        rate_He4_V46_to_Mn50_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_V46_to_Mn50_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_V46_to_Mn50_reaclib) = drate_dT;
        }
    }

    {
        // p_V47_to_Cr48_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_V47);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V47);
        }
#endif
        rate_p_V47_to_Cr48_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_V47_to_Cr48_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_V47_to_Cr48_reaclib) = drate_dT;
        }
    }

    {
        // He4_V47_to_Mn51_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_V47);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V47);
        }
#endif
        rate_He4_V47_to_Mn51_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_V47_to_Mn51_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_V47_to_Mn51_reaclib) = drate_dT;
        }
    }

    {
        // He4_Cr48_to_Fe52_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Cr48);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr48);
        }
#endif
        rate_He4_Cr48_to_Fe52_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Cr48_to_Fe52_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Cr48_to_Fe52_reaclib) = drate_dT;
        }
    }

    {
        // p_Cr49_to_Mn50_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Cr49);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr49);
        }
#endif
        rate_p_Cr49_to_Mn50_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Cr49_to_Mn50_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Cr49_to_Mn50_reaclib) = drate_dT;
        }
    }

    {
        // p_Mn51_to_Fe52_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Mn51);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn51);
        }
#endif
        rate_p_Mn51_to_Fe52_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Mn51_to_Fe52_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Mn51_to_Fe52_reaclib) = drate_dT;
        }
    }

    {
        // He4_Mn51_to_Co55_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Mn51);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn51);
        }
#endif
        rate_He4_Mn51_to_Co55_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Mn51_to_Co55_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Mn51_to_Co55_reaclib) = drate_dT;
        }
    }

    {
        // He4_Fe52_to_Ni56_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Fe52);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe52);
        }
#endif
        rate_He4_Fe52_to_Ni56_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Fe52_to_Ni56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Fe52_to_Ni56_reaclib) = drate_dT;
        }
    }

    {
        // p_Co55_to_Ni56_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Co55);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co55);
        }
#endif
        rate_p_Co55_to_Ni56_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Co55_to_Ni56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Co55_to_Ni56_reaclib) = drate_dT;
        }
    }

    {
        // C12_C12_to_p_Na23_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_C12_C12);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_C12);
        }
#endif
        rate_C12_C12_to_p_Na23_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_C12_C12_to_p_Na23_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_C12_C12_to_p_Na23_reaclib) = drate_dT;
        }
    }

    {
        // C12_C12_to_He4_Ne20_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_C12_C12);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_C12);
        }
#endif
        rate_C12_C12_to_He4_Ne20_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_C12_C12_to_He4_Ne20_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_C12_C12_to_He4_Ne20_reaclib) = drate_dT;
        }
    }

    {
        // He4_N13_to_p_O16_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_N13);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_N13);
        }
#endif
        rate_He4_N13_to_p_O16_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_N13_to_p_O16_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_N13_to_p_O16_reaclib) = drate_dT;
        }
    }

    {
        // p_N15_to_He4_C12_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_N15);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_N15);
        }
#endif
        rate_p_N15_to_He4_C12_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_N15_to_He4_C12_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_N15_to_He4_C12_reaclib) = drate_dT;
        }
    }

    {
        // C12_O16_to_p_Al27_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_C12_O16);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_O16);
        }
#endif
        rate_C12_O16_to_p_Al27_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_C12_O16_to_p_Al27_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_C12_O16_to_p_Al27_reaclib) = drate_dT;
        }
    }

    {
        // C12_O16_to_He4_Mg24_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_C12_O16);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_O16);
        }
#endif
        rate_C12_O16_to_He4_Mg24_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_C12_O16_to_He4_Mg24_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_C12_O16_to_He4_Mg24_reaclib) = drate_dT;
        }
    }

    {
        // O16_O16_to_p_P31_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_O16_O16);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_O16_O16);
        }
#endif
        rate_O16_O16_to_p_P31_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_O16_O16_to_p_P31_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_O16_O16_to_p_P31_reaclib) = drate_dT;
        }
    }

    {
        // O16_O16_to_He4_Si28_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_O16_O16);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_O16_O16);
        }
#endif
        rate_O16_O16_to_He4_Si28_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_O16_O16_to_He4_Si28_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_O16_O16_to_He4_Si28_reaclib) = drate_dT;
        }
    }

    {
        // p_O17_to_He4_N14_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_O17);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_O17);
        }
#endif
        rate_p_O17_to_He4_N14_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_O17_to_He4_N14_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_O17_to_He4_N14_reaclib) = drate_dT;
        }
    }

    {
        // p_O18_to_He4_N15_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_O18);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_O18);
        }
#endif
        rate_p_O18_to_He4_N15_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_O18_to_He4_N15_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_O18_to_He4_N15_reaclib) = drate_dT;
        }
    }

    {
        // He4_F18_to_p_Ne21_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_F18);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_F18);
        }
#endif
        rate_He4_F18_to_p_Ne21_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_F18_to_p_Ne21_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_F18_to_p_Ne21_reaclib) = drate_dT;
        }
    }

    {
        // p_F19_to_He4_O16_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_F19);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_F19);
        }
#endif
        rate_p_F19_to_He4_O16_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_F19_to_He4_O16_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_F19_to_He4_O16_reaclib) = drate_dT;
        }
    }

    {
        // He4_F19_to_p_Ne22_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_F19);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_F19);
        }
#endif
        rate_He4_F19_to_p_Ne22_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_F19_to_p_Ne22_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_F19_to_p_Ne22_reaclib) = drate_dT;
        }
    }

    {
        // C12_Ne20_to_p_P31_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_C12_Ne20);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_Ne20);
        }
#endif
        rate_C12_Ne20_to_p_P31_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_C12_Ne20_to_p_P31_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_C12_Ne20_to_p_P31_reaclib) = drate_dT;
        }
    }

    {
        // C12_Ne20_to_He4_Si28_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_C12_Ne20);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_Ne20);
        }
#endif
        rate_C12_Ne20_to_He4_Si28_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_C12_Ne20_to_He4_Si28_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_C12_Ne20_to_He4_Si28_reaclib) = drate_dT;
        }
    }

    {
        // He4_Na22_to_p_Mg25_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Na22);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Na22);
        }
#endif
        rate_He4_Na22_to_p_Mg25_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Na22_to_p_Mg25_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Na22_to_p_Mg25_reaclib) = drate_dT;
        }
    }

    {
        // p_Na23_to_He4_Ne20_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Na23);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Na23);
        }
#endif
        rate_p_Na23_to_He4_Ne20_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Na23_to_He4_Ne20_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Na23_to_He4_Ne20_reaclib) = drate_dT;
        }
    }

    {
        // He4_Al26_to_p_Si29_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Al26);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Al26);
        }
#endif
        rate_He4_Al26_to_p_Si29_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Al26_to_p_Si29_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Al26_to_p_Si29_reaclib) = drate_dT;
        }
    }

    {
        // p_Al27_to_He4_Mg24_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Al27);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Al27);
        }
#endif
        rate_p_Al27_to_He4_Mg24_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Al27_to_He4_Mg24_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Al27_to_He4_Mg24_reaclib) = drate_dT;
        }
    }

    {
        // He4_P30_to_p_S33_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_P30);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_P30);
        }
#endif
        rate_He4_P30_to_p_S33_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_P30_to_p_S33_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_P30_to_p_S33_reaclib) = drate_dT;
        }
    }

    {
        // p_P31_to_He4_Si28_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_P31);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_P31);
        }
#endif
        rate_p_P31_to_He4_Si28_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_P31_to_He4_Si28_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_P31_to_He4_Si28_reaclib) = drate_dT;
        }
    }

    {
        // He4_Cl34_to_p_Ar37_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Cl34);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cl34);
        }
#endif
        rate_He4_Cl34_to_p_Ar37_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Cl34_to_p_Ar37_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Cl34_to_p_Ar37_reaclib) = drate_dT;
        }
    }

    {
        // p_Cl35_to_He4_S32_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Cl35);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cl35);
        }
#endif
        rate_p_Cl35_to_He4_S32_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Cl35_to_He4_S32_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Cl35_to_He4_S32_reaclib) = drate_dT;
        }
    }

    {
        // He4_K38_to_p_Ca41_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_K38);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_K38);
        }
#endif
        rate_He4_K38_to_p_Ca41_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_K38_to_p_Ca41_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_K38_to_p_Ca41_reaclib) = drate_dT;
        }
    }

    {
        // p_K39_to_He4_Ar36_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_K39);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_K39);
        }
#endif
        rate_p_K39_to_He4_Ar36_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_K39_to_He4_Ar36_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_K39_to_He4_Ar36_reaclib) = drate_dT;
        }
    }

    {
        // He4_Sc42_to_p_Ti45_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Sc42);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc42);
        }
#endif
        rate_He4_Sc42_to_p_Ti45_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Sc42_to_p_Ti45_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Sc42_to_p_Ti45_reaclib) = drate_dT;
        }
    }

    {
        // p_Sc43_to_He4_Ca40_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Sc43);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc43);
        }
#endif
        rate_p_Sc43_to_He4_Ca40_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Sc43_to_He4_Ca40_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Sc43_to_He4_Ca40_reaclib) = drate_dT;
        }
    }

    {
        // He4_Ti44_to_p_V47_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Ti44);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti44);
        }
#endif
        rate_He4_Ti44_to_p_V47_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Ti44_to_p_V47_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Ti44_to_p_V47_reaclib) = drate_dT;
        }
    }

    {
        // He4_V46_to_p_Cr49_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_V46);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V46);
        }
#endif
        rate_He4_V46_to_p_Cr49_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_V46_to_p_Cr49_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_V46_to_p_Cr49_reaclib) = drate_dT;
        }
    }

    {
        // He4_Cr48_to_p_Mn51_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Cr48);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr48);
        }
#endif
        rate_He4_Cr48_to_p_Mn51_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Cr48_to_p_Mn51_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Cr48_to_p_Mn51_reaclib) = drate_dT;
        }
    }

    {
        // He4_Cr49_to_p_Mn52_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Cr49);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr49);
        }
#endif
        rate_He4_Cr49_to_p_Mn52_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Cr49_to_p_Mn52_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Cr49_to_p_Mn52_reaclib) = drate_dT;
        }
    }

    {
        // He4_Fe52_to_p_Co55_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Fe52);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe52);
        }
#endif
        rate_He4_Fe52_to_p_Co55_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Fe52_to_p_Co55_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Fe52_to_p_Co55_reaclib) = drate_dT;
        }
    }

    {
        // He4_He4_He4_to_C12_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_He4) + rate_eval.log_screen(k_He4_Be8);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_He4) + rate_eval.dlog_screen_dT(k_He4_Be8);
        }
#endif
        rate_He4_He4_He4_to_C12_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_He4_He4_to_C12_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_He4_He4_to_C12_reaclib) = drate_dT;
        }
    }

    {
        // n_Mn50_to_Mn51_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_n_Mn50_to_Mn51_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_n_Mn50_to_Mn51_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_n_Mn50_to_Mn51_reaclib) = drate_dT;
        }
    }

    {
        // n_Mn51_to_Mn52_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_n_Mn51_to_Mn52_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_n_Mn51_to_Mn52_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_n_Mn51_to_Mn52_reaclib) = drate_dT;
        }
    }

    {
        // p_Mn52_to_Fe53_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Mn52);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn52);
        }
#endif
        rate_p_Mn52_to_Fe53_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Mn52_to_Fe53_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Mn52_to_Fe53_reaclib) = drate_dT;
        }
    }

    {
        // He4_Mn52_to_Co56_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Mn52);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn52);
        }
#endif
        rate_He4_Mn52_to_Co56_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Mn52_to_Co56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Mn52_to_Co56_reaclib) = drate_dT;
        }
    }

    {
        // n_Fe52_to_Fe53_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_n_Fe52_to_Fe53_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_n_Fe52_to_Fe53_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_n_Fe52_to_Fe53_reaclib) = drate_dT;
        }
    }

    {
        // n_Fe53_to_Fe54_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_n_Fe53_to_Fe54_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_n_Fe53_to_Fe54_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_n_Fe53_to_Fe54_reaclib) = drate_dT;
        }
    }

    {
        // He4_Fe53_to_Ni57_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Fe53);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe53);
        }
#endif
        rate_He4_Fe53_to_Ni57_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Fe53_to_Ni57_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Fe53_to_Ni57_reaclib) = drate_dT;
        }
    }

    {
        // n_Fe54_to_Fe55_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_n_Fe54_to_Fe55_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_n_Fe54_to_Fe55_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_n_Fe54_to_Fe55_reaclib) = drate_dT;
        }
    }

    {
        // p_Fe54_to_Co55_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Fe54);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe54);
        }
#endif
        rate_p_Fe54_to_Co55_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Fe54_to_Co55_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Fe54_to_Co55_reaclib) = drate_dT;
        }
    }

    {
        // He4_Fe54_to_Ni58_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Fe54);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe54);
        }
#endif
        rate_He4_Fe54_to_Ni58_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Fe54_to_Ni58_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Fe54_to_Ni58_reaclib) = drate_dT;
        }
    }

    {
        // n_Fe55_to_Fe56_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_n_Fe55_to_Fe56_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_n_Fe55_to_Fe56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_n_Fe55_to_Fe56_reaclib) = drate_dT;
        }
    }

    {
        // p_Fe55_to_Co56_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Fe55);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe55);
        }
#endif
        rate_p_Fe55_to_Co56_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Fe55_to_Co56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Fe55_to_Co56_reaclib) = drate_dT;
        }
    }

    {
        // p_Fe56_to_Co57_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Fe56);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe56);
        }
#endif
        rate_p_Fe56_to_Co57_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Fe56_to_Co57_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Fe56_to_Co57_reaclib) = drate_dT;
        }
    }

    {
        // n_Co55_to_Co56_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_n_Co55_to_Co56_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_n_Co55_to_Co56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_n_Co55_to_Co56_reaclib) = drate_dT;
        }
    }

    {
        // He4_Co55_to_Cu59_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Co55);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co55);
        }
#endif
        rate_He4_Co55_to_Cu59_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Co55_to_Cu59_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Co55_to_Cu59_reaclib) = drate_dT;
        }
    }

    {
        // n_Co56_to_Co57_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_n_Co56_to_Co57_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_n_Co56_to_Co57_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_n_Co56_to_Co57_reaclib) = drate_dT;
        }
    }

    {
        // p_Co56_to_Ni57_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Co56);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co56);
        }
#endif
        rate_p_Co56_to_Ni57_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Co56_to_Ni57_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Co56_to_Ni57_reaclib) = drate_dT;
        }
    }

    {
        // p_Co57_to_Ni58_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Co57);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co57);
        }
#endif
        rate_p_Co57_to_Ni58_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Co57_to_Ni58_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Co57_to_Ni58_reaclib) = drate_dT;
        }
    }

    {
        // n_Ni56_to_Ni57_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_n_Ni56_to_Ni57_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_n_Ni56_to_Ni57_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_n_Ni56_to_Ni57_reaclib) = drate_dT;
        }
    }

    {
        // He4_Ni56_to_Zn60_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Ni56);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni56);
        }
#endif
        rate_He4_Ni56_to_Zn60_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Ni56_to_Zn60_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Ni56_to_Zn60_reaclib) = drate_dT;
        }
    }

    {
        // n_Ni57_to_Ni58_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_n_Ni57_to_Ni58_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_n_Ni57_to_Ni58_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_n_Ni57_to_Ni58_reaclib) = drate_dT;
        }
    }

    {
        // p_Ni58_to_Cu59_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Ni58);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni58);
        }
#endif
        rate_p_Ni58_to_Cu59_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Ni58_to_Cu59_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Ni58_to_Cu59_reaclib) = drate_dT;
        }
    }

    {
        // p_Cu59_to_Zn60_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Cu59);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cu59);
        }
#endif
        rate_p_Cu59_to_Zn60_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Cu59_to_Zn60_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Cu59_to_Zn60_reaclib) = drate_dT;
        }
    }

    {
        // He4_Mn50_to_p_Fe53_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Mn50);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn50);
        }
#endif
        rate_He4_Mn50_to_p_Fe53_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Mn50_to_p_Fe53_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Mn50_to_p_Fe53_reaclib) = drate_dT;
        }
    }

    {
        // He4_Mn51_to_p_Fe54_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Mn51);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn51);
        }
#endif
        rate_He4_Mn51_to_p_Fe54_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Mn51_to_p_Fe54_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Mn51_to_p_Fe54_reaclib) = drate_dT;
        }
    }

    {
        // He4_Mn52_to_p_Fe55_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Mn52);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn52);
        }
#endif
        rate_He4_Mn52_to_p_Fe55_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Mn52_to_p_Fe55_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Mn52_to_p_Fe55_reaclib) = drate_dT;
        }
    }

    {
        // n_Fe52_to_p_Mn52_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_n_Fe52_to_p_Mn52_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_n_Fe52_to_p_Mn52_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_n_Fe52_to_p_Mn52_reaclib) = drate_dT;
        }
    }

    {
        // He4_Fe53_to_p_Co56_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Fe53);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe53);
        }
#endif
        rate_He4_Fe53_to_p_Co56_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Fe53_to_p_Co56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Fe53_to_p_Co56_reaclib) = drate_dT;
        }
    }

    {
        // n_Co55_to_p_Fe55_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_n_Co55_to_p_Fe55_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_n_Co55_to_p_Fe55_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_n_Co55_to_p_Fe55_reaclib) = drate_dT;
        }
    }

    {
        // n_Co55_to_He4_Mn52_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_n_Co55_to_He4_Mn52_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_n_Co55_to_He4_Mn52_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_n_Co55_to_He4_Mn52_reaclib) = drate_dT;
        }
    }

    {
        // He4_Co55_to_p_Ni58_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Co55);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co55);
        }
#endif
        rate_He4_Co55_to_p_Ni58_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_He4_Co55_to_p_Ni58_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Co55_to_p_Ni58_reaclib) = drate_dT;
        }
    }

    {
        // n_Co56_to_p_Fe56_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_n_Co56_to_p_Fe56_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_n_Co56_to_p_Fe56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_n_Co56_to_p_Fe56_reaclib) = drate_dT;
        }
    }

    {
        // p_Co57_to_He4_Fe54_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Co57);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co57);
        }
#endif
        rate_p_Co57_to_He4_Fe54_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Co57_to_He4_Fe54_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Co57_to_He4_Fe54_reaclib) = drate_dT;
        }
    }

    {
        // n_Ni56_to_p_Co56_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_n_Ni56_to_p_Co56_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_n_Ni56_to_p_Co56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_n_Ni56_to_p_Co56_reaclib) = drate_dT;
        }
    }

    {
        // n_Ni56_to_He4_Fe53_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_n_Ni56_to_He4_Fe53_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_n_Ni56_to_He4_Fe53_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_n_Ni56_to_He4_Fe53_reaclib) = drate_dT;
        }
    }

    {
        // n_Ni57_to_p_Co57_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_n_Ni57_to_p_Co57_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_n_Ni57_to_p_Co57_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_n_Ni57_to_p_Co57_reaclib) = drate_dT;
        }
    }

    {
        // n_Ni57_to_He4_Fe54_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_n_Ni57_to_He4_Fe54_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_n_Ni57_to_He4_Fe54_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_n_Ni57_to_He4_Fe54_reaclib) = drate_dT;
        }
    }

    {
        // n_Ni58_to_He4_Fe55_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_n_Ni58_to_He4_Fe55_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_n_Ni58_to_He4_Fe55_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_n_Ni58_to_He4_Fe55_reaclib) = drate_dT;
        }
    }

    {
        // n_Cu59_to_He4_Co56_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_n_Cu59_to_He4_Co56_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_n_Cu59_to_He4_Co56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_n_Cu59_to_He4_Co56_reaclib) = drate_dT;
        }
    }

    {
        // p_Cu59_to_He4_Ni56_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Cu59);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cu59);
        }
#endif
        rate_p_Cu59_to_He4_Ni56_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_p_Cu59_to_He4_Ni56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Cu59_to_He4_Ni56_reaclib) = drate_dT;
        }
    }

    {
        // n_Zn60_to_He4_Ni57_reaclib

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_n_Zn60_to_He4_Ni57_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
        rate_eval.screened_rates(k_n_Zn60_to_He4_Ni57_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_n_Zn60_to_He4_Ni57_reaclib) = drate_dT;
        }
    }


}

#endif
