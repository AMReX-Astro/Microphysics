#ifndef DERIVED_RATES_H
#define DERIVED_RATES_H

#include <tfactors.H>
#include <actual_network.H>
#include <partition_functions.H>
#include <temperature_table_rates.H>
#include <interp_tools.H>

using namespace Rates;

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_N13_to_p_C12_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // N13 --> p + C12

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real C12_log_pf, dC12_log_pf_dT9;

    // setting C12 log(partition function) to 0.0 by default, independent of T
    C12_log_pf = 0.0_rt;
    dC12_log_pf_dT9 = 0.0_rt;

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real N13_log_pf, dN13_log_pf_dT9;

    // setting N13 log(partition function) to 0.0 by default, independent of T
    N13_log_pf = 0.0_rt;
    dN13_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = p_log_pf + C12_log_pf - N13_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dC12_log_pf_dT9 - dN13_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ls09n
    ln_set_rate =  40.059128683693466 + -22.553277271250543 * tfactors.T9i + -13.692 * tfactors.T913i + -0.230881 * tfactors.T913
                         + 4.44362 * tfactors.T9 + -3.15898 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  22.553277271250543 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -13.692 * tfactors.T943i + (1.0/3.0) * -0.230881 * tfactors.T923i
                                  + 4.44362 + (5.0/3.0) * -3.15898 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from ls09r
    ln_set_rate =  40.45372868369347 + -26.331767271250545 * tfactors.T9i + -5.10735 * tfactors.T913i + -2.24111 * tfactors.T913
                         + 0.148883 * tfactors.T9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  26.331767271250545 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -5.10735 * tfactors.T943i + (1.0/3.0) * -2.24111 * tfactors.T923i
                                  + 0.148883;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_N14_to_p_C13_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // N14 --> p + C13

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real C13_log_pf, dC13_log_pf_dT9;

    // setting C13 log(partition function) to 0.0 by default, independent of T
    C13_log_pf = 0.0_rt;
    dC13_log_pf_dT9 = 0.0_rt;

    amrex::Real N14_log_pf, dN14_log_pf_dT9;

    // setting N14 log(partition function) to 0.0 by default, independent of T
    N14_log_pf = 0.0_rt;
    dN14_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = p_log_pf + C13_log_pf - N14_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dC13_log_pf_dT9 - dN14_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from nacrn
    ln_set_rate =  41.72421240281699 + -87.62065170634826 * tfactors.T9i + -13.72 * tfactors.T913i + -0.450018 * tfactors.T913
                         + 3.70823 * tfactors.T9 + -1.70545 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  87.62065170634826 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -13.72 * tfactors.T943i + (1.0/3.0) * -0.450018 * tfactors.T923i
                                  + 3.70823 + (5.0/3.0) * -1.70545 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from nacrr
    ln_set_rate =  37.17241240281699 + -93.40212170634825 * tfactors.T9i + -0.196703 * tfactors.T913
                         + 0.142126 * tfactors.T9 + -0.0238912 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  93.40212170634825 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.196703 * tfactors.T923i
                                  + 0.142126 + (5.0/3.0) * -0.0238912 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from nacrr
    ln_set_rate =  38.39121240281698 + -101.17495170634825 * tfactors.T9i;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  101.17495170634825 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_N15_to_p_C14_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // N15 --> p + C14

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real N15_log_pf, dN15_log_pf_dT9;

    // setting N15 log(partition function) to 0.0 by default, independent of T
    N15_log_pf = 0.0_rt;
    dN15_log_pf_dT9 = 0.0_rt;

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real C14_log_pf, dC14_log_pf_dT9;

    // setting C14 log(partition function) to 0.0 by default, independent of T
    C14_log_pf = 0.0_rt;
    dC14_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = p_log_pf + C14_log_pf - N15_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dC14_log_pf_dT9 - dN15_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from il10r
    ln_set_rate =  40.031567482849354 + -119.97566957362164 * tfactors.T9i + -10.658 * tfactors.T913i + 1.73644 * tfactors.T913
                         + -0.350498 * tfactors.T9 + 0.0279902 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  119.97566957362164 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -10.658 * tfactors.T943i + (1.0/3.0) * 1.73644 * tfactors.T923i
                                  + -0.350498 + (5.0/3.0) * 0.0279902 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10n
    ln_set_rate =  43.04816748284936 + -118.45225957362165 * tfactors.T9i + -13.9619 * tfactors.T913i + -4.34315 * tfactors.T913
                         + 6.64922 * tfactors.T9 + -3.22592 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  118.45225957362165 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -13.9619 * tfactors.T943i + (1.0/3.0) * -4.34315 * tfactors.T923i
                                  + 6.64922 + (5.0/3.0) * -3.22592 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_O16_to_p_N15_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // O16 --> p + N15

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real N15_log_pf, dN15_log_pf_dT9;

    // setting N15 log(partition function) to 0.0 by default, independent of T
    N15_log_pf = 0.0_rt;
    dN15_log_pf_dT9 = 0.0_rt;

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real O16_log_pf, dO16_log_pf_dT9;

    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors.T9, pf_cache, O16_log_pf, dO16_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + N15_log_pf - O16_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dN15_log_pf_dT9 - dO16_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from li10r
    ln_set_rate =  38.86679552635226 + -150.96226378057287 * tfactors.T9i
                         + 0.0459037 * tfactors.T9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  150.96226378057287 * tfactors.T9i * tfactors.T9i
                                  + 0.0459037;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from li10r
    ln_set_rate =  30.912955526352267 + -143.65591378057286 * tfactors.T9i;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  143.65591378057286 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from li10n
    ln_set_rate =  44.33999552635227 + -140.73276378057287 * tfactors.T9i + -15.24 * tfactors.T913i + 0.334926 * tfactors.T913
                         + 4.59088 * tfactors.T9 + -4.78468 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  140.73276378057287 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -15.24 * tfactors.T943i + (1.0/3.0) * 0.334926 * tfactors.T923i
                                  + 4.59088 + (5.0/3.0) * -4.78468 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_O16_to_He4_C12_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // O16 --> He4 + C12

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real C12_log_pf, dC12_log_pf_dT9;

    // setting C12 log(partition function) to 0.0 by default, independent of T
    C12_log_pf = 0.0_rt;
    dC12_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real O16_log_pf, dO16_log_pf_dT9;

    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors.T9, pf_cache, O16_log_pf, dO16_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + C12_log_pf - O16_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dC12_log_pf_dT9 - dO16_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from nac2
    ln_set_rate =  279.29694929711803 + -84.95157686791683 * tfactors.T9i + 103.411 * tfactors.T913i + -420.567 * tfactors.T913
                         + 64.0874 * tfactors.T9 + -12.4624 * tfactors.T953 + 138.803 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  84.95157686791683 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 103.411 * tfactors.T943i + (1.0/3.0) * -420.567 * tfactors.T923i
                                  + 64.0874 + (5.0/3.0) * -12.4624 * tfactors.T923 + 138.803 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from nac2
    ln_set_rate =  94.31554929711803 + -84.50314686791683 * tfactors.T9i + 58.9128 * tfactors.T913i + -148.273 * tfactors.T913
                         + 9.08324 * tfactors.T9 + -0.541041 * tfactors.T953 + 71.8554 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  84.50314686791683 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 58.9128 * tfactors.T943i + (1.0/3.0) * -148.273 * tfactors.T923i
                                  + 9.08324 + (5.0/3.0) * -0.541041 * tfactors.T923 + 71.8554 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_O18_to_He4_C14_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // O18 --> He4 + C14

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real O18_log_pf, dO18_log_pf_dT9;

    // interpolating O18 partition function
    get_partition_function_cached(O18, tfactors.T9, pf_cache, O18_log_pf, dO18_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real C14_log_pf, dC14_log_pf_dT9;

    // setting C14 log(partition function) to 0.0 by default, independent of T
    C14_log_pf = 0.0_rt;
    dC14_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = He4_log_pf + C14_log_pf - O18_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dC14_log_pf_dT9 - dO18_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from il10r
    ln_set_rate =  36.54830162906396 + -82.66688913950057 * tfactors.T9i + -3.83188 * tfactors.T913
                         + 1.64358 * tfactors.T9 + -0.177785 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  82.66688913950057 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -3.83188 * tfactors.T923i
                                  + 1.64358 + (5.0/3.0) * -0.177785 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  0.9124016290639609 + -74.33734913950056 * tfactors.T9i;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  74.33734913950056 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10n
    ln_set_rate =  43.20510162906396 + -72.26858913950056 * tfactors.T9i + -31.7222 * tfactors.T913i + 11.3923 * tfactors.T913
                         + -9.92249 * tfactors.T9 + -2.0 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  72.26858913950056 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -31.7222 * tfactors.T943i + (1.0/3.0) * 11.3923 * tfactors.T923i
                                  + -9.92249 + (5.0/3.0) * -2.0 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_F18_to_p_O17_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // F18 --> p + O17

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real F18_log_pf, dF18_log_pf_dT9;

    // interpolating F18 partition function
    get_partition_function_cached(F18, tfactors.T9, pf_cache, F18_log_pf, dF18_log_pf_dT9);

    amrex::Real O17_log_pf, dO17_log_pf_dT9;

    // interpolating O17 partition function
    get_partition_function_cached(O17, tfactors.T9, pf_cache, O17_log_pf, dO17_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + O17_log_pf - F18_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dO17_log_pf_dT9 - dF18_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from il10r
    ln_set_rate =  33.72287495567065 + -71.29605321275191 * tfactors.T9i + 2.31435 * tfactors.T913
                         + -0.302835 * tfactors.T9 + 0.020133 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  71.29605321275191 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 2.31435 * tfactors.T923i
                                  + -0.302835 + (5.0/3.0) * 0.020133 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  11.255394955670651 + -65.81406921275192 * tfactors.T9i;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  65.81406921275192 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10n
    ln_set_rate =  40.22529495567065 + -65.06777321275192 * tfactors.T9i + -16.4035 * tfactors.T913i + 4.31885 * tfactors.T913
                         + -0.709921 * tfactors.T9 + -2.0 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  65.06777321275192 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -16.4035 * tfactors.T943i + (1.0/3.0) * 4.31885 * tfactors.T923i
                                  + -0.709921 + (5.0/3.0) * -2.0 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_F18_to_He4_N14_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // F18 --> He4 + N14

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real F18_log_pf, dF18_log_pf_dT9;

    // interpolating F18 partition function
    get_partition_function_cached(F18, tfactors.T9, pf_cache, F18_log_pf, dF18_log_pf_dT9);

    amrex::Real N14_log_pf, dN14_log_pf_dT9;

    // setting N14 log(partition function) to 0.0 by default, independent of T
    N14_log_pf = 0.0_rt;
    dN14_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = He4_log_pf + N14_log_pf - F18_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dN14_log_pf_dT9 - dF18_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from il10r
    ln_set_rate =  38.61662473666887 + -62.20224752987261 * tfactors.T9i + -5.6227 * tfactors.T913i;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  62.20224752987261 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -5.6227 * tfactors.T943i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  24.91396273666887 + -56.396987529872604 * tfactors.T9i;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  56.396987529872604 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10n
    ln_set_rate =  46.25102473666887 + -51.236647529872606 * tfactors.T9i + -36.2504 * tfactors.T913i
                         + -5.0 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  51.236647529872606 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -36.2504 * tfactors.T943i
                                  + (5.0/3.0) * -5.0 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_F19_to_p_O18_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // F19 --> p + O18

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real O18_log_pf, dO18_log_pf_dT9;

    // interpolating O18 partition function
    get_partition_function_cached(O18, tfactors.T9, pf_cache, O18_log_pf, dO18_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real F19_log_pf, dF19_log_pf_dT9;

    // interpolating F19 partition function
    get_partition_function_cached(F19, tfactors.T9, pf_cache, F19_log_pf, dF19_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + O18_log_pf - F19_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dO18_log_pf_dT9 - dF19_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from il10n
    ln_set_rate =  42.868088644181995 + -92.7618744782197 * tfactors.T9i + -16.7246 * tfactors.T913i
                         + -3.0 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  92.7618744782197 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -16.7246 * tfactors.T943i
                                  + (5.0/3.0) * -3.0 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  30.219848644181997 + -99.4871744782197 * tfactors.T9i + 3.99059 * tfactors.T913
                         + -0.593127 * tfactors.T9 + 0.0877534 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  99.4871744782197 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 3.99059 * tfactors.T923i
                                  + -0.593127 + (5.0/3.0) * 0.0877534 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  28.027568644181997 + -94.41868447821969 * tfactors.T9i;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  94.41868447821969 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  -12.056811355818002 + -93.0066174782197 * tfactors.T9i;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  93.0066174782197 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_F19_to_He4_N15_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // F19 --> He4 + N15

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real N15_log_pf, dN15_log_pf_dT9;

    // setting N15 log(partition function) to 0.0 by default, independent of T
    N15_log_pf = 0.0_rt;
    dN15_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real F19_log_pf, dF19_log_pf_dT9;

    // interpolating F19 partition function
    get_partition_function_cached(F19, tfactors.T9, pf_cache, F19_log_pf, dF19_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + N15_log_pf - F19_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dN15_log_pf_dT9 - dF19_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from il10r
    ln_set_rate =  -4.0595772096034 + -50.77806404408662 * tfactors.T9i + 35.4292 * tfactors.T913
                         + -5.5767 * tfactors.T9 + 0.441293 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  50.77806404408662 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 35.4292 * tfactors.T923i
                                  + -5.5767 + (5.0/3.0) * 0.441293 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  28.273522790396598 + -53.56282404408662 * tfactors.T9i;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  53.56282404408662 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  15.3204027903966 + -50.75615404408662 * tfactors.T9i;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  50.75615404408662 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10n
    ln_set_rate =  50.130922790396596 + -46.57820404408662 * tfactors.T9i + -36.2324 * tfactors.T913i
                         + -2.0 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  46.57820404408662 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -36.2324 * tfactors.T943i
                                  + (5.0/3.0) * -2.0 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ne20_to_p_F19_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ne20 --> p + F19

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real F19_log_pf, dF19_log_pf_dT9;

    // interpolating F19 partition function
    get_partition_function_cached(F19, tfactors.T9, pf_cache, F19_log_pf, dF19_log_pf_dT9);

    amrex::Real Ne20_log_pf, dNe20_log_pf_dT9;

    // interpolating Ne20 partition function
    get_partition_function_cached(Ne20, tfactors.T9, pf_cache, Ne20_log_pf, dNe20_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + F19_log_pf - Ne20_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dF19_log_pf_dT9 - dNe20_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from nacrr
    ln_set_rate =  18.711255537457795 + -156.78628300750532 * tfactors.T9i + 31.6442 * tfactors.T913i + -58.6563 * tfactors.T913
                         + 67.7365 * tfactors.T9 + -22.9721 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  156.78628300750532 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 31.6442 * tfactors.T943i + (1.0/3.0) * -58.6563 * tfactors.T923i
                                  + 67.7365 + (5.0/3.0) * -22.9721 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from nacrr
    ln_set_rate =  36.72378553745779 + -150.75597300750533 * tfactors.T9i + -11.3832 * tfactors.T913i + 5.47872 * tfactors.T913
                         + -1.07203 * tfactors.T9 + 0.11196 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  150.75597300750533 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -11.3832 * tfactors.T943i + (1.0/3.0) * 5.47872 * tfactors.T923i
                                  + -1.07203 + (5.0/3.0) * 0.11196 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from nacrn
    ln_set_rate =  42.62288553745779 + -149.04214300750533 * tfactors.T9i + -18.116 * tfactors.T913i + -1.4622 * tfactors.T913
                         + 6.95113 * tfactors.T9 + -2.90366 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  149.04214300750533 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -18.116 * tfactors.T943i + (1.0/3.0) * -1.4622 * tfactors.T923i
                                  + 6.95113 + (5.0/3.0) * -2.90366 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ne20_to_He4_O16_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ne20 --> He4 + O16

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real O16_log_pf, dO16_log_pf_dT9;

    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors.T9, pf_cache, O16_log_pf, dO16_log_pf_dT9);

    amrex::Real Ne20_log_pf, dNe20_log_pf_dT9;

    // interpolating Ne20 partition function
    get_partition_function_cached(Ne20, tfactors.T9, pf_cache, Ne20_log_pf, dNe20_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + O16_log_pf - Ne20_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dO16_log_pf_dT9 - dNe20_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from co10r
    ln_set_rate =  28.644822801502126 + -65.24608327101186 * tfactors.T9i;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  65.24608327101186 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from co10n
    ln_set_rate =  48.662112801502126 + -54.887583271011856 * tfactors.T9i + -39.7262 * tfactors.T913i + -0.210799 * tfactors.T913
                         + 0.442879 * tfactors.T9 + -0.0797753 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  54.887583271011856 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -39.7262 * tfactors.T943i + (1.0/3.0) * -0.210799 * tfactors.T923i
                                  + 0.442879 + (5.0/3.0) * -0.0797753 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from co10r
    ln_set_rate =  34.267592801502126 + -67.65188327101185 * tfactors.T9i + -3.65925 * tfactors.T913
                         + 0.714224 * tfactors.T9 + -0.00107508 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  67.65188327101185 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -3.65925 * tfactors.T923i
                                  + 0.714224 + (5.0/3.0) * -0.00107508 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ne21_to_He4_O17_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ne21 --> He4 + O17

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Ne21_log_pf, dNe21_log_pf_dT9;

    // interpolating Ne21 partition function
    get_partition_function_cached(Ne21, tfactors.T9, pf_cache, Ne21_log_pf, dNe21_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real O17_log_pf, dO17_log_pf_dT9;

    // interpolating O17 partition function
    get_partition_function_cached(O17, tfactors.T9, pf_cache, O17_log_pf, dO17_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + O17_log_pf - Ne21_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dO17_log_pf_dT9 - dNe21_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from be13r
    ln_set_rate =  0.09298411543850094 + -90.7784662204353 * tfactors.T9i + 123.363 * tfactors.T913i + -87.4351 * tfactors.T913
                         + -3.40974e-06 * tfactors.T9 + -57.0469 * tfactors.T953 + 83.7218 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  90.7784662204353 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 123.363 * tfactors.T943i + (1.0/3.0) * -87.4351 * tfactors.T923i
                                  + -3.40974e-06 + (5.0/3.0) * -57.0469 * tfactors.T923 + 83.7218 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from be13r
    ln_set_rate =  -91.95121588456149 + -98.9451062204353 * tfactors.T9i + 3.31162e-08 * tfactors.T913i + 130.258 * tfactors.T913
                         + -7.92551e-05 * tfactors.T9 + -4.13772 * tfactors.T953 + -41.2753 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  98.9451062204353 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 3.31162e-08 * tfactors.T943i + (1.0/3.0) * 130.258 * tfactors.T923i
                                  + -7.92551e-05 + (5.0/3.0) * -4.13772 * tfactors.T923 + -41.2753 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from be13r
    ln_set_rate =  27.322784115438502 + -91.26872622043531 * tfactors.T9i + 2.87641 * tfactors.T913i + -3.54489 * tfactors.T913
                         + -2.11222e-08 * tfactors.T9 + -3.90649e-09 * tfactors.T953 + 6.25778 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  91.26872622043531 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 2.87641 * tfactors.T943i + (1.0/3.0) * -3.54489 * tfactors.T923i
                                  + -2.11222e-08 + (5.0/3.0) * -3.90649e-09 * tfactors.T923 + 6.25778 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ne22_to_He4_O18_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ne22 --> He4 + O18

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real O18_log_pf, dO18_log_pf_dT9;

    // interpolating O18 partition function
    get_partition_function_cached(O18, tfactors.T9, pf_cache, O18_log_pf, dO18_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ne22_log_pf, dNe22_log_pf_dT9;

    // interpolating Ne22 partition function
    get_partition_function_cached(Ne22, tfactors.T9, pf_cache, Ne22_log_pf, dNe22_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + O18_log_pf - Ne22_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dO18_log_pf_dT9 - dNe22_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from il10r
    ln_set_rate =  -7.119065921580848 + -114.18179645973349 * tfactors.T9i;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  114.18179645973349 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  -56.510065921580846 + -112.85484845973349 * tfactors.T9i;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  112.85484845973349 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  39.76833407841915 + -143.2255364597335 * tfactors.T9i;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  143.2255364597335 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  106.99883407841915 + -113.76407645973349 * tfactors.T9i + -44.3823 * tfactors.T913i + -46.6617 * tfactors.T913
                         + 7.88059 * tfactors.T9 + -0.590829 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  113.76407645973349 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -44.3823 * tfactors.T943i + (1.0/3.0) * -46.6617 * tfactors.T923i
                                  + 7.88059 + (5.0/3.0) * -0.590829 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Na22_to_p_Ne21_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Na22 --> p + Ne21

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Ne21_log_pf, dNe21_log_pf_dT9;

    // interpolating Ne21 partition function
    get_partition_function_cached(Ne21, tfactors.T9, pf_cache, Ne21_log_pf, dNe21_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Na22_log_pf, dNa22_log_pf_dT9;

    // interpolating Na22 partition function
    get_partition_function_cached(Na22, tfactors.T9, pf_cache, Na22_log_pf, dNa22_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Ne21_log_pf - Na22_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dNe21_log_pf_dT9 - dNa22_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from il10n
    ln_set_rate =  42.165346980750385 + -78.19798607071328 * tfactors.T9i + -19.2096 * tfactors.T913i
                         + -1.0 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  78.19798607071328 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -19.2096 * tfactors.T943i
                                  + (5.0/3.0) * -1.0 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  -16.390453019249613 + -82.41183607071328 * tfactors.T9i + 21.1176 * tfactors.T913i + 34.0411 * tfactors.T913
                         + -4.45593 * tfactors.T9 + 0.328613 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  82.41183607071328 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 21.1176 * tfactors.T943i + (1.0/3.0) * 34.0411 * tfactors.T923i
                                  + -4.45593 + (5.0/3.0) * 0.328613 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  24.852786980750384 + -79.59755607071328 * tfactors.T9i;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  79.59755607071328 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  -24.559653019249616 + -78.39416607071328 * tfactors.T9i;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  78.39416607071328 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Na22_to_He4_F18_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Na22 --> He4 + F18

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real F18_log_pf, dF18_log_pf_dT9;

    // interpolating F18 partition function
    get_partition_function_cached(F18, tfactors.T9, pf_cache, F18_log_pf, dF18_log_pf_dT9);

    amrex::Real Na22_log_pf, dNa22_log_pf_dT9;

    // interpolating Na22 partition function
    get_partition_function_cached(Na22, tfactors.T9, pf_cache, Na22_log_pf, dNa22_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + F18_log_pf - Na22_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dF18_log_pf_dT9 - dNa22_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from rpsmr
    ln_set_rate =  59.32473614051823 + -100.22898907841348 * tfactors.T9i + 18.8956 * tfactors.T913i + -65.6134 * tfactors.T913
                         + 1.71114 * tfactors.T9 + -0.0260999 * tfactors.T953 + 39.3396 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  100.22898907841348 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 18.8956 * tfactors.T943i + (1.0/3.0) * -65.6134 * tfactors.T923i
                                  + 1.71114 + (5.0/3.0) * -0.0260999 * tfactors.T923 + 39.3396 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Na23_to_p_Ne22_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Na23 --> p + Ne22

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Na23_log_pf, dNa23_log_pf_dT9;

    // interpolating Na23 partition function
    get_partition_function_cached(Na23, tfactors.T9, pf_cache, Na23_log_pf, dNa23_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Ne22_log_pf, dNe22_log_pf_dT9;

    // interpolating Ne22 partition function
    get_partition_function_cached(Ne22, tfactors.T9, pf_cache, Ne22_log_pf, dNe22_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Ne22_log_pf - Na23_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dNe22_log_pf_dT9 - dNa23_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ke17r
    ln_set_rate =  21.673073381706992 + -103.77146219384292 * tfactors.T9i
                         + 1.189235 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  103.77146219384292 * tfactors.T9i * tfactors.T9i
                                  + 1.189235 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from ke17r
    ln_set_rate =  0.8379043817069913 + -102.46235419384291 * tfactors.T9i
                         + 0.009810000000000096 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  102.46235419384291 * tfactors.T9i * tfactors.T9i
                                  + 0.009810000000000096 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from ke17r
    ln_set_rate =  18.18205438170699 + -102.85114819384292 * tfactors.T9i
                         + 4.735580000000001 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  102.85114819384292 * tfactors.T9i * tfactors.T9i
                                  + 4.735580000000001 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from ke17r
    ln_set_rate =  36.30970438170699 + -110.77516219384292 * tfactors.T9i
                         + 0.732533 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  110.77516219384292 * tfactors.T9i * tfactors.T9i
                                  + 0.732533 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from ke17r
    ln_set_rate =  33.91320438170699 + -106.65075219384292 * tfactors.T9i
                         + 1.656226 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  106.65075219384292 * tfactors.T9i * tfactors.T9i
                                  + 1.656226 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from ke17r
    ln_set_rate =  18.266434381706993 + -104.66929219384292 * tfactors.T9i
                         + -2.79964 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  104.66929219384292 * tfactors.T9i * tfactors.T9i
                                  + -2.79964 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Na23_to_He4_F19_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Na23 --> He4 + F19

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Na23_log_pf, dNa23_log_pf_dT9;

    // interpolating Na23 partition function
    get_partition_function_cached(Na23, tfactors.T9, pf_cache, Na23_log_pf, dNa23_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real F19_log_pf, dF19_log_pf_dT9;

    // interpolating F19 partition function
    get_partition_function_cached(F19, tfactors.T9, pf_cache, F19_log_pf, dF19_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + F19_log_pf - Na23_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dF19_log_pf_dT9 - dNa23_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from rpsmr
    ln_set_rate =  76.90044981594414 + -123.58233417534949 * tfactors.T9i + 39.7219 * tfactors.T913i + -100.401 * tfactors.T913
                         + 3.15808 * tfactors.T9 + -0.0629822 * tfactors.T953 + 55.9823 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  123.58233417534949 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 39.7219 * tfactors.T943i + (1.0/3.0) * -100.401 * tfactors.T923i
                                  + 3.15808 + (5.0/3.0) * -0.0629822 * tfactors.T923 + 55.9823 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Mg24_to_p_Na23_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mg24 --> p + Na23

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Na23_log_pf, dNa23_log_pf_dT9;

    // interpolating Na23 partition function
    get_partition_function_cached(Na23, tfactors.T9, pf_cache, Na23_log_pf, dNa23_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Mg24_log_pf, dMg24_log_pf_dT9;

    // interpolating Mg24 partition function
    get_partition_function_cached(Mg24, tfactors.T9, pf_cache, Mg24_log_pf, dMg24_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Na23_log_pf - Mg24_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dNa23_log_pf_dT9 - dMg24_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from il10r
    ln_set_rate =  34.107544561463015 + -138.96838756220697 * tfactors.T9i + -0.360588 * tfactors.T913
                         + 1.4187 * tfactors.T9 + -0.184061 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  138.96838756220697 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.360588 * tfactors.T923i
                                  + 1.4187 + (5.0/3.0) * -0.184061 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  20.02229456146302 + -137.30028756220696 * tfactors.T9i;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  137.30028756220696 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10n
    ln_set_rate =  43.95564456146302 + -135.68809756220696 * tfactors.T9i + -20.6428 * tfactors.T913i + 1.52954 * tfactors.T913
                         + 2.7487 * tfactors.T9 + -1.0 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  135.68809756220696 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -20.6428 * tfactors.T943i + (1.0/3.0) * 1.52954 * tfactors.T923i
                                  + 2.7487 + (5.0/3.0) * -1.0 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Mg24_to_He4_Ne20_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mg24 --> He4 + Ne20

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Mg24_log_pf, dMg24_log_pf_dT9;

    // interpolating Mg24 partition function
    get_partition_function_cached(Mg24, tfactors.T9, pf_cache, Mg24_log_pf, dMg24_log_pf_dT9);

    amrex::Real Ne20_log_pf, dNe20_log_pf_dT9;

    // interpolating Ne20 partition function
    get_partition_function_cached(Ne20, tfactors.T9, pf_cache, Ne20_log_pf, dNe20_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Ne20_log_pf - Mg24_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dNe20_log_pf_dT9 - dMg24_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from il10r
    ln_set_rate =  16.022538839949377 + -120.89510873006793 * tfactors.T9i + 16.9229 * tfactors.T913
                         + -2.57325 * tfactors.T9 + 0.208997 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  120.89510873006793 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 16.9229 * tfactors.T923i
                                  + -2.57325 + (5.0/3.0) * 0.208997 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  26.803878839949377 + -117.33446873006793 * tfactors.T9i;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  117.33446873006793 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  -13.884691160050625 + -110.62025873006793 * tfactors.T9i;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  110.62025873006793 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10n
    ln_set_rate =  49.32660883994937 + -108.11420873006793 * tfactors.T9i + -46.2525 * tfactors.T913i + 5.58901 * tfactors.T913
                         + 7.61843 * tfactors.T9 + -3.683 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  108.11420873006793 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -46.2525 * tfactors.T943i + (1.0/3.0) * 5.58901 * tfactors.T923i
                                  + 7.61843 + (5.0/3.0) * -3.683 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Mg25_to_He4_Ne21_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mg25 --> He4 + Ne21

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Ne21_log_pf, dNe21_log_pf_dT9;

    // interpolating Ne21 partition function
    get_partition_function_cached(Ne21, tfactors.T9, pf_cache, Ne21_log_pf, dNe21_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Mg25_log_pf, dMg25_log_pf_dT9;

    // interpolating Mg25 partition function
    get_partition_function_cached(Mg25, tfactors.T9, pf_cache, Mg25_log_pf, dMg25_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Ne21_log_pf - Mg25_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dNe21_log_pf_dT9 - dMg25_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from cf88r
    ln_set_rate =  50.67026655911652 + -136.77028986154426 * tfactors.T9i + -29.4583 * tfactors.T913
                         + 14.6328 * tfactors.T9 + -3.47392 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  136.77028986154426 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -29.4583 * tfactors.T923i
                                  + 14.6328 + (5.0/3.0) * -3.47392 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from cf88n
    ln_set_rate =  61.12006655911652 + -114.72128986154426 * tfactors.T9i + -46.89 * tfactors.T913i + -0.72642 * tfactors.T913
                         + -0.76406 * tfactors.T9 + 0.0797483 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  114.72128986154426 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -46.89 * tfactors.T943i + (1.0/3.0) * -0.72642 * tfactors.T923i
                                  + -0.76406 + (5.0/3.0) * 0.0797483 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Al26_to_p_Mg25_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Al26 --> p + Mg25

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Al26_log_pf, dAl26_log_pf_dT9;

    // interpolating Al26 partition function
    get_partition_function_cached(Al26, tfactors.T9, pf_cache, Al26_log_pf, dAl26_log_pf_dT9);

    amrex::Real Mg25_log_pf, dMg25_log_pf_dT9;

    // interpolating Mg25 partition function
    get_partition_function_cached(Mg25, tfactors.T9, pf_cache, Mg25_log_pf, dMg25_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Mg25_log_pf - Al26_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dMg25_log_pf_dT9 - dAl26_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from il10r
    ln_set_rate =  27.27838039741493 + -73.90176311269589 * tfactors.T9i + -88.9297 * tfactors.T913
                         + 302.948 * tfactors.T9 + -346.461 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  73.90176311269589 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -88.9297 * tfactors.T923i
                                  + 302.948 + (5.0/3.0) * -346.461 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  -14.13617960258507 + -73.61129911269589 * tfactors.T9i;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  73.61129911269589 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  25.28790039741493 + -76.40546311269588 * tfactors.T9i + 8.46334 * tfactors.T913
                         + -0.907024 * tfactors.T9 + 0.0642981 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  76.40546311269588 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 8.46334 * tfactors.T923i
                                  + -0.907024 + (5.0/3.0) * 0.0642981 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Al26_to_He4_Na22_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Al26 --> He4 + Na22

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Al26_log_pf, dAl26_log_pf_dT9;

    // interpolating Al26 partition function
    get_partition_function_cached(Al26, tfactors.T9, pf_cache, Al26_log_pf, dAl26_log_pf_dT9);

    amrex::Real Na22_log_pf, dNa22_log_pf_dT9;

    // interpolating Na22 partition function
    get_partition_function_cached(Na22, tfactors.T9, pf_cache, Na22_log_pf, dNa22_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Na22_log_pf - Al26_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dNa22_log_pf_dT9 - dAl26_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  60.771439975781064 + -109.70523690351963 * tfactors.T9i + -50.0924 * tfactors.T913i + -0.390826 * tfactors.T913
                         + -0.99531 * tfactors.T9 + 0.101354 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  109.70523690351963 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -50.0924 * tfactors.T943i + (1.0/3.0) * -0.390826 * tfactors.T923i
                                  + -0.99531 + (5.0/3.0) * 0.101354 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Al27_to_He4_Na23_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Al27 --> He4 + Na23

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Na23_log_pf, dNa23_log_pf_dT9;

    // interpolating Na23 partition function
    get_partition_function_cached(Na23, tfactors.T9, pf_cache, Na23_log_pf, dNa23_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Al27_log_pf, dAl27_log_pf_dT9;

    // interpolating Al27 partition function
    get_partition_function_cached(Al27, tfactors.T9, pf_cache, Al27_log_pf, dAl27_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Na23_log_pf - Al27_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dNa23_log_pf_dT9 - dAl27_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  69.2206963145215 + -117.11189602392653 * tfactors.T9i + -50.2042 * tfactors.T913i + -1.64239 * tfactors.T913
                         + -1.59995 * tfactors.T9 + 0.184933 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  117.11189602392653 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -50.2042 * tfactors.T943i + (1.0/3.0) * -1.64239 * tfactors.T923i
                                  + -1.59995 + (5.0/3.0) * 0.184933 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Si28_to_p_Al27_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Si28 --> p + Al27

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Al27_log_pf, dAl27_log_pf_dT9;

    // interpolating Al27 partition function
    get_partition_function_cached(Al27, tfactors.T9, pf_cache, Al27_log_pf, dAl27_log_pf_dT9);

    amrex::Real Si28_log_pf, dSi28_log_pf_dT9;

    // interpolating Si28 partition function
    get_partition_function_cached(Si28, tfactors.T9, pf_cache, Si28_log_pf, dSi28_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Al27_log_pf - Si28_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dAl27_log_pf_dT9 - dSi28_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from il10r
    ln_set_rate =  11.796402252999236 + -136.3412367256539 * tfactors.T9i + 23.8634 * tfactors.T913
                         + -3.70135 * tfactors.T9 + 0.28964 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  136.3412367256539 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 23.8634 * tfactors.T923i
                                  + -3.70135 + (5.0/3.0) * 0.28964 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  111.48620225299923 + -134.8245897256539 * tfactors.T9i + -26.8327 * tfactors.T913i + -116.137 * tfactors.T913
                         + 0.00950567 * tfactors.T9 + 0.00999755 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  134.8245897256539 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -26.8327 * tfactors.T943i + (1.0/3.0) * -116.137 * tfactors.T923i
                                  + 0.00950567 + (5.0/3.0) * 0.00999755 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10n
    ln_set_rate =  46.569302252999236 + -134.4372767256539 * tfactors.T9i + -23.2205 * tfactors.T913i
                         + -2.0 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  134.4372767256539 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -23.2205 * tfactors.T943i
                                  + (5.0/3.0) * -2.0 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Si28_to_He4_Mg24_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Si28 --> He4 + Mg24

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Mg24_log_pf, dMg24_log_pf_dT9;

    // interpolating Mg24 partition function
    get_partition_function_cached(Mg24, tfactors.T9, pf_cache, Mg24_log_pf, dMg24_log_pf_dT9);

    amrex::Real Si28_log_pf, dSi28_log_pf_dT9;

    // interpolating Si28 partition function
    get_partition_function_cached(Si28, tfactors.T9, pf_cache, Si28_log_pf, dSi28_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Mg24_log_pf - Si28_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dMg24_log_pf_dT9 - dSi28_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from st08r
    ln_set_rate =  32.902724006057724 + -131.4900751873663 * tfactors.T9i;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  131.4900751873663 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from st08r
    ln_set_rate =  -25.68644599394228 + -128.6942751873663 * tfactors.T9i + 21.3721 * tfactors.T913i + 37.7649 * tfactors.T913
                         + -4.10635 * tfactors.T9 + 0.249618 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  128.6942751873663 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 21.3721 * tfactors.T943i + (1.0/3.0) * 37.7649 * tfactors.T923i
                                  + -4.10635 + (5.0/3.0) * 0.249618 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Si29_to_He4_Mg25_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Si29 --> He4 + Mg25

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Si29_log_pf, dSi29_log_pf_dT9;

    // interpolating Si29 partition function
    get_partition_function_cached(Si29, tfactors.T9, pf_cache, Si29_log_pf, dSi29_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Mg25_log_pf, dMg25_log_pf_dT9;

    // interpolating Mg25 partition function
    get_partition_function_cached(Mg25, tfactors.T9, pf_cache, Mg25_log_pf, dMg25_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Mg25_log_pf - Si29_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dMg25_log_pf_dT9 - dSi29_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from cf88n
    ln_set_rate =  66.34177040898071 + -129.1259956125848 * tfactors.T9i + -53.41 * tfactors.T913i + -1.83266 * tfactors.T913
                         + -0.573073 * tfactors.T9 + 0.0462678 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  129.1259956125848 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -53.41 * tfactors.T943i + (1.0/3.0) * -1.83266 * tfactors.T923i
                                  + -0.573073 + (5.0/3.0) * 0.0462678 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_P30_to_p_Si29_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // P30 --> p + Si29

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Si29_log_pf, dSi29_log_pf_dT9;

    // interpolating Si29 partition function
    get_partition_function_cached(Si29, tfactors.T9, pf_cache, Si29_log_pf, dSi29_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real P30_log_pf, dP30_log_pf_dT9;

    // interpolating P30 partition function
    get_partition_function_cached(P30, tfactors.T9, pf_cache, P30_log_pf, dP30_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Si29_log_pf - P30_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dSi29_log_pf_dT9 - dP30_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from il10r
    ln_set_rate =  22.020796544549516 + -68.26366270095757 * tfactors.T9i + 14.0921 * tfactors.T913
                         + -3.92096 * tfactors.T9 + 0.447706 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  68.26366270095757 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 14.0921 * tfactors.T923i
                                  + -3.92096 + (5.0/3.0) * 0.447706 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  9.798606544549516 + -66.17463270095756 * tfactors.T9i;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  66.17463270095756 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10n
    ln_set_rate =  39.78690654454952 + -64.92437270095756 * tfactors.T9i + -23.9101 * tfactors.T913i + 10.7796 * tfactors.T913
                         + -3.04181 * tfactors.T9 + 0.274565 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  64.92437270095756 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -23.9101 * tfactors.T943i + (1.0/3.0) * 10.7796 * tfactors.T923i
                                  + -3.04181 + (5.0/3.0) * 0.274565 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_P30_to_He4_Al26_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // P30 --> He4 + Al26

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Al26_log_pf, dAl26_log_pf_dT9;

    // interpolating Al26 partition function
    get_partition_function_cached(Al26, tfactors.T9, pf_cache, Al26_log_pf, dAl26_log_pf_dT9);

    amrex::Real P30_log_pf, dP30_log_pf_dT9;

    // interpolating P30 partition function
    get_partition_function_cached(P30, tfactors.T9, pf_cache, P30_log_pf, dP30_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Al26_log_pf - P30_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dAl26_log_pf_dT9 - dP30_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  69.1566565561153 + -120.86843520086329 * tfactors.T9i + -56.4422 * tfactors.T913i + -2.44848 * tfactors.T913
                         + -1.17578 * tfactors.T9 + 0.150757 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  120.86843520086329 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -56.4422 * tfactors.T943i + (1.0/3.0) * -2.44848 * tfactors.T923i
                                  + -1.17578 + (5.0/3.0) * 0.150757 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_P31_to_He4_Al27_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // P31 --> He4 + Al27

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Al27_log_pf, dAl27_log_pf_dT9;

    // interpolating Al27 partition function
    get_partition_function_cached(Al27, tfactors.T9, pf_cache, Al27_log_pf, dAl27_log_pf_dT9);

    amrex::Real P31_log_pf, dP31_log_pf_dT9;

    // interpolating P31 partition function
    get_partition_function_cached(P31, tfactors.T9, pf_cache, P31_log_pf, dP31_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Al27_log_pf - P31_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dAl27_log_pf_dT9 - dP31_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  73.21877134266595 + -112.19944704325766 * tfactors.T9i + -56.5351 * tfactors.T913i + -0.896208 * tfactors.T913
                         + -1.72024 * tfactors.T9 + 0.185409 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  112.19944704325766 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -56.5351 * tfactors.T943i + (1.0/3.0) * -0.896208 * tfactors.T923i
                                  + -1.72024 + (5.0/3.0) * 0.185409 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_S32_to_p_P31_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // S32 --> p + P31

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real S32_log_pf, dS32_log_pf_dT9;

    // interpolating S32 partition function
    get_partition_function_cached(S32, tfactors.T9, pf_cache, S32_log_pf, dS32_log_pf_dT9);

    amrex::Real P31_log_pf, dP31_log_pf_dT9;

    // interpolating P31 partition function
    get_partition_function_cached(P31, tfactors.T9, pf_cache, P31_log_pf, dP31_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + P31_log_pf - S32_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dP31_log_pf_dT9 - dS32_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from il10r
    ln_set_rate =  25.192389422303183 + -106.63906812813134 * tfactors.T9i + 8.09341 * tfactors.T913
                         + -0.615971 * tfactors.T9 + 0.031159 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  106.63906812813134 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 8.09341 * tfactors.T923i
                                  + -0.615971 + (5.0/3.0) * 0.031159 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  21.702443422303183 + -105.12160812813134 * tfactors.T9i;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  105.12160812813134 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10n
    ln_set_rate =  43.630433422303184 + -102.86202812813134 * tfactors.T9i + -25.3278 * tfactors.T913i + 6.4931 * tfactors.T913
                         + -9.27513 * tfactors.T9 + -0.610439 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  102.86202812813134 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -25.3278 * tfactors.T943i + (1.0/3.0) * 6.4931 * tfactors.T923i
                                  + -9.27513 + (5.0/3.0) * -0.610439 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_S32_to_He4_Si28_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // S32 --> He4 + Si28

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real S32_log_pf, dS32_log_pf_dT9;

    // interpolating S32 partition function
    get_partition_function_cached(S32, tfactors.T9, pf_cache, S32_log_pf, dS32_log_pf_dT9);

    amrex::Real Si28_log_pf, dSi28_log_pf_dT9;

    // interpolating Si28 partition function
    get_partition_function_cached(Si28, tfactors.T9, pf_cache, Si28_log_pf, dSi28_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Si28_log_pf - S32_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dSi28_log_pf_dT9 - dS32_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  72.8147025119699 + -80.62419844573513 * tfactors.T9i + -59.4896 * tfactors.T913i + 4.47205 * tfactors.T913
                         + -4.78989 * tfactors.T9 + 0.557201 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  80.62419844573513 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -59.4896 * tfactors.T943i + (1.0/3.0) * 4.47205 * tfactors.T923i
                                  + -4.78989 + (5.0/3.0) * 0.557201 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_S33_to_He4_Si29_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // S33 --> He4 + Si29

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Si29_log_pf, dSi29_log_pf_dT9;

    // interpolating Si29 partition function
    get_partition_function_cached(Si29, tfactors.T9, pf_cache, Si29_log_pf, dSi29_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real S33_log_pf, dS33_log_pf_dT9;

    // interpolating S33 partition function
    get_partition_function_cached(S33, tfactors.T9, pf_cache, S33_log_pf, dS33_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Si29_log_pf - S33_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dSi29_log_pf_dT9 - dS33_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  73.77254890676622 + -82.57418453637936 * tfactors.T9i + -59.5755 * tfactors.T913i + 1.06274 * tfactors.T913
                         + -3.07529 * tfactors.T9 + 0.372011 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  82.57418453637936 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -59.5755 * tfactors.T943i + (1.0/3.0) * 1.06274 * tfactors.T923i
                                  + -3.07529 + (5.0/3.0) * 0.372011 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Cl34_to_p_S33_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Cl34 --> p + S33

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Cl34_log_pf, dCl34_log_pf_dT9;

    // interpolating Cl34 partition function
    get_partition_function_cached(Cl34, tfactors.T9, pf_cache, Cl34_log_pf, dCl34_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real S33_log_pf, dS33_log_pf_dT9;

    // interpolating S33 partition function
    get_partition_function_cached(S33, tfactors.T9, pf_cache, S33_log_pf, dS33_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + S33_log_pf - Cl34_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dS33_log_pf_dT9 - dCl34_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  61.557325011533386 + -59.684389677674176 * tfactors.T9i + -26.777 * tfactors.T913i + -5.96882 * tfactors.T913
                         + -1.0706 * tfactors.T9 + 0.19692 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  59.684389677674176 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -26.777 * tfactors.T943i + (1.0/3.0) * -5.96882 * tfactors.T923i
                                  + -1.0706 + (5.0/3.0) * 0.19692 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Cl34_to_He4_P30_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Cl34 --> He4 + P30

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Cl34_log_pf, dCl34_log_pf_dT9;

    // interpolating Cl34 partition function
    get_partition_function_cached(Cl34, tfactors.T9, pf_cache, Cl34_log_pf, dCl34_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real P30_log_pf, dP30_log_pf_dT9;

    // interpolating P30 partition function
    get_partition_function_cached(P30, tfactors.T9, pf_cache, P30_log_pf, dP30_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + P30_log_pf - Cl34_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dP30_log_pf_dT9 - dCl34_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  71.33676737375008 + -77.33420151312718 * tfactors.T9i + -62.4643 * tfactors.T913i + -3.19028 * tfactors.T913
                         + -0.832633 * tfactors.T9 + 0.0987525 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  77.33420151312718 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -62.4643 * tfactors.T943i + (1.0/3.0) * -3.19028 * tfactors.T923i
                                  + -0.832633 + (5.0/3.0) * 0.0987525 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Cl35_to_He4_P31_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Cl35 --> He4 + P31

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Cl35_log_pf, dCl35_log_pf_dT9;

    // interpolating Cl35 partition function
    get_partition_function_cached(Cl35, tfactors.T9, pf_cache, Cl35_log_pf, dCl35_log_pf_dT9);

    amrex::Real P31_log_pf, dP31_log_pf_dT9;

    // interpolating P31 partition function
    get_partition_function_cached(P31, tfactors.T9, pf_cache, P31_log_pf, dP31_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + P31_log_pf - Cl35_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dP31_log_pf_dT9 - dCl35_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  74.6695981567728 + -81.20727442147438 * tfactors.T9i + -62.5433 * tfactors.T913i + -2.95026 * tfactors.T913
                         + -0.89652 * tfactors.T9 + 0.0774126 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  81.20727442147438 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -62.5433 * tfactors.T943i + (1.0/3.0) * -2.95026 * tfactors.T923i
                                  + -0.89652 + (5.0/3.0) * 0.0774126 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ar36_to_p_Cl35_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ar36 --> p + Cl35

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Cl35_log_pf, dCl35_log_pf_dT9;

    // interpolating Cl35 partition function
    get_partition_function_cached(Cl35, tfactors.T9, pf_cache, Cl35_log_pf, dCl35_log_pf_dT9);

    amrex::Real Ar36_log_pf, dAr36_log_pf_dT9;

    // interpolating Ar36 partition function
    get_partition_function_cached(Ar36, tfactors.T9, pf_cache, Ar36_log_pf, dAr36_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Cl35_log_pf - Ar36_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dCl35_log_pf_dT9 - dAr36_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from il10r
    ln_set_rate =  -17.455589662878502 + -99.28409012589869 * tfactors.T9i;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  99.28409012589869 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10n
    ln_set_rate =  60.7561103371215 + -98.71943912589869 * tfactors.T9i + -27.8971 * tfactors.T913i + -16.2304 * tfactors.T913
                         + 35.255 * tfactors.T9 + -25.8411 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  98.71943912589869 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -27.8971 * tfactors.T943i + (1.0/3.0) * -16.2304 * tfactors.T923i
                                  + 35.255 + (5.0/3.0) * -25.8411 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  17.222320337121502 + -102.37035912589869 * tfactors.T9i + 18.0179 * tfactors.T913
                         + -2.86304 * tfactors.T9 + 0.250854 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  102.37035912589869 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 18.0179 * tfactors.T923i
                                  + -2.86304 + (5.0/3.0) * 0.250854 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  16.0363703371215 + -100.7293991258987 * tfactors.T9i;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  100.7293991258987 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ar36_to_He4_S32_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ar36 --> He4 + S32

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real S32_log_pf, dS32_log_pf_dT9;

    // interpolating S32 partition function
    get_partition_function_cached(S32, tfactors.T9, pf_cache, S32_log_pf, dS32_log_pf_dT9);

    amrex::Real Ar36_log_pf, dAr36_log_pf_dT9;

    // interpolating Ar36 partition function
    get_partition_function_cached(Ar36, tfactors.T9, pf_cache, Ar36_log_pf, dAr36_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + S32_log_pf - Ar36_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dS32_log_pf_dT9 - dAr36_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  73.81807507159112 + -77.06468541928012 * tfactors.T9i + -65.3709 * tfactors.T913i + 5.68294 * tfactors.T913
                         + -5.00388 * tfactors.T9 + 0.571407 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  77.06468541928012 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -65.3709 * tfactors.T943i + (1.0/3.0) * 5.68294 * tfactors.T923i
                                  + -5.00388 + (5.0/3.0) * 0.571407 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ar37_to_He4_S33_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ar37 --> He4 + S33

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Ar37_log_pf, dAr37_log_pf_dT9;

    // interpolating Ar37 partition function
    get_partition_function_cached(Ar37, tfactors.T9, pf_cache, Ar37_log_pf, dAr37_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real S33_log_pf, dS33_log_pf_dT9;

    // interpolating S33 partition function
    get_partition_function_cached(S33, tfactors.T9, pf_cache, S33_log_pf, dS33_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + S33_log_pf - Ar37_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dS33_log_pf_dT9 - dAr37_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  74.85364487556038 + -78.75681911725914 * tfactors.T9i + -65.4446 * tfactors.T913i + 3.59607 * tfactors.T913
                         + -3.40501 * tfactors.T9 + 0.363961 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  78.75681911725914 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -65.4446 * tfactors.T943i + (1.0/3.0) * 3.59607 * tfactors.T923i
                                  + -3.40501 + (5.0/3.0) * 0.363961 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_K38_to_p_Ar37_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // K38 --> p + Ar37

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Ar37_log_pf, dAr37_log_pf_dT9;

    // interpolating Ar37 partition function
    get_partition_function_cached(Ar37, tfactors.T9, pf_cache, Ar37_log_pf, dAr37_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real K38_log_pf, dK38_log_pf_dT9;

    // interpolating K38 partition function
    get_partition_function_cached(K38, tfactors.T9, pf_cache, K38_log_pf, dK38_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Ar37_log_pf - K38_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dAr37_log_pf_dT9 - dK38_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  60.85081421447474 + -59.67102475410763 * tfactors.T9i + -28.9954 * tfactors.T913i + -5.57966 * tfactors.T913
                         + -1.27516 * tfactors.T9 + 0.217568 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  59.67102475410763 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -28.9954 * tfactors.T943i + (1.0/3.0) * -5.57966 * tfactors.T923i
                                  + -1.27516 + (5.0/3.0) * 0.217568 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_K38_to_He4_Cl34_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // K38 --> He4 + Cl34

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Cl34_log_pf, dCl34_log_pf_dT9;

    // interpolating Cl34 partition function
    get_partition_function_cached(Cl34, tfactors.T9, pf_cache, Cl34_log_pf, dCl34_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real K38_log_pf, dK38_log_pf_dT9;

    // interpolating K38 partition function
    get_partition_function_cached(K38, tfactors.T9, pf_cache, K38_log_pf, dK38_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Cl34_log_pf - K38_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCl34_log_pf_dT9 - dK38_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  70.93993407850172 + -78.74345419369978 * tfactors.T9i + -68.2165 * tfactors.T913i + -2.18792 * tfactors.T913
                         + -1.25322 * tfactors.T9 + 0.140901 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  78.74345419369978 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -68.2165 * tfactors.T943i + (1.0/3.0) * -2.18792 * tfactors.T923i
                                  + -1.25322 + (5.0/3.0) * 0.140901 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_K39_to_He4_Cl35_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // K39 --> He4 + Cl35

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Cl35_log_pf, dCl35_log_pf_dT9;

    // interpolating Cl35 partition function
    get_partition_function_cached(Cl35, tfactors.T9, pf_cache, Cl35_log_pf, dCl35_log_pf_dT9);

    amrex::Real K39_log_pf, dK39_log_pf_dT9;

    // interpolating K39 partition function
    get_partition_function_cached(K39, tfactors.T9, pf_cache, K39_log_pf, dK39_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Cl35_log_pf - K39_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCl35_log_pf_dT9 - dK39_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  77.64945910439207 + -83.76815251770985 * tfactors.T9i + -68.2848 * tfactors.T913i + 0.0178545 * tfactors.T913
                         + -2.06783 * tfactors.T9 + 0.199659 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  83.76815251770985 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -68.2848 * tfactors.T943i + (1.0/3.0) * 0.0178545 * tfactors.T923i
                                  + -2.06783 + (5.0/3.0) * 0.199659 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ca40_to_p_K39_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ca40 --> p + K39

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Ca40_log_pf, dCa40_log_pf_dT9;

    // interpolating Ca40 partition function
    get_partition_function_cached(Ca40, tfactors.T9, pf_cache, Ca40_log_pf, dCa40_log_pf_dT9);

    amrex::Real K39_log_pf, dK39_log_pf_dT9;

    // interpolating K39 partition function
    get_partition_function_cached(K39, tfactors.T9, pf_cache, K39_log_pf, dK39_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + K39_log_pf - Ca40_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dK39_log_pf_dT9 - dCa40_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from lo18r
    ln_set_rate =  2786.453572485239 + -101.86683326322635 * tfactors.T9i + 802.18 * tfactors.T913i + -4010.27 * tfactors.T913
                         + 1137.69 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  101.86683326322635 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 802.18 * tfactors.T943i + (1.0/3.0) * -4010.27 * tfactors.T923i
                                  + 1137.69 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from lo18r
    ln_set_rate =  613.1725724852387 + -109.20919326322635 * tfactors.T9i + 641.844 * tfactors.T913i + -1248.49 * tfactors.T913
                         + 566.426 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  109.20919326322635 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 641.844 * tfactors.T943i + (1.0/3.0) * -1248.49 * tfactors.T923i
                                  + 566.426 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from lo18r
    ln_set_rate =  127.32557248523868 + -98.30957326322635 * tfactors.T9i + 41.1723 * tfactors.T913i + -149.299 * tfactors.T913
                         + 10.5229 * tfactors.T9 + -0.68208 * tfactors.T953 + 60.7367 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  98.30957326322635 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 41.1723 * tfactors.T943i + (1.0/3.0) * -149.299 * tfactors.T923i
                                  + 10.5229 + (5.0/3.0) * -0.68208 * tfactors.T923 + 60.7367 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ca40_to_He4_Ar36_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ca40 --> He4 + Ar36

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ca40_log_pf, dCa40_log_pf_dT9;

    // interpolating Ca40 partition function
    get_partition_function_cached(Ca40, tfactors.T9, pf_cache, Ca40_log_pf, dCa40_log_pf_dT9);

    amrex::Real Ar36_log_pf, dAr36_log_pf_dT9;

    // interpolating Ar36 partition function
    get_partition_function_cached(Ar36, tfactors.T9, pf_cache, Ar36_log_pf, dAr36_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Ar36_log_pf - Ca40_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dAr36_log_pf_dT9 - dCa40_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  77.28432125250924 + -81.69320665503751 * tfactors.T9i + -71.0046 * tfactors.T913i + 4.0656 * tfactors.T913
                         + -5.26509 * tfactors.T9 + 0.683546 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  81.69320665503751 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -71.0046 * tfactors.T943i + (1.0/3.0) * 4.0656 * tfactors.T923i
                                  + -5.26509 + (5.0/3.0) * 0.683546 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ca41_to_He4_Ar37_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ca41 --> He4 + Ar37

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Ar37_log_pf, dAr37_log_pf_dT9;

    // interpolating Ar37 partition function
    get_partition_function_cached(Ar37, tfactors.T9, pf_cache, Ar37_log_pf, dAr37_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ca41_log_pf, dCa41_log_pf_dT9;

    // interpolating Ca41 partition function
    get_partition_function_cached(Ca41, tfactors.T9, pf_cache, Ca41_log_pf, dCa41_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Ar37_log_pf - Ca41_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dAr37_log_pf_dT9 - dCa41_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  76.99880461067016 + -76.76558012507475 * tfactors.T9i + -71.0688 * tfactors.T913i + 2.71847 * tfactors.T913
                         + -3.17333 * tfactors.T9 + 0.335323 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  76.76558012507475 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -71.0688 * tfactors.T943i + (1.0/3.0) * 2.71847 * tfactors.T923i
                                  + -3.17333 + (5.0/3.0) * 0.335323 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Sc42_to_p_Ca41_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Sc42 --> p + Ca41

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Sc42_log_pf, dSc42_log_pf_dT9;

    // interpolating Sc42 partition function
    get_partition_function_cached(Sc42, tfactors.T9, pf_cache, Sc42_log_pf, dSc42_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Ca41_log_pf, dCa41_log_pf_dT9;

    // interpolating Ca41 partition function
    get_partition_function_cached(Ca41, tfactors.T9, pf_cache, Ca41_log_pf, dCa41_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Ca41_log_pf - Sc42_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dCa41_log_pf_dT9 - dSc42_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  65.08329070961724 + -49.5752100335526 * tfactors.T9i + -31.1321 * tfactors.T913i + -6.33393 * tfactors.T913
                         + -1.39316 * tfactors.T9 + 0.255024 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  49.5752100335526 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -31.1321 * tfactors.T943i + (1.0/3.0) * -6.33393 * tfactors.T923i
                                  + -1.39316 + (5.0/3.0) * 0.255024 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Sc42_to_He4_K38_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Sc42 --> He4 + K38

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Sc42_log_pf, dSc42_log_pf_dT9;

    // interpolating Sc42 partition function
    get_partition_function_cached(Sc42, tfactors.T9, pf_cache, Sc42_log_pf, dSc42_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real K38_log_pf, dK38_log_pf_dT9;

    // interpolating K38 partition function
    get_partition_function_cached(K38, tfactors.T9, pf_cache, K38_log_pf, dK38_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + K38_log_pf - Sc42_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dK38_log_pf_dT9 - dSc42_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  80.14878110581267 + -66.66976540451253 * tfactors.T9i + -73.7407 * tfactors.T913i + -7.08848 * tfactors.T913
                         + 0.322473 * tfactors.T9 + -0.0523557 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  66.66976540451253 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -73.7407 * tfactors.T943i + (1.0/3.0) * -7.08848 * tfactors.T923i
                                  + 0.322473 + (5.0/3.0) * -0.0523557 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Sc43_to_He4_K39_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Sc43 --> He4 + K39

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real K39_log_pf, dK39_log_pf_dT9;

    // interpolating K39 partition function
    get_partition_function_cached(K39, tfactors.T9, pf_cache, K39_log_pf, dK39_log_pf_dT9);

    amrex::Real Sc43_log_pf, dSc43_log_pf_dT9;

    // interpolating Sc43 partition function
    get_partition_function_cached(Sc43, tfactors.T9, pf_cache, Sc43_log_pf, dSc43_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + K39_log_pf - Sc43_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dK39_log_pf_dT9 - dSc43_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  78.37420492735454 + -55.76923537484507 * tfactors.T9i + -73.8006 * tfactors.T913i + 1.87885 * tfactors.T913
                         + -2.75862 * tfactors.T9 + 0.279678 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  55.76923537484507 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -73.8006 * tfactors.T943i + (1.0/3.0) * 1.87885 * tfactors.T923i
                                  + -2.75862 + (5.0/3.0) * 0.279678 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ti44_to_p_Sc43_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ti44 --> p + Sc43

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Ti44_log_pf, dTi44_log_pf_dT9;

    // interpolating Ti44 partition function
    get_partition_function_cached(Ti44, tfactors.T9, pf_cache, Ti44_log_pf, dTi44_log_pf_dT9);

    amrex::Real Sc43_log_pf, dSc43_log_pf_dT9;

    // interpolating Sc43 partition function
    get_partition_function_cached(Sc43, tfactors.T9, pf_cache, Sc43_log_pf, dSc43_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Sc43_log_pf - Ti44_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dSc43_log_pf_dT9 - dTi44_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  62.61343220623275 + -100.37294370395585 * tfactors.T9i + -32.1734 * tfactors.T913i + -1.77078 * tfactors.T913
                         + -2.21706 * tfactors.T9 + 0.298499 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  100.37294370395585 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -32.1734 * tfactors.T943i + (1.0/3.0) * -1.77078 * tfactors.T923i
                                  + -2.21706 + (5.0/3.0) * 0.298499 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ti44_to_He4_Ca40_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ti44 --> He4 + Ca40

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ca40_log_pf, dCa40_log_pf_dT9;

    // interpolating Ca40 partition function
    get_partition_function_cached(Ca40, tfactors.T9, pf_cache, Ca40_log_pf, dCa40_log_pf_dT9);

    amrex::Real Ti44_log_pf, dTi44_log_pf_dT9;

    // interpolating Ti44 partition function
    get_partition_function_cached(Ti44, tfactors.T9, pf_cache, Ti44_log_pf, dTi44_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Ca40_log_pf - Ti44_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCa40_log_pf_dT9 - dTi44_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from chw0
    ln_set_rate =  78.7006646483486 + -59.49768581557457 * tfactors.T9i + -76.4273 * tfactors.T913i + 3.87451 * tfactors.T913
                         + -3.61477 * tfactors.T9 + 0.367451 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  59.49768581557457 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -76.4273 * tfactors.T943i + (1.0/3.0) * 3.87451 * tfactors.T923i
                                  + -3.61477 + (5.0/3.0) * 0.367451 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ti45_to_He4_Ca41_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ti45 --> He4 + Ca41

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ca41_log_pf, dCa41_log_pf_dT9;

    // interpolating Ca41 partition function
    get_partition_function_cached(Ca41, tfactors.T9, pf_cache, Ca41_log_pf, dCa41_log_pf_dT9);

    amrex::Real Ti45_log_pf, dTi45_log_pf_dT9;

    // interpolating Ti45 partition function
    get_partition_function_cached(Ti45, tfactors.T9, pf_cache, Ti45_log_pf, dTi45_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Ca41_log_pf - Ti45_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCa41_log_pf_dT9 - dTi45_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  79.58456857742877 + -73.07720008533106 * tfactors.T9i + -76.4839 * tfactors.T913i + 3.03748 * tfactors.T913
                         + -2.59814 * tfactors.T9 + 0.210582 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  73.07720008533106 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -76.4839 * tfactors.T943i + (1.0/3.0) * 3.03748 * tfactors.T923i
                                  + -2.59814 + (5.0/3.0) * 0.210582 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_V46_to_p_Ti45_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // V46 --> p + Ti45

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real V46_log_pf, dV46_log_pf_dT9;

    // interpolating V46 partition function
    get_partition_function_cached(V46, tfactors.T9, pf_cache, V46_log_pf, dV46_log_pf_dT9);

    amrex::Real Ti45_log_pf, dTi45_log_pf_dT9;

    // interpolating Ti45 partition function
    get_partition_function_cached(Ti45, tfactors.T9, pf_cache, Ti45_log_pf, dTi45_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Ti45_log_pf - V46_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dTi45_log_pf_dT9 - dV46_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  63.21835285453139 + -62.13721694532568 * tfactors.T9i + -33.1981 * tfactors.T913i + -1.66837 * tfactors.T913
                         + -2.50285 * tfactors.T9 + 0.349152 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  62.13721694532568 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -33.1981 * tfactors.T943i + (1.0/3.0) * -1.66837 * tfactors.T923i
                                  + -2.50285 + (5.0/3.0) * 0.349152 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_V46_to_He4_Sc42_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // V46 --> He4 + Sc42

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Sc42_log_pf, dSc42_log_pf_dT9;

    // interpolating Sc42 partition function
    get_partition_function_cached(Sc42, tfactors.T9, pf_cache, Sc42_log_pf, dSc42_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real V46_log_pf, dV46_log_pf_dT9;

    // interpolating V46 partition function
    get_partition_function_cached(V46, tfactors.T9, pf_cache, V46_log_pf, dV46_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Sc42_log_pf - V46_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dSc42_log_pf_dT9 - dV46_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  76.73003072234292 + -85.63920699713294 * tfactors.T9i + -79.0687 * tfactors.T913i + -4.25493 * tfactors.T913
                         + 0.184703 * tfactors.T9 + -0.0785452 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  85.63920699713294 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -79.0687 * tfactors.T943i + (1.0/3.0) * -4.25493 * tfactors.T923i
                                  + 0.184703 + (5.0/3.0) * -0.0785452 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_V47_to_He4_Sc43_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // V47 --> He4 + Sc43

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real V47_log_pf, dV47_log_pf_dT9;

    // interpolating V47 partition function
    get_partition_function_cached(V47, tfactors.T9, pf_cache, V47_log_pf, dV47_log_pf_dT9);

    amrex::Real Sc43_log_pf, dSc43_log_pf_dT9;

    // interpolating Sc43 partition function
    get_partition_function_cached(Sc43, tfactors.T9, pf_cache, Sc43_log_pf, dSc43_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Sc43_log_pf - V47_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dSc43_log_pf_dT9 - dV47_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  84.67304382470803 + -95.66632297040007 * tfactors.T9i + -79.122 * tfactors.T913i + -7.07006 * tfactors.T913
                         + 0.424183 * tfactors.T9 + -0.0665231 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  95.66632297040007 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -79.122 * tfactors.T943i + (1.0/3.0) * -7.07006 * tfactors.T923i
                                  + 0.424183 + (5.0/3.0) * -0.0665231 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Cr48_to_p_V47_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Cr48 --> p + V47

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Cr48_log_pf, dCr48_log_pf_dT9;

    // interpolating Cr48 partition function
    get_partition_function_cached(Cr48, tfactors.T9, pf_cache, Cr48_log_pf, dCr48_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real V47_log_pf, dV47_log_pf_dT9;

    // interpolating V47 partition function
    get_partition_function_cached(V47, tfactors.T9, pf_cache, V47_log_pf, dV47_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + V47_log_pf - Cr48_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dV47_log_pf_dT9 - dCr48_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from nfisn
    ln_set_rate =  67.75975303984967 + -100.08956223422913 * tfactors.T9i + -34.0548 * tfactors.T913i + -3.41973 * tfactors.T913
                         + 1.16501 * tfactors.T9 + -0.105543 * tfactors.T953 + -6.20886 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  100.08956223422913 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -34.0548 * tfactors.T943i + (1.0/3.0) * -3.41973 * tfactors.T923i
                                  + 1.16501 + (5.0/3.0) * -0.105543 * tfactors.T923 + -6.20886 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from nfisn
    ln_set_rate =  536.5429530398497 + -99.32517223422914 * tfactors.T9i + 317.171 * tfactors.T913i + -911.679 * tfactors.T913
                         + 94.4245 * tfactors.T9 + -10.1973 * tfactors.T953 + 332.227 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  99.32517223422914 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 317.171 * tfactors.T943i + (1.0/3.0) * -911.679 * tfactors.T923i
                                  + 94.4245 + (5.0/3.0) * -10.1973 * tfactors.T923 + 332.227 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from nfisn
    ln_set_rate =  48.911453039849675 + -93.78359723422915 * tfactors.T9i + -45.9868 * tfactors.T913i + 13.6822 * tfactors.T913
                         + -0.376902 * tfactors.T9 + -0.0194875 * tfactors.T953 + -6.9232499999999995 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  93.78359723422915 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -45.9868 * tfactors.T943i + (1.0/3.0) * 13.6822 * tfactors.T923i
                                  + -0.376902 + (5.0/3.0) * -0.0194875 * tfactors.T923 + -6.9232499999999995 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from nfisn
    ln_set_rate =  65.64255303984967 + -94.54467623422914 * tfactors.T9i + -110.655 * tfactors.T913i + 83.0232 * tfactors.T913
                         + -19.7762 * tfactors.T9 + 3.03961 * tfactors.T953 + -47.9742 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  94.54467623422914 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -110.655 * tfactors.T943i + (1.0/3.0) * 83.0232 * tfactors.T923i
                                  + -19.7762 + (5.0/3.0) * 3.03961 * tfactors.T923 + -47.9742 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Cr48_to_He4_Ti44_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Cr48 --> He4 + Ti44

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Cr48_log_pf, dCr48_log_pf_dT9;

    // interpolating Cr48 partition function
    get_partition_function_cached(Cr48, tfactors.T9, pf_cache, Cr48_log_pf, dCr48_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ti44_log_pf, dTi44_log_pf_dT9;

    // interpolating Ti44 partition function
    get_partition_function_cached(Ti44, tfactors.T9, pf_cache, Ti44_log_pf, dTi44_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Ti44_log_pf - Cr48_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dTi44_log_pf_dT9 - dCr48_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  89.75906465832496 + -89.32364150067336 * tfactors.T9i + -81.667 * tfactors.T913i + -10.6333 * tfactors.T913
                         + -0.672613 * tfactors.T9 + 0.161209 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  89.32364150067336 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -81.667 * tfactors.T943i + (1.0/3.0) * -10.6333 * tfactors.T923i
                                  + -0.672613 + (5.0/3.0) * 0.161209 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Cr49_to_He4_Ti45_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Cr49 --> He4 + Ti45

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Cr49_log_pf, dCr49_log_pf_dT9;

    // interpolating Cr49 partition function
    get_partition_function_cached(Cr49, tfactors.T9, pf_cache, Cr49_log_pf, dCr49_log_pf_dT9);

    amrex::Real Ti45_log_pf, dTi45_log_pf_dT9;

    // interpolating Ti45 partition function
    get_partition_function_cached(Ti45, tfactors.T9, pf_cache, Ti45_log_pf, dTi45_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Ti45_log_pf - Cr49_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dTi45_log_pf_dT9 - dCr49_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  89.2767795487736 + -101.50490417292947 * tfactors.T9i + -81.7175 * tfactors.T913i + -10.1755 * tfactors.T913
                         + 0.364167 * tfactors.T9 + -0.000317666 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  101.50490417292947 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -81.7175 * tfactors.T943i + (1.0/3.0) * -10.1755 * tfactors.T923i
                                  + 0.364167 + (5.0/3.0) * -0.000317666 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Mn50_to_p_Cr49_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mn50 --> p + Cr49

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Mn50_log_pf, dMn50_log_pf_dT9;

    // interpolating Mn50 partition function
    get_partition_function_cached(Mn50, tfactors.T9, pf_cache, Mn50_log_pf, dMn50_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Cr49_log_pf, dCr49_log_pf_dT9;

    // interpolating Cr49 partition function
    get_partition_function_cached(Cr49, tfactors.T9, pf_cache, Cr49_log_pf, dCr49_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Cr49_log_pf - Mn50_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dCr49_log_pf_dT9 - dMn50_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  61.87490359919837 + -53.188856976553964 * tfactors.T9i + -35.2018 * tfactors.T913i + 0.168579 * tfactors.T913
                         + -2.87983 * tfactors.T9 + 0.378768 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  53.188856976553964 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -35.2018 * tfactors.T943i + (1.0/3.0) * 0.168579 * tfactors.T923i
                                  + -2.87983 + (5.0/3.0) * 0.378768 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Mn50_to_He4_V46_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mn50 --> He4 + V46

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Mn50_log_pf, dMn50_log_pf_dT9;

    // interpolating Mn50 partition function
    get_partition_function_cached(Mn50, tfactors.T9, pf_cache, Mn50_log_pf, dMn50_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real V46_log_pf, dV46_log_pf_dT9;

    // interpolating V46 partition function
    get_partition_function_cached(V46, tfactors.T9, pf_cache, V46_log_pf, dV46_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + V46_log_pf - Mn50_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dV46_log_pf_dT9 - dMn50_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  79.2797302934406 + -92.55654420415776 * tfactors.T9i + -84.2255 * tfactors.T913i + -4.85634 * tfactors.T913
                         + 0.0528515 * tfactors.T9 + -0.0425496 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  92.55654420415776 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -84.2255 * tfactors.T943i + (1.0/3.0) * -4.85634 * tfactors.T923i
                                  + 0.0528515 + (5.0/3.0) * -0.0425496 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Mn51_to_He4_V47_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mn51 --> He4 + V47

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Mn51_log_pf, dMn51_log_pf_dT9;

    // interpolating Mn51 partition function
    get_partition_function_cached(Mn51, tfactors.T9, pf_cache, Mn51_log_pf, dMn51_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real V47_log_pf, dV47_log_pf_dT9;

    // interpolating V47 partition function
    get_partition_function_cached(V47, tfactors.T9, pf_cache, V47_log_pf, dV47_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + V47_log_pf - Mn51_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dV47_log_pf_dT9 - dMn51_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  81.42763571008429 + -100.50726374990671 * tfactors.T9i + -84.2732 * tfactors.T913i + -2.98061 * tfactors.T913
                         + -0.531361 * tfactors.T9 + 0.023612 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  100.50726374990671 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -84.2732 * tfactors.T943i + (1.0/3.0) * -2.98061 * tfactors.T923i
                                  + -0.531361 + (5.0/3.0) * 0.023612 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Fe52_to_p_Mn51_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe52 --> p + Mn51

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Mn51_log_pf, dMn51_log_pf_dT9;

    // interpolating Mn51 partition function
    get_partition_function_cached(Mn51, tfactors.T9, pf_cache, Mn51_log_pf, dMn51_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Fe52_log_pf, dFe52_log_pf_dT9;

    // interpolating Fe52 partition function
    get_partition_function_cached(Fe52, tfactors.T9, pf_cache, Fe52_log_pf, dFe52_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Mn51_log_pf - Fe52_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dMn51_log_pf_dT9 - dFe52_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  61.74743132228039 + -85.61663846070292 * tfactors.T9i + -36.1825 * tfactors.T913i + 0.873042 * tfactors.T913
                         + -2.89731 * tfactors.T9 + 0.364394 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  85.61663846070292 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -36.1825 * tfactors.T943i + (1.0/3.0) * 0.873042 * tfactors.T923i
                                  + -2.89731 + (5.0/3.0) * 0.364394 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Fe52_to_He4_Cr48_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe52 --> He4 + Cr48

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Cr48_log_pf, dCr48_log_pf_dT9;

    // interpolating Cr48 partition function
    get_partition_function_cached(Cr48, tfactors.T9, pf_cache, Cr48_log_pf, dCr48_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Fe52_log_pf, dFe52_log_pf_dT9;

    // interpolating Fe52 partition function
    get_partition_function_cached(Fe52, tfactors.T9, pf_cache, Fe52_log_pf, dFe52_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Cr48_log_pf - Fe52_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCr48_log_pf_dT9 - dFe52_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  90.149113992515 + -92.09363997638769 * tfactors.T9i + -86.7459 * tfactors.T913i + -9.79373 * tfactors.T913
                         + -0.772169 * tfactors.T9 + 0.155883 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  92.09363997638769 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -86.7459 * tfactors.T943i + (1.0/3.0) * -9.79373 * tfactors.T923i
                                  + -0.772169 + (5.0/3.0) * 0.155883 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Co55_to_He4_Mn51_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co55 --> He4 + Mn51

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Mn51_log_pf, dMn51_log_pf_dT9;

    // interpolating Mn51 partition function
    get_partition_function_cached(Mn51, tfactors.T9, pf_cache, Mn51_log_pf, dMn51_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Co55_log_pf, dCo55_log_pf_dT9;

    // interpolating Co55 partition function
    get_partition_function_cached(Co55, tfactors.T9, pf_cache, Co55_log_pf, dCo55_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Mn51_log_pf - Co55_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dMn51_log_pf_dT9 - dCo55_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  90.61473213109255 + -95.29300562248025 * tfactors.T9i + -89.274 * tfactors.T913i + -10.4373 * tfactors.T913
                         + 1.00492 * tfactors.T9 + -0.125548 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  95.29300562248025 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -89.274 * tfactors.T943i + (1.0/3.0) * -10.4373 * tfactors.T923i
                                  + 1.00492 + (5.0/3.0) * -0.125548 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ni56_to_p_Co55_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni56 --> p + Co55

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Ni56_log_pf, dNi56_log_pf_dT9;

    // interpolating Ni56 partition function
    get_partition_function_cached(Ni56, tfactors.T9, pf_cache, Ni56_log_pf, dNi56_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Co55_log_pf, dCo55_log_pf_dT9;

    // interpolating Co55 partition function
    get_partition_function_cached(Co55, tfactors.T9, pf_cache, Co55_log_pf, dCo55_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Co55_log_pf - Ni56_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dCo55_log_pf_dT9 - dNi56_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  63.15120404192561 + -83.16460378149378 * tfactors.T9i + -38.1053 * tfactors.T913i + -0.210947 * tfactors.T913
                         + -2.68377 * tfactors.T9 + 0.355814 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  83.16460378149378 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -38.1053 * tfactors.T943i + (1.0/3.0) * -0.210947 * tfactors.T923i
                                  + -2.68377 + (5.0/3.0) * 0.355814 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ni56_to_He4_Fe52_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni56 --> He4 + Fe52

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Ni56_log_pf, dNi56_log_pf_dT9;

    // interpolating Ni56 partition function
    get_partition_function_cached(Ni56, tfactors.T9, pf_cache, Ni56_log_pf, dNi56_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Fe52_log_pf, dFe52_log_pf_dT9;

    // interpolating Fe52 partition function
    get_partition_function_cached(Fe52, tfactors.T9, pf_cache, Fe52_log_pf, dFe52_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Fe52_log_pf - Ni56_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dFe52_log_pf_dT9 - dNi56_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  91.62430485073777 + -92.84097094326391 * tfactors.T9i + -91.6819 * tfactors.T913i + -9.51885 * tfactors.T913
                         + -0.533014 * tfactors.T9 + 0.0892607 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  92.84097094326391 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -91.6819 * tfactors.T943i + (1.0/3.0) * -9.51885 * tfactors.T923i
                                  + -0.533014 + (5.0/3.0) * 0.0892607 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_C12_to_He4_He4_He4_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // C12 --> 3 He4

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real C12_log_pf, dC12_log_pf_dT9;

    // setting C12 log(partition function) to 0.0 by default, independent of T
    C12_log_pf = 0.0_rt;
    dC12_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = He4_log_pf + He4_log_pf + He4_log_pf - C12_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dHe4_log_pf_dT9 + dHe4_log_pf_dT9 - dC12_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from fy05r
    ln_set_rate =  34.96090397991297 + -85.44440046993657 * tfactors.T9i + -23.57 * tfactors.T913i + 20.4886 * tfactors.T913
                         + -12.9882 * tfactors.T9 + -20.0 * tfactors.T953 + 0.8333300000000001 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  85.44440046993657 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -23.57 * tfactors.T943i + (1.0/3.0) * 20.4886 * tfactors.T923i
                                  + -12.9882 + (5.0/3.0) * -20.0 * tfactors.T923 + 0.8333300000000001 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from fy05n
    ln_set_rate =  45.77825197991297 + -84.41994046993656 * tfactors.T9i + -37.06 * tfactors.T913i + 29.3493 * tfactors.T913
                         + -115.507 * tfactors.T9 + -10.0 * tfactors.T953 + 1.66667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  84.41994046993656 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -37.06 * tfactors.T943i + (1.0/3.0) * 29.3493 * tfactors.T923i
                                  + -115.507 + (5.0/3.0) * -10.0 * tfactors.T923 + 1.66667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from fy05r
    ln_set_rate =  22.398803979912973 + -88.54650046993656 * tfactors.T9i + -13.49 * tfactors.T913i + 21.4259 * tfactors.T913
                         + -1.34769 * tfactors.T9 + 0.0879816 * tfactors.T953 + -10.1653 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  88.54650046993656 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -13.49 * tfactors.T943i + (1.0/3.0) * 21.4259 * tfactors.T923i
                                  + -1.34769 + (5.0/3.0) * 0.0879816 * tfactors.T923 + -10.1653 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_C12_to_p_N15_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // C12 + He4 --> p + N15

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real N15_log_pf, dN15_log_pf_dT9;

    // setting N15 log(partition function) to 0.0 by default, independent of T
    N15_log_pf = 0.0_rt;
    dN15_log_pf_dT9 = 0.0_rt;

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real C12_log_pf, dC12_log_pf_dT9;

    // setting C12 log(partition function) to 0.0 by default, independent of T
    C12_log_pf = 0.0_rt;
    dC12_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = p_log_pf + N15_log_pf - He4_log_pf - C12_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dN15_log_pf_dT9 - dHe4_log_pf_dT9 - dC12_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from nacrn
    ln_set_rate =  27.135846229234243 + -57.62215691264642 * tfactors.T9i + -15.253 * tfactors.T913i + 1.59318 * tfactors.T913
                         + 2.4479 * tfactors.T9 + -2.19708 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  57.62215691264642 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -15.253 * tfactors.T943i + (1.0/3.0) * 1.59318 * tfactors.T923i
                                  + 2.4479 + (5.0/3.0) * -2.19708 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from nacrr
    ln_set_rate =  -6.915773770765759 + -58.78595691264642 * tfactors.T9i + 22.7105 * tfactors.T913
                         + -2.90707 * tfactors.T9 + 0.205754 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  58.78595691264642 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 22.7105 * tfactors.T923i
                                  + -2.90707 + (5.0/3.0) * 0.205754 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from nacrr
    ln_set_rate =  20.556646229234243 + -65.02815691264642 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  65.02815691264642 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from nacrr
    ln_set_rate =  -5.21402377076576 + -59.643326912646415 * tfactors.T9i + 30.8497 * tfactors.T913
                         + -8.50433 * tfactors.T9 + -1.54426 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  59.643326912646415 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 30.8497 * tfactors.T923i
                                  + -8.50433 + (5.0/3.0) * -1.54426 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_N14_to_p_O17_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // N14 + He4 --> p + O17

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real O17_log_pf, dO17_log_pf_dT9;

    // interpolating O17 partition function
    get_partition_function_cached(O17, tfactors.T9, pf_cache, O17_log_pf, dO17_log_pf_dT9);

    amrex::Real N14_log_pf, dN14_log_pf_dT9;

    // setting N14 log(partition function) to 0.0 by default, independent of T
    N14_log_pf = 0.0_rt;
    dN14_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = p_log_pf + O17_log_pf - He4_log_pf - N14_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dO17_log_pf_dT9 - dHe4_log_pf_dT9 - dN14_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from il10r
    ln_set_rate =  -7.592359780998223 + -14.584520682879308 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  14.584520682879308 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10n
    ln_set_rate =  19.194270219001776 + -13.831125682879309 * tfactors.T9i + -16.9078 * tfactors.T913i
                         + -2.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  13.831125682879309 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -16.9078 * tfactors.T943i
                                  + (5.0/3.0) * -2.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  9.789270219001777 + -18.78977568287931 * tfactors.T9i + 5.10182 * tfactors.T913
                         + 0.379373 * tfactors.T9 + -0.0672515 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  18.78977568287931 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 5.10182 * tfactors.T923i
                                  + 0.379373 + (5.0/3.0) * -0.0672515 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  5.148870219001777 + -15.945895682879309 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  15.945895682879309 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_N15_to_p_O18_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // N15 + He4 --> p + O18

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real O18_log_pf, dO18_log_pf_dT9;

    // interpolating O18 partition function
    get_partition_function_cached(O18, tfactors.T9, pf_cache, O18_log_pf, dO18_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real N15_log_pf, dN15_log_pf_dT9;

    // setting N15 log(partition function) to 0.0 by default, independent of T
    N15_log_pf = 0.0_rt;
    dN15_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = p_log_pf + O18_log_pf - He4_log_pf - N15_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dO18_log_pf_dT9 - dHe4_log_pf_dT9 - dN15_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from il10r
    ln_set_rate =  -29.6926341462146 + -46.42955443414268 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  46.42955443414268 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10n
    ln_set_rate =  25.178865853785396 + -46.183670434142684 * tfactors.T9i + -16.6979 * tfactors.T913i
                         + -3.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  46.183670434142684 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -16.6979 * tfactors.T943i
                                  + (5.0/3.0) * -3.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  7.155285853785398 + -51.50702043414268 * tfactors.T9i + 11.6568 * tfactors.T913
                         + -2.16303 * tfactors.T9 + 0.209965 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  51.50702043414268 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 11.6568 * tfactors.T923i
                                  + -2.16303 + (5.0/3.0) * 0.209965 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  8.4842658537854 + -47.84667043414268 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  47.84667043414268 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_O16_to_He4_N13_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // O16 + p --> He4 + N13

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real O16_log_pf, dO16_log_pf_dT9;

    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors.T9, pf_cache, O16_log_pf, dO16_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real N13_log_pf, dN13_log_pf_dT9;

    // setting N13 log(partition function) to 0.0 by default, independent of T
    N13_log_pf = 0.0_rt;
    dN13_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = He4_log_pf + N13_log_pf - p_log_pf - O16_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dN13_log_pf_dT9 - dp_log_pf_dT9 - dO16_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from cf88n
    ln_set_rate =  42.21642061342456 + -60.557329596678294 * tfactors.T9i + -35.829 * tfactors.T913i + -0.530275 * tfactors.T913
                         + -0.982462 * tfactors.T9 + 0.0808059 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  60.557329596678294 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -35.829 * tfactors.T943i + (1.0/3.0) * -0.530275 * tfactors.T923i
                                  + -0.982462 + (5.0/3.0) * 0.0808059 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_O16_to_p_F19_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // O16 + He4 --> p + F19

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real O16_log_pf, dO16_log_pf_dT9;

    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors.T9, pf_cache, O16_log_pf, dO16_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real F19_log_pf, dF19_log_pf_dT9;

    // interpolating F19 partition function
    get_partition_function_cached(F19, tfactors.T9, pf_cache, F19_log_pf, dF19_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = p_log_pf + F19_log_pf - He4_log_pf - O16_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dF19_log_pf_dT9 - dHe4_log_pf_dT9 - dO16_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from nacr
    ln_set_rate =  -53.121227264044336 + -94.28220973651266 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  94.28220973651266 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from nacr
    ln_set_rate =  25.874672735955667 + -94.15455973651265 * tfactors.T9i + -18.116 * tfactors.T913i
                         + 1.86674 * tfactors.T9 + -7.5666 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  94.15455973651265 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -18.116 * tfactors.T943i
                                  + 1.86674 + (5.0/3.0) * -7.5666 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from nacrr
    ln_set_rate =  13.941672735955667 + -97.44055973651265 * tfactors.T9i
                         + -0.21103 * tfactors.T9 + 2.87702 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  97.44055973651265 * tfactors.T9i * tfactors.T9i
                                  + -0.21103 + 2.87702 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from nacr
    ln_set_rate =  14.778572735955667 + -97.90640973651266 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  97.90640973651266 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from nacr
    ln_set_rate =  7.822072735955667 + -96.62283973651265 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  96.62283973651265 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Ne20_to_p_Na23_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ne20 + He4 --> p + Na23

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Na23_log_pf, dNa23_log_pf_dT9;

    // interpolating Na23 partition function
    get_partition_function_cached(Na23, tfactors.T9, pf_cache, Na23_log_pf, dNa23_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ne20_log_pf, dNe20_log_pf_dT9;

    // interpolating Ne20 partition function
    get_partition_function_cached(Ne20, tfactors.T9, pf_cache, Ne20_log_pf, dNe20_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Na23_log_pf - He4_log_pf - Ne20_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dNa23_log_pf_dT9 - dHe4_log_pf_dT9 - dNe20_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from il10r
    ln_set_rate =  -6.360024278486355 + -29.889658832139027 * tfactors.T9i + 19.7297 * tfactors.T913
                         + -2.20987 * tfactors.T9 + 0.153374 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  29.889658832139027 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 19.7297 * tfactors.T923i
                                  + -2.20987 + (5.0/3.0) * 0.153374 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  0.24516522151364528 + -29.434918832139026 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  29.434918832139026 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10n
    ln_set_rate =  19.202935721513644 + -27.573888832139026 * tfactors.T9i + -20.0024 * tfactors.T913i + 11.5988 * tfactors.T913
                         + -1.37398 * tfactors.T9 + -1.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  27.573888832139026 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -20.0024 * tfactors.T943i + (1.0/3.0) * 11.5988 * tfactors.T923i
                                  + -1.37398 + (5.0/3.0) * -1.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Ne20_to_C12_C12_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ne20 + He4 --> C12 + C12

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real C12_log_pf, dC12_log_pf_dT9;

    // setting C12 log(partition function) to 0.0 by default, independent of T
    C12_log_pf = 0.0_rt;
    dC12_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ne20_log_pf, dNe20_log_pf_dT9;

    // interpolating Ne20 partition function
    get_partition_function_cached(Ne20, tfactors.T9, pf_cache, Ne20_log_pf, dNe20_log_pf_dT9);

    amrex::Real net_log_pf = C12_log_pf + C12_log_pf - He4_log_pf - Ne20_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dC12_log_pf_dT9 + dC12_log_pf_dT9 - dHe4_log_pf_dT9 - dNe20_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from cf88r
    ln_set_rate =  61.474151468919175 + -53.57824966896331 * tfactors.T9i + -84.165 * tfactors.T913i + -1.56627 * tfactors.T913
                         + -0.0736084 * tfactors.T9 + -0.072797 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  53.57824966896331 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -84.165 * tfactors.T943i + (1.0/3.0) * -1.56627 * tfactors.T923i
                                  + -0.0736084 + (5.0/3.0) * -0.072797 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Ne21_to_He4_F18_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ne21 + p --> He4 + F18

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Ne21_log_pf, dNe21_log_pf_dT9;

    // interpolating Ne21 partition function
    get_partition_function_cached(Ne21, tfactors.T9, pf_cache, Ne21_log_pf, dNe21_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real F18_log_pf, dF18_log_pf_dT9;

    // interpolating F18 partition function
    get_partition_function_cached(F18, tfactors.T9, pf_cache, F18_log_pf, dF18_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + F18_log_pf - p_log_pf - Ne21_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dF18_log_pf_dT9 - dp_log_pf_dT9 - dNe21_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from rpsmr
    ln_set_rate =  50.63668915976785 + -22.04702300768339 * tfactors.T9i + 21.4461 * tfactors.T913i + -73.252 * tfactors.T913
                         + 2.42329 * tfactors.T9 + -0.077278 * tfactors.T953 + 40.7604 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  22.04702300768339 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 21.4461 * tfactors.T943i + (1.0/3.0) * -73.252 * tfactors.T923i
                                  + 2.42329 + (5.0/3.0) * -0.077278 * tfactors.T923 + 40.7604 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Ne22_to_He4_F19_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ne22 + p --> He4 + F19

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real F19_log_pf, dF19_log_pf_dT9;

    // interpolating F19 partition function
    get_partition_function_cached(F19, tfactors.T9, pf_cache, F19_log_pf, dF19_log_pf_dT9);

    amrex::Real Ne22_log_pf, dNe22_log_pf_dT9;

    // interpolating Ne22 partition function
    get_partition_function_cached(Ne22, tfactors.T9, pf_cache, Ne22_log_pf, dNe22_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + F19_log_pf - p_log_pf - Ne22_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dF19_log_pf_dT9 - dp_log_pf_dT9 - dNe22_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from da18r
    ln_set_rate =  54.774145434237155 + -22.26126198150658 * tfactors.T9i + -38.7722 * tfactors.T913i + -13.3654 * tfactors.T913
                         + 0.863648 * tfactors.T9 + -0.0451491 * tfactors.T953 + 1.33333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  22.26126198150658 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -38.7722 * tfactors.T943i + (1.0/3.0) * -13.3654 * tfactors.T923i
                                  + 0.863648 + (5.0/3.0) * -0.0451491 * tfactors.T923 + 1.33333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from da18r
    ln_set_rate =  53.51334543423716 + -65.19766198150657 * tfactors.T9i + -34.5008 * tfactors.T913i + 56.9316 * tfactors.T913
                         + 2.09613 * tfactors.T9 + -32.496 * tfactors.T953 + 0.333333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  65.19766198150657 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -34.5008 * tfactors.T943i + (1.0/3.0) * 56.9316 * tfactors.T923i
                                  + 2.09613 + (5.0/3.0) * -32.496 * tfactors.T923 + 0.333333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from da18r
    ln_set_rate =  29432.442445434237 + -152.44286198150658 * tfactors.T9i + 12625.1 * tfactors.T913i + -49107.1 * tfactors.T913
                         + 9227.53 * tfactors.T9 + -2086.65 * tfactors.T953 + 14520.2 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  152.44286198150658 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 12625.1 * tfactors.T943i + (1.0/3.0) * -49107.1 * tfactors.T923i
                                  + 9227.53 + (5.0/3.0) * -2086.65 * tfactors.T923 + 14520.2 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Na23_to_C12_C12_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Na23 + p --> C12 + C12

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real C12_log_pf, dC12_log_pf_dT9;

    // setting C12 log(partition function) to 0.0 by default, independent of T
    C12_log_pf = 0.0_rt;
    dC12_log_pf_dT9 = 0.0_rt;

    amrex::Real Na23_log_pf, dNa23_log_pf_dT9;

    // interpolating Na23 partition function
    get_partition_function_cached(Na23, tfactors.T9, pf_cache, Na23_log_pf, dNa23_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = C12_log_pf + C12_log_pf - p_log_pf - Na23_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dC12_log_pf_dT9 + dC12_log_pf_dT9 - dp_log_pf_dT9 - dNa23_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from cf88r
    ln_set_rate =  60.92541574740554 + -26.004360836807496 * tfactors.T9i + -84.165 * tfactors.T913i + -1.4191 * tfactors.T913
                         + -0.114619 * tfactors.T9 + -0.070307 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  26.004360836807496 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -84.165 * tfactors.T943i + (1.0/3.0) * -1.4191 * tfactors.T923i
                                  + -0.114619 + (5.0/3.0) * -0.070307 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Mg24_to_p_Al27_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mg24 + He4 --> p + Al27

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Al27_log_pf, dAl27_log_pf_dT9;

    // interpolating Al27 partition function
    get_partition_function_cached(Al27, tfactors.T9, pf_cache, Al27_log_pf, dAl27_log_pf_dT9);

    amrex::Real Mg24_log_pf, dMg24_log_pf_dT9;

    // interpolating Mg24 partition function
    get_partition_function_cached(Mg24, tfactors.T9, pf_cache, Mg24_log_pf, dMg24_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = p_log_pf + Al27_log_pf - He4_log_pf - Mg24_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dAl27_log_pf_dT9 - dHe4_log_pf_dT9 - dMg24_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from il10r
    ln_set_rate =  -6.4280417530584835 + -22.818701538287602 * tfactors.T9i + 18.0416 * tfactors.T913
                         + -1.54137 * tfactors.T9 + 0.0847506 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  22.818701538287602 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 18.0416 * tfactors.T923i
                                  + -1.54137 + (5.0/3.0) * 0.0847506 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  -26.268451753058486 + -19.5392135382876 * tfactors.T9i + 5.18642 * tfactors.T913i + -34.7936 * tfactors.T913
                         + 168.225 * tfactors.T9 + -115.825 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  19.5392135382876 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 5.18642 * tfactors.T943i + (1.0/3.0) * -34.7936 * tfactors.T923i
                                  + 168.225 + (5.0/3.0) * -115.825 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10n
    ln_set_rate =  30.057448246941515 + -18.576201538287602 * tfactors.T9i + -26.4162 * tfactors.T913i
                         + -2.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  18.576201538287602 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -26.4162 * tfactors.T943i
                                  + (5.0/3.0) * -2.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Mg24_to_C12_O16_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mg24 + He4 --> C12 + O16

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real C12_log_pf, dC12_log_pf_dT9;

    // setting C12 log(partition function) to 0.0 by default, independent of T
    C12_log_pf = 0.0_rt;
    dC12_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Mg24_log_pf, dMg24_log_pf_dT9;

    // interpolating Mg24 partition function
    get_partition_function_cached(Mg24, tfactors.T9, pf_cache, Mg24_log_pf, dMg24_log_pf_dT9);

    amrex::Real O16_log_pf, dO16_log_pf_dT9;

    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors.T9, pf_cache, O16_log_pf, dO16_log_pf_dT9);

    amrex::Real net_log_pf = C12_log_pf + O16_log_pf - He4_log_pf - Mg24_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dC12_log_pf_dT9 + dO16_log_pf_dT9 - dHe4_log_pf_dT9 - dMg24_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from cf88r
    ln_set_rate =  49.57295819231048 + -78.2098115311048 * tfactors.T9i + -133.413 * tfactors.T913i + 50.1572 * tfactors.T913
                         + -3.15987 * tfactors.T9 + 0.0178251 * tfactors.T953 + -23.7027 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  78.2098115311048 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -133.413 * tfactors.T943i + (1.0/3.0) * 50.1572 * tfactors.T923i
                                  + -3.15987 + (5.0/3.0) * 0.0178251 * tfactors.T923 + -23.7027 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Mg25_to_He4_Na22_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mg25 + p --> He4 + Na22

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Mg25_log_pf, dMg25_log_pf_dT9;

    // interpolating Mg25 partition function
    get_partition_function_cached(Mg25, tfactors.T9, pf_cache, Mg25_log_pf, dMg25_log_pf_dT9);

    amrex::Real Na22_log_pf, dNa22_log_pf_dT9;

    // interpolating Na22 partition function
    get_partition_function_cached(Na22, tfactors.T9, pf_cache, Na22_log_pf, dNa22_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Na22_log_pf - p_log_pf - Mg25_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dNa22_log_pf_dT9 - dp_log_pf_dT9 - dMg25_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  46.304619578366136 + -36.52330379082377 * tfactors.T9i + -50.0924 * tfactors.T913i + 0.807739 * tfactors.T913
                         + -0.956029 * tfactors.T9 + 0.0793321 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  36.52330379082377 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -50.0924 * tfactors.T943i + (1.0/3.0) * 0.807739 * tfactors.T923i
                                  + -0.956029 + (5.0/3.0) * 0.0793321 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Al27_to_C12_O16_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Al27 + p --> C12 + O16

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real C12_log_pf, dC12_log_pf_dT9;

    // setting C12 log(partition function) to 0.0 by default, independent of T
    C12_log_pf = 0.0_rt;
    dC12_log_pf_dT9 = 0.0_rt;

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Al27_log_pf, dAl27_log_pf_dT9;

    // interpolating Al27 partition function
    get_partition_function_cached(Al27, tfactors.T9, pf_cache, Al27_log_pf, dAl27_log_pf_dT9);

    amrex::Real O16_log_pf, dO16_log_pf_dT9;

    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors.T9, pf_cache, O16_log_pf, dO16_log_pf_dT9);

    amrex::Real net_log_pf = C12_log_pf + O16_log_pf - p_log_pf - Al27_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dC12_log_pf_dT9 + dO16_log_pf_dT9 - dp_log_pf_dT9 - dAl27_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from cf88r
    ln_set_rate =  68.96430994536895 + -59.8005159928172 * tfactors.T9i + -119.242 * tfactors.T913i + 13.3667 * tfactors.T913
                         + 0.295425 * tfactors.T9 + -0.267288 * tfactors.T953 + -9.91729 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  59.8005159928172 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -119.242 * tfactors.T943i + (1.0/3.0) * 13.3667 * tfactors.T923i
                                  + 0.295425 + (5.0/3.0) * -0.267288 * tfactors.T923 + -9.91729 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Si28_to_p_P31_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Si28 + He4 --> p + P31

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real P31_log_pf, dP31_log_pf_dT9;

    // interpolating P31 partition function
    get_partition_function_cached(P31, tfactors.T9, pf_cache, P31_log_pf, dP31_log_pf_dT9);

    amrex::Real Si28_log_pf, dSi28_log_pf_dT9;

    // interpolating Si28 partition function
    get_partition_function_cached(Si28, tfactors.T9, pf_cache, Si28_log_pf, dSi28_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + P31_log_pf - He4_log_pf - Si28_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dP31_log_pf_dT9 - dHe4_log_pf_dT9 - dSi28_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from il10r
    ln_set_rate =  -11.415669089666713 + -25.663579682396225 * tfactors.T9i + 21.521 * tfactors.T913
                         + -1.90355 * tfactors.T9 + 0.092724 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  25.663579682396225 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 21.521 * tfactors.T923i
                                  + -1.90355 + (5.0/3.0) * 0.092724 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  -13.441669089666712 + -24.114989682396224 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  24.114989682396224 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10n
    ln_set_rate =  60.360230910333286 + -22.237829682396224 * tfactors.T9i + -31.932 * tfactors.T913i + -77.0334 * tfactors.T913
                         + -43.6847 * tfactors.T9 + -4.28955 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  22.237829682396224 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -31.932 * tfactors.T943i + (1.0/3.0) * -77.0334 * tfactors.T923i
                                  + -43.6847 + (5.0/3.0) * -4.28955 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Si28_to_C12_Ne20_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Si28 + He4 --> C12 + Ne20

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real C12_log_pf, dC12_log_pf_dT9;

    // setting C12 log(partition function) to 0.0 by default, independent of T
    C12_log_pf = 0.0_rt;
    dC12_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Si28_log_pf, dSi28_log_pf_dT9;

    // interpolating Si28 partition function
    get_partition_function_cached(Si28, tfactors.T9, pf_cache, Si28_log_pf, dSi28_log_pf_dT9);

    amrex::Real Ne20_log_pf, dNe20_log_pf_dT9;

    // interpolating Ne20 partition function
    get_partition_function_cached(Ne20, tfactors.T9, pf_cache, Ne20_log_pf, dNe20_log_pf_dT9);

    amrex::Real net_log_pf = C12_log_pf + Ne20_log_pf - He4_log_pf - Si28_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dC12_log_pf_dT9 + dNe20_log_pf_dT9 - dHe4_log_pf_dT9 - dSi28_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from rolfr
    ln_set_rate =  -307.7623006031339 + -186.77284344745925 * tfactors.T9i + 514.197 * tfactors.T913i + -200.896 * tfactors.T913
                         + -6.42713 * tfactors.T9 + 0.758256 * tfactors.T953 + 236.359 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  186.77284344745925 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 514.197 * tfactors.T943i + (1.0/3.0) * -200.896 * tfactors.T923i
                                  + -6.42713 + (5.0/3.0) * 0.758256 * tfactors.T923 + 236.359 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Si28_to_O16_O16_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Si28 + He4 --> O16 + O16

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real O16_log_pf, dO16_log_pf_dT9;

    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors.T9, pf_cache, O16_log_pf, dO16_log_pf_dT9);

    amrex::Real Si28_log_pf, dSi28_log_pf_dT9;

    // interpolating Si28 partition function
    get_partition_function_cached(Si28, tfactors.T9, pf_cache, Si28_log_pf, dSi28_log_pf_dT9);

    amrex::Real net_log_pf = O16_log_pf + O16_log_pf - He4_log_pf - Si28_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dO16_log_pf_dT9 + dO16_log_pf_dT9 - dHe4_log_pf_dT9 - dSi28_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from cf88r
    ln_set_rate =  97.78921572069022 + -111.60083385054465 * tfactors.T9i + -119.324 * tfactors.T913i + -32.2497 * tfactors.T913
                         + 1.46214 * tfactors.T9 + -0.200893 * tfactors.T953 + 13.2148 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  111.60083385054465 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -119.324 * tfactors.T943i + (1.0/3.0) * -32.2497 * tfactors.T923i
                                  + 1.46214 + (5.0/3.0) * -0.200893 * tfactors.T923 + 13.2148 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Si29_to_He4_Al26_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Si29 + p --> He4 + Al26

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Al26_log_pf, dAl26_log_pf_dT9;

    // interpolating Al26 partition function
    get_partition_function_cached(Al26, tfactors.T9, pf_cache, Al26_log_pf, dAl26_log_pf_dT9);

    amrex::Real Si29_log_pf, dSi29_log_pf_dT9;

    // interpolating Si29 partition function
    get_partition_function_cached(Si29, tfactors.T9, pf_cache, Si29_log_pf, dSi29_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Al26_log_pf - p_log_pf - Si29_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dAl26_log_pf_dT9 - dp_log_pf_dT9 - dSi29_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  50.61935001156578 + -55.94406249988892 * tfactors.T9i + -56.4422 * tfactors.T913i + 0.705353 * tfactors.T913
                         + -0.957427 * tfactors.T9 + 0.0756045 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  55.94406249988892 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -56.4422 * tfactors.T943i + (1.0/3.0) * 0.705353 * tfactors.T923i
                                  + -0.957427 + (5.0/3.0) * 0.0756045 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_P31_to_C12_Ne20_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // P31 + p --> C12 + Ne20

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real C12_log_pf, dC12_log_pf_dT9;

    // setting C12 log(partition function) to 0.0 by default, independent of T
    C12_log_pf = 0.0_rt;
    dC12_log_pf_dT9 = 0.0_rt;

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real P31_log_pf, dP31_log_pf_dT9;

    // interpolating P31 partition function
    get_partition_function_cached(P31, tfactors.T9, pf_cache, P31_log_pf, dP31_log_pf_dT9);

    amrex::Real Ne20_log_pf, dNe20_log_pf_dT9;

    // interpolating Ne20 partition function
    get_partition_function_cached(Ne20, tfactors.T9, pf_cache, Ne20_log_pf, dNe20_log_pf_dT9);

    amrex::Real net_log_pf = C12_log_pf + Ne20_log_pf - p_log_pf - P31_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dC12_log_pf_dT9 + dNe20_log_pf_dT9 - dp_log_pf_dT9 - dP31_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from rolfr
    ln_set_rate =  -266.47063151346725 + -156.07991376506303 * tfactors.T9i + 361.154 * tfactors.T913i + -92.643 * tfactors.T913
                         + -9.98738 * tfactors.T9 + 0.892737 * tfactors.T953 + 161.042 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  156.07991376506303 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 361.154 * tfactors.T943i + (1.0/3.0) * -92.643 * tfactors.T923i
                                  + -9.98738 + (5.0/3.0) * 0.892737 * tfactors.T923 + 161.042 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_P31_to_O16_O16_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // P31 + p --> O16 + O16

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real O16_log_pf, dO16_log_pf_dT9;

    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors.T9, pf_cache, O16_log_pf, dO16_log_pf_dT9);

    amrex::Real P31_log_pf, dP31_log_pf_dT9;

    // interpolating P31 partition function
    get_partition_function_cached(P31, tfactors.T9, pf_cache, P31_log_pf, dP31_log_pf_dT9);

    amrex::Real net_log_pf = O16_log_pf + O16_log_pf - p_log_pf - P31_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dO16_log_pf_dT9 + dO16_log_pf_dT9 - dp_log_pf_dT9 - dP31_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from cf88r
    ln_set_rate =  86.33118481035693 + -88.87103716814843 * tfactors.T9i + -145.844 * tfactors.T913i + 8.72612 * tfactors.T913
                         + -0.554035 * tfactors.T9 + -0.137562 * tfactors.T953 + -6.88807 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  88.87103716814843 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -145.844 * tfactors.T943i + (1.0/3.0) * 8.72612 * tfactors.T923i
                                  + -0.554035 + (5.0/3.0) * -0.137562 * tfactors.T923 + -6.88807 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_S32_to_p_Cl35_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // S32 + He4 --> p + Cl35

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Cl35_log_pf, dCl35_log_pf_dT9;

    // interpolating Cl35 partition function
    get_partition_function_cached(Cl35, tfactors.T9, pf_cache, Cl35_log_pf, dCl35_log_pf_dT9);

    amrex::Real S32_log_pf, dS32_log_pf_dT9;

    // interpolating S32 partition function
    get_partition_function_cached(S32, tfactors.T9, pf_cache, S32_log_pf, dS32_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = p_log_pf + Cl35_log_pf - He4_log_pf - S32_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dCl35_log_pf_dT9 - dHe4_log_pf_dT9 - dS32_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from il10r
    ln_set_rate =  -0.8597847344696254 + -25.589703706656966 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  25.589703706656966 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  -57.37716473446963 + -22.187684706656967 * tfactors.T9i + 25.5338 * tfactors.T913
                         + 6.45824 * tfactors.T9 + -0.950294 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  22.187684706656967 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 25.5338 * tfactors.T923i
                                  + 6.45824 + (5.0/3.0) * -0.950294 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10n
    ln_set_rate =  32.27223526553037 + -21.654753706656965 * tfactors.T9i + -30.9147 * tfactors.T913i + -1.2345 * tfactors.T913
                         + 22.5118 * tfactors.T9 + -33.0589 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  21.654753706656965 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -30.9147 * tfactors.T943i + (1.0/3.0) * -1.2345 * tfactors.T923i
                                  + 22.5118 + (5.0/3.0) * -33.0589 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  2.4434452655303747 + -27.664513706656965 * tfactors.T9i + 5.33756 * tfactors.T913
                         + 1.64418 * tfactors.T9 + -0.246167 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  27.664513706656965 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 5.33756 * tfactors.T923i
                                  + 1.64418 + (5.0/3.0) * -0.246167 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_S33_to_He4_P30_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // S33 + p --> He4 + P30

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real S33_log_pf, dS33_log_pf_dT9;

    // interpolating S33 partition function
    get_partition_function_cached(S33, tfactors.T9, pf_cache, S33_log_pf, dS33_log_pf_dT9);

    amrex::Real P30_log_pf, dP30_log_pf_dT9;

    // interpolating P30 partition function
    get_partition_function_cached(P30, tfactors.T9, pf_cache, P30_log_pf, dP30_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + P30_log_pf - p_log_pf - S33_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dP30_log_pf_dT9 - dp_log_pf_dT9 - dS33_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  50.751142362216704 + -17.649811835453008 * tfactors.T9i + -62.4643 * tfactors.T913i + 0.492934 * tfactors.T913
                         + -0.841855 * tfactors.T9 + 0.059263 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  17.649811835453008 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -62.4643 * tfactors.T943i + (1.0/3.0) * 0.492934 * tfactors.T923i
                                  + -0.841855 + (5.0/3.0) * 0.059263 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Ar36_to_p_K39_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ar36 + He4 --> p + K39

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real K39_log_pf, dK39_log_pf_dT9;

    // interpolating K39 partition function
    get_partition_function_cached(K39, tfactors.T9, pf_cache, K39_log_pf, dK39_log_pf_dT9);

    amrex::Real Ar36_log_pf, dAr36_log_pf_dT9;

    // interpolating Ar36 partition function
    get_partition_function_cached(Ar36, tfactors.T9, pf_cache, Ar36_log_pf, dAr36_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = p_log_pf + K39_log_pf - He4_log_pf - Ar36_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dK39_log_pf_dT9 - dHe4_log_pf_dT9 - dAr36_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  20.654451232729432 + -14.951286608188834 * tfactors.T9i + -30.0732 * tfactors.T913i + 7.03263 * tfactors.T913
                         + -1.10085 * tfactors.T9 + 0.133768 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  14.951286608188834 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -30.0732 * tfactors.T943i + (1.0/3.0) * 7.03263 * tfactors.T923i
                                  + -1.10085 + (5.0/3.0) * 0.133768 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Ar37_to_He4_Cl34_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ar37 + p --> He4 + Cl34

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Cl34_log_pf, dCl34_log_pf_dT9;

    // interpolating Cl34 partition function
    get_partition_function_cached(Cl34, tfactors.T9, pf_cache, Cl34_log_pf, dCl34_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ar37_log_pf, dAr37_log_pf_dT9;

    // interpolating Ar37 partition function
    get_partition_function_cached(Ar37, tfactors.T9, pf_cache, Ar37_log_pf, dAr37_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Cl34_log_pf - p_log_pf - Ar37_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCl34_log_pf_dT9 - dp_log_pf_dT9 - dAr37_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  52.65061986402699 + -19.07242943961376 * tfactors.T9i + -68.2165 * tfactors.T913i + 0.330057 * tfactors.T913
                         + -0.873334 * tfactors.T9 + 0.0592127 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  19.07242943961376 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -68.2165 * tfactors.T943i + (1.0/3.0) * 0.330057 * tfactors.T923i
                                  + -0.873334 + (5.0/3.0) * 0.0592127 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Ca40_to_p_Sc43_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ca40 + He4 --> p + Sc43

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Ca40_log_pf, dCa40_log_pf_dT9;

    // interpolating Ca40 partition function
    get_partition_function_cached(Ca40, tfactors.T9, pf_cache, Ca40_log_pf, dCa40_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Sc43_log_pf, dSc43_log_pf_dT9;

    // interpolating Sc43 partition function
    get_partition_function_cached(Sc43, tfactors.T9, pf_cache, Sc43_log_pf, dSc43_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Sc43_log_pf - He4_log_pf - Ca40_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dSc43_log_pf_dT9 - dHe4_log_pf_dT9 - dCa40_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  35.67546755788414 + -40.87525788838128 * tfactors.T9i + -32.1734 * tfactors.T913i + 0.0296879 * tfactors.T913
                         + -0.95232 * tfactors.T9 + 0.129022 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  40.87525788838128 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -32.1734 * tfactors.T943i + (1.0/3.0) * 0.0296879 * tfactors.T923i
                                  + -0.95232 + (5.0/3.0) * 0.129022 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Ca41_to_He4_K38_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ca41 + p --> He4 + K38

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ca41_log_pf, dCa41_log_pf_dT9;

    // interpolating Ca41 partition function
    get_partition_function_cached(Ca41, tfactors.T9, pf_cache, Ca41_log_pf, dCa41_log_pf_dT9);

    amrex::Real K38_log_pf, dK38_log_pf_dT9;

    // interpolating K38 partition function
    get_partition_function_cached(K38, tfactors.T9, pf_cache, K38_log_pf, dK38_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + K38_log_pf - p_log_pf - Ca41_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dK38_log_pf_dT9 - dp_log_pf_dT9 - dCa41_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  55.541590396195424 + -17.09455537095993 * tfactors.T9i + -73.7407 * tfactors.T913i + 0.600725 * tfactors.T913
                         + -0.957475 * tfactors.T9 + 0.0689286 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  17.09455537095993 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -73.7407 * tfactors.T943i + (1.0/3.0) * 0.600725 * tfactors.T923i
                                  + -0.957475 + (5.0/3.0) * 0.0689286 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Ti45_to_He4_Sc42_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ti45 + p --> He4 + Sc42

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Sc42_log_pf, dSc42_log_pf_dT9;

    // interpolating Sc42 partition function
    get_partition_function_cached(Sc42, tfactors.T9, pf_cache, Sc42_log_pf, dSc42_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ti45_log_pf, dTi45_log_pf_dT9;

    // interpolating Ti45 partition function
    get_partition_function_cached(Ti45, tfactors.T9, pf_cache, Ti45_log_pf, dTi45_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Sc42_log_pf - p_log_pf - Ti45_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dSc42_log_pf_dT9 - dp_log_pf_dT9 - dTi45_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  56.88137786781152 + -23.50199005180725 * tfactors.T9i + -79.0687 * tfactors.T913i + -0.981235 * tfactors.T913
                         + -0.48951 * tfactors.T9 + 0.0114272 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  23.50199005180725 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -79.0687 * tfactors.T943i + (1.0/3.0) * -0.981235 * tfactors.T923i
                                  + -0.48951 + (5.0/3.0) * 0.0114272 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_V47_to_He4_Ti44_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // V47 + p --> He4 + Ti44

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real V47_log_pf, dV47_log_pf_dT9;

    // interpolating V47 partition function
    get_partition_function_cached(V47, tfactors.T9, pf_cache, V47_log_pf, dV47_log_pf_dT9);

    amrex::Real Ti44_log_pf, dTi44_log_pf_dT9;

    // interpolating Ti44 partition function
    get_partition_function_cached(Ti44, tfactors.T9, pf_cache, Ti44_log_pf, dTi44_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Ti44_log_pf - p_log_pf - V47_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dTi44_log_pf_dT9 - dp_log_pf_dT9 - dV47_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from chw0r
    ln_set_rate =  -76.63208838152472 + -6.086479266444236 * tfactors.T9i + 70.2835 * tfactors.T913
                         + -7.99061 * tfactors.T9 + 0.486213 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  6.086479266444236 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 70.2835 * tfactors.T923i
                                  + -7.99061 + (5.0/3.0) * 0.486213 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Cr49_to_He4_V46_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Cr49 + p --> He4 + V46

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real V46_log_pf, dV46_log_pf_dT9;

    // interpolating V46 partition function
    get_partition_function_cached(V46, tfactors.T9, pf_cache, V46_log_pf, dV46_log_pf_dT9);

    amrex::Real Cr49_log_pf, dCr49_log_pf_dT9;

    // interpolating Cr49 partition function
    get_partition_function_cached(Cr49, tfactors.T9, pf_cache, Cr49_log_pf, dCr49_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + V46_log_pf - p_log_pf - Cr49_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dV46_log_pf_dT9 - dp_log_pf_dT9 - dCr49_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  59.280226694242224 + -39.36768722760379 * tfactors.T9i + -84.2255 * tfactors.T913i + -1.09776 * tfactors.T913
                         + -0.495216 * tfactors.T9 + 0.0124719 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  39.36768722760379 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -84.2255 * tfactors.T943i + (1.0/3.0) * -1.09776 * tfactors.T923i
                                  + -0.495216 + (5.0/3.0) * 0.0124719 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Mn51_to_He4_Cr48_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mn51 + p --> He4 + Cr48

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Cr48_log_pf, dCr48_log_pf_dT9;

    // interpolating Cr48 partition function
    get_partition_function_cached(Cr48, tfactors.T9, pf_cache, Cr48_log_pf, dCr48_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Mn51_log_pf, dMn51_log_pf_dT9;

    // interpolating Mn51 partition function
    get_partition_function_cached(Mn51, tfactors.T9, pf_cache, Mn51_log_pf, dMn51_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Cr48_log_pf - p_log_pf - Mn51_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCr48_log_pf_dT9 - dp_log_pf_dT9 - dMn51_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  58.71348267023461 + -6.47700151570639 * tfactors.T9i + -86.7459 * tfactors.T913i + 1.05653 * tfactors.T913
                         + -1.15757 * tfactors.T9 + 0.0877546 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  6.47700151570639 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -86.7459 * tfactors.T943i + (1.0/3.0) * 1.05653 * tfactors.T923i
                                  + -1.15757 + (5.0/3.0) * 0.0877546 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Mn52_to_He4_Cr49_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mn52 + p --> He4 + Cr49

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Mn52_log_pf, dMn52_log_pf_dT9;

    // interpolating Mn52 partition function
    get_partition_function_cached(Mn52, tfactors.T9, pf_cache, Mn52_log_pf, dMn52_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Cr49_log_pf, dCr49_log_pf_dT9;

    // interpolating Cr49 partition function
    get_partition_function_cached(Cr49, tfactors.T9, pf_cache, Cr49_log_pf, dCr49_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Cr49_log_pf - p_log_pf - Mn52_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCr49_log_pf_dT9 - dp_log_pf_dT9 - dMn52_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  53.01374392177941 + -5.975570287637617 * tfactors.T9i + -86.7913 * tfactors.T913i + 10.6892 * tfactors.T913
                         + -4.01586 * tfactors.T9 + 0.411275 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  5.975570287637617 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -86.7913 * tfactors.T943i + (1.0/3.0) * 10.6892 * tfactors.T923i
                                  + -4.01586 + (5.0/3.0) * 0.411275 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Co55_to_He4_Fe52_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co55 + p --> He4 + Fe52

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Co55_log_pf, dCo55_log_pf_dT9;

    // interpolating Co55 partition function
    get_partition_function_cached(Co55, tfactors.T9, pf_cache, Co55_log_pf, dCo55_log_pf_dT9);

    amrex::Real Fe52_log_pf, dFe52_log_pf_dT9;

    // interpolating Fe52 partition function
    get_partition_function_cached(Fe52, tfactors.T9, pf_cache, Fe52_log_pf, dFe52_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Fe52_log_pf - p_log_pf - Co55_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dFe52_log_pf_dT9 - dp_log_pf_dT9 - dCo55_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  61.42570080881217 + -9.676367161770125 * tfactors.T9i + -91.6819 * tfactors.T913i + -0.329235 * tfactors.T913
                         + -0.780924 * tfactors.T9 + 0.0425179 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  9.676367161770125 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -91.6819 * tfactors.T943i + (1.0/3.0) * -0.329235 * tfactors.T923i
                                  + -0.780924 + (5.0/3.0) * 0.0425179 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Mn51_to_n_Mn50_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mn51 --> n + Mn50

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Mn50_log_pf, dMn50_log_pf_dT9;

    // interpolating Mn50 partition function
    get_partition_function_cached(Mn50, tfactors.T9, pf_cache, Mn50_log_pf, dMn50_log_pf_dT9);

    amrex::Real Mn51_log_pf, dMn51_log_pf_dT9;

    // interpolating Mn51 partition function
    get_partition_function_cached(Mn51, tfactors.T9, pf_cache, Mn51_log_pf, dMn51_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Mn50_log_pf - Mn51_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dMn50_log_pf_dT9 - dMn51_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  37.34951332478154 + -158.83832832749732 * tfactors.T9i + -0.538879 * tfactors.T913
                         + 0.284528 * tfactors.T9 + -0.0459849 * tfactors.T953 + 1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  158.83832832749732 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.538879 * tfactors.T923i
                                  + 0.284528 + (5.0/3.0) * -0.0459849 * tfactors.T923 + 1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Mn52_to_n_Mn51_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mn52 --> n + Mn51

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Mn51_log_pf, dMn51_log_pf_dT9;

    // interpolating Mn51 partition function
    get_partition_function_cached(Mn51, tfactors.T9, pf_cache, Mn51_log_pf, dMn51_log_pf_dT9);

    amrex::Real Mn52_log_pf, dMn52_log_pf_dT9;

    // interpolating Mn52 partition function
    get_partition_function_cached(Mn52, tfactors.T9, pf_cache, Mn52_log_pf, dMn52_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Mn51_log_pf - Mn52_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dMn51_log_pf_dT9 - dMn52_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  40.130587431823315 + -122.30591273875024 * tfactors.T9i + -1.83611 * tfactors.T913
                         + 0.460384 * tfactors.T9 + -0.0584947 * tfactors.T953 + 1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  122.30591273875024 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.83611 * tfactors.T923i
                                  + 0.460384 + (5.0/3.0) * -0.0584947 * tfactors.T923 + 1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Fe53_to_n_Fe52_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe53 --> n + Fe52

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Fe52_log_pf, dFe52_log_pf_dT9;

    // interpolating Fe52 partition function
    get_partition_function_cached(Fe52, tfactors.T9, pf_cache, Fe52_log_pf, dFe52_log_pf_dT9);

    amrex::Real Fe53_log_pf, dFe53_log_pf_dT9;

    // interpolating Fe53 partition function
    get_partition_function_cached(Fe53, tfactors.T9, pf_cache, Fe53_log_pf, dFe53_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Fe52_log_pf - Fe53_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dFe52_log_pf_dT9 - dFe53_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  37.50789742709652 + -124.01421385132892 * tfactors.T9i + -0.344319 * tfactors.T913
                         + 0.178277 * tfactors.T9 + -0.0334326 * tfactors.T953 + 1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  124.01421385132892 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.344319 * tfactors.T923i
                                  + 0.178277 + (5.0/3.0) * -0.0334326 * tfactors.T923 + 1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Fe53_to_p_Mn52_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe53 --> p + Mn52

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Mn52_log_pf, dMn52_log_pf_dT9;

    // interpolating Mn52 partition function
    get_partition_function_cached(Mn52, tfactors.T9, pf_cache, Mn52_log_pf, dMn52_log_pf_dT9);

    amrex::Real Fe53_log_pf, dFe53_log_pf_dT9;

    // interpolating Fe53 partition function
    get_partition_function_cached(Fe53, tfactors.T9, pf_cache, Fe53_log_pf, dFe53_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Mn52_log_pf - Fe53_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dMn52_log_pf_dT9 - dFe53_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  59.6140413175536 + -87.32493957326004 * tfactors.T9i + -36.187 * tfactors.T913i + 1.79171 * tfactors.T913
                         + -1.77786 * tfactors.T9 + 0.179936 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  87.32493957326004 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -36.187 * tfactors.T943i + (1.0/3.0) * 1.79171 * tfactors.T923i
                                  + -1.77786 + (5.0/3.0) * 0.179936 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Fe54_to_n_Fe53_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe54 --> n + Fe53

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Fe53_log_pf, dFe53_log_pf_dT9;

    // interpolating Fe53 partition function
    get_partition_function_cached(Fe53, tfactors.T9, pf_cache, Fe53_log_pf, dFe53_log_pf_dT9);

    amrex::Real Fe54_log_pf, dFe54_log_pf_dT9;

    // interpolating Fe54 partition function
    get_partition_function_cached(Fe54, tfactors.T9, pf_cache, Fe54_log_pf, dFe54_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Fe53_log_pf - Fe54_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dFe53_log_pf_dT9 - dFe54_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  42.23234536763494 + -155.25009527915816 * tfactors.T9i + -1.10421 * tfactors.T913
                         + 0.379905 * tfactors.T9 + -0.0581878 * tfactors.T953 + 1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  155.25009527915816 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.10421 * tfactors.T923i
                                  + 0.379905 + (5.0/3.0) * -0.0581878 * tfactors.T923 + 1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Fe55_to_n_Fe54_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe55 --> n + Fe54

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Fe55_log_pf, dFe55_log_pf_dT9;

    // interpolating Fe55 partition function
    get_partition_function_cached(Fe55, tfactors.T9, pf_cache, Fe55_log_pf, dFe55_log_pf_dT9);

    amrex::Real Fe54_log_pf, dFe54_log_pf_dT9;

    // interpolating Fe54 partition function
    get_partition_function_cached(Fe54, tfactors.T9, pf_cache, Fe54_log_pf, dFe54_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Fe54_log_pf - Fe55_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dFe54_log_pf_dT9 - dFe55_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ks03
    ln_set_rate =  21.50487677145662 + -107.84100838776591 * tfactors.T9i + -8.66617 * tfactors.T913i + 26.4472 * tfactors.T913
                         + -1.9222 * tfactors.T9 + 0.0986404 * tfactors.T953 + -8.28317 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  107.84100838776591 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -8.66617 * tfactors.T943i + (1.0/3.0) * 26.4472 * tfactors.T923i
                                  + -1.9222 + (5.0/3.0) * 0.0986404 * tfactors.T923 + -8.28317 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Fe56_to_n_Fe55_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe56 --> n + Fe55

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Fe56_log_pf, dFe56_log_pf_dT9;

    // interpolating Fe56 partition function
    get_partition_function_cached(Fe56, tfactors.T9, pf_cache, Fe56_log_pf, dFe56_log_pf_dT9);

    amrex::Real Fe55_log_pf, dFe55_log_pf_dT9;

    // interpolating Fe55 partition function
    get_partition_function_cached(Fe55, tfactors.T9, pf_cache, Fe55_log_pf, dFe55_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Fe55_log_pf - Fe56_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dFe55_log_pf_dT9 - dFe56_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ks03
    ln_set_rate =  46.806892234523204 + -130.0322634199283 * tfactors.T9i + 8.06062 * tfactors.T913i + -14.4809 * tfactors.T913
                         + 0.94252 * tfactors.T9 + -0.0776007 * tfactors.T953 + 7.97093 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  130.0322634199283 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 8.06062 * tfactors.T943i + (1.0/3.0) * -14.4809 * tfactors.T923i
                                  + 0.94252 + (5.0/3.0) * -0.0776007 * tfactors.T923 + 7.97093 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Co55_to_p_Fe54_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co55 --> p + Fe54

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Co55_log_pf, dCo55_log_pf_dT9;

    // interpolating Co55 partition function
    get_partition_function_cached(Co55, tfactors.T9, pf_cache, Co55_log_pf, dCo55_log_pf_dT9);

    amrex::Real Fe54_log_pf, dFe54_log_pf_dT9;

    // interpolating Fe54 partition function
    get_partition_function_cached(Fe54, tfactors.T9, pf_cache, Fe54_log_pf, dFe54_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Fe54_log_pf - Co55_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dFe54_log_pf_dT9 - dCo55_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  57.84851844810452 + -58.76958578644693 * tfactors.T9i + -37.1544 * tfactors.T913i + 0.950364 * tfactors.T913
                         + -1.77529 * tfactors.T9 + 0.198562 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  58.76958578644693 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -37.1544 * tfactors.T943i + (1.0/3.0) * 0.950364 * tfactors.T923i
                                  + -1.77529 + (5.0/3.0) * 0.198562 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Co56_to_n_Co55_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co56 --> n + Co55

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Co55_log_pf, dCo55_log_pf_dT9;

    // interpolating Co55 partition function
    get_partition_function_cached(Co55, tfactors.T9, pf_cache, Co55_log_pf, dCo55_log_pf_dT9);

    amrex::Real Co56_log_pf, dCo56_log_pf_dT9;

    // interpolating Co56 partition function
    get_partition_function_cached(Co56, tfactors.T9, pf_cache, Co56_log_pf, dCo56_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Co55_log_pf - Co56_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dCo55_log_pf_dT9 - dCo56_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  40.647564333711024 + -116.994640839646 * tfactors.T9i + -1.86357 * tfactors.T913
                         + 0.616591 * tfactors.T9 + -0.0839313 * tfactors.T953 + 1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  116.994640839646 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.86357 * tfactors.T923i
                                  + 0.616591 + (5.0/3.0) * -0.0839313 * tfactors.T923 + 1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Co56_to_p_Fe55_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co56 --> p + Fe55

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Fe55_log_pf, dFe55_log_pf_dT9;

    // interpolating Fe55 partition function
    get_partition_function_cached(Fe55, tfactors.T9, pf_cache, Fe55_log_pf, dFe55_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Co56_log_pf, dCo56_log_pf_dT9;

    // interpolating Co56 partition function
    get_partition_function_cached(Co56, tfactors.T9, pf_cache, Co56_log_pf, dCo56_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Fe55_log_pf - Co56_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dFe55_log_pf_dT9 - dCo56_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  58.620466010358925 + -67.86404663834863 * tfactors.T9i + -37.1585 * tfactors.T913i + 1.66198 * tfactors.T913
                         + -1.60842 * tfactors.T9 + 0.148916 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  67.86404663834863 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -37.1585 * tfactors.T943i + (1.0/3.0) * 1.66198 * tfactors.T923i
                                  + -1.60842 + (5.0/3.0) * 0.148916 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Co56_to_He4_Mn52_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co56 --> He4 + Mn52

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Mn52_log_pf, dMn52_log_pf_dT9;

    // interpolating Mn52 partition function
    get_partition_function_cached(Mn52, tfactors.T9, pf_cache, Mn52_log_pf, dMn52_log_pf_dT9);

    amrex::Real Co56_log_pf, dCo56_log_pf_dT9;

    // interpolating Co56 partition function
    get_partition_function_cached(Co56, tfactors.T9, pf_cache, Co56_log_pf, dCo56_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Mn52_log_pf - Co56_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dMn52_log_pf_dT9 - dCo56_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  88.91600903298027 + -89.98173372334725 * tfactors.T9i + -89.3157 * tfactors.T913i + -1.95982 * tfactors.T913
                         + -2.45671 * tfactors.T9 + 0.296095 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  89.98173372334725 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -89.3157 * tfactors.T943i + (1.0/3.0) * -1.95982 * tfactors.T923i
                                  + -2.45671 + (5.0/3.0) * 0.296095 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Co57_to_n_Co56_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co57 --> n + Co56

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Co56_log_pf, dCo56_log_pf_dT9;

    // interpolating Co56 partition function
    get_partition_function_cached(Co56, tfactors.T9, pf_cache, Co56_log_pf, dCo56_log_pf_dT9);

    amrex::Real Co57_log_pf, dCo57_log_pf_dT9;

    // interpolating Co57 partition function
    get_partition_function_cached(Co57, tfactors.T9, pf_cache, Co57_log_pf, dCo57_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Co56_log_pf - Co57_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dCo56_log_pf_dT9 - dCo57_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  41.17386809654295 + -132.01901045161702 * tfactors.T9i + -1.37855 * tfactors.T913
                         + 0.299896 * tfactors.T9 + -0.04382 * tfactors.T953 + 1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  132.01901045161702 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.37855 * tfactors.T923i
                                  + 0.299896 + (5.0/3.0) * -0.04382 * tfactors.T923 + 1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Co57_to_p_Fe56_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co57 --> p + Fe56

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Fe56_log_pf, dFe56_log_pf_dT9;

    // interpolating Fe56 partition function
    get_partition_function_cached(Fe56, tfactors.T9, pf_cache, Fe56_log_pf, dFe56_log_pf_dT9);

    amrex::Real Co57_log_pf, dCo57_log_pf_dT9;

    // interpolating Co57 partition function
    get_partition_function_cached(Co57, tfactors.T9, pf_cache, Co57_log_pf, dCo57_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Fe56_log_pf - Co57_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dFe56_log_pf_dT9 - dCo57_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  57.68564187237867 + -69.94636137001572 * tfactors.T9i + -37.1625 * tfactors.T913i + 1.06776 * tfactors.T913
                         + -1.31689 * tfactors.T9 + 0.122089 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  69.94636137001572 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -37.1625 * tfactors.T943i + (1.0/3.0) * 1.06776 * tfactors.T923i
                                  + -1.31689 + (5.0/3.0) * 0.122089 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ni57_to_n_Ni56_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni57 --> n + Ni56

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Ni56_log_pf, dNi56_log_pf_dT9;

    // interpolating Ni56 partition function
    get_partition_function_cached(Ni56, tfactors.T9, pf_cache, Ni56_log_pf, dNi56_log_pf_dT9);

    amrex::Real Ni57_log_pf, dNi57_log_pf_dT9;

    // interpolating Ni57 partition function
    get_partition_function_cached(Ni57, tfactors.T9, pf_cache, Ni57_log_pf, dNi57_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Ni56_log_pf - Ni57_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dNi56_log_pf_dT9 - dNi57_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  38.391039283996996 + -118.91983039605458 * tfactors.T9i + -1.19665 * tfactors.T913
                         + 0.507179 * tfactors.T9 + -0.074604 * tfactors.T953 + 1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  118.91983039605458 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.19665 * tfactors.T923i
                                  + 0.507179 + (5.0/3.0) * -0.074604 * tfactors.T923 + 1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ni57_to_p_Co56_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni57 --> p + Co56

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Ni57_log_pf, dNi57_log_pf_dT9;

    // interpolating Ni57 partition function
    get_partition_function_cached(Ni57, tfactors.T9, pf_cache, Ni57_log_pf, dNi57_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Co56_log_pf, dCo56_log_pf_dT9;

    // interpolating Co56 partition function
    get_partition_function_cached(Co56, tfactors.T9, pf_cache, Co56_log_pf, dCo56_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Co56_log_pf - Ni57_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dCo56_log_pf_dT9 - dNi57_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  60.10327899221158 + -85.08979333790236 * tfactors.T9i + -38.1094 * tfactors.T913i + 2.57091 * tfactors.T913
                         + -2.07795 * tfactors.T9 + 0.20757 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  85.08979333790236 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -38.1094 * tfactors.T943i + (1.0/3.0) * 2.57091 * tfactors.T923i
                                  + -2.07795 + (5.0/3.0) * 0.20757 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ni57_to_He4_Fe53_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni57 --> He4 + Fe53

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Ni57_log_pf, dNi57_log_pf_dT9;

    // interpolating Ni57 partition function
    get_partition_function_cached(Ni57, tfactors.T9, pf_cache, Ni57_log_pf, dNi57_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Fe53_log_pf, dFe53_log_pf_dT9;

    // interpolating Fe53 partition function
    get_partition_function_cached(Fe53, tfactors.T9, pf_cache, Fe53_log_pf, dFe53_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Fe53_log_pf - Ni57_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dFe53_log_pf_dT9 - dNi57_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  94.74134670763823 + -87.74658748798956 * tfactors.T9i + -91.7231 * tfactors.T913i + -10.1816 * tfactors.T913
                         + -0.0406303 * tfactors.T9 + 0.0345056 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  87.74658748798956 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -91.7231 * tfactors.T943i + (1.0/3.0) * -10.1816 * tfactors.T923i
                                  + -0.0406303 + (5.0/3.0) * 0.0345056 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ni58_to_n_Ni57_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni58 --> n + Ni57

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Ni57_log_pf, dNi57_log_pf_dT9;

    // interpolating Ni57 partition function
    get_partition_function_cached(Ni57, tfactors.T9, pf_cache, Ni57_log_pf, dNi57_log_pf_dT9);

    amrex::Real Ni58_log_pf, dNi58_log_pf_dT9;

    // interpolating Ni58 partition function
    get_partition_function_cached(Ni58, tfactors.T9, pf_cache, Ni58_log_pf, dNi58_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Ni57_log_pf - Ni58_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dNi57_log_pf_dT9 - dNi58_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  42.36077897558323 + -141.76332431826057 * tfactors.T9i + -1.90814 * tfactors.T913
                         + 0.493188 * tfactors.T9 + -0.0684633 * tfactors.T953 + 1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  141.76332431826057 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.90814 * tfactors.T923i
                                  + 0.493188 + (5.0/3.0) * -0.0684633 * tfactors.T923 + 1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ni58_to_p_Co57_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni58 --> p + Co57

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Ni58_log_pf, dNi58_log_pf_dT9;

    // interpolating Ni58 partition function
    get_partition_function_cached(Ni58, tfactors.T9, pf_cache, Ni58_log_pf, dNi58_log_pf_dT9);

    amrex::Real Co57_log_pf, dCo57_log_pf_dT9;

    // interpolating Co57 partition function
    get_partition_function_cached(Co57, tfactors.T9, pf_cache, Co57_log_pf, dCo57_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Co57_log_pf - Ni58_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dCo57_log_pf_dT9 - dNi58_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  61.794489871251855 + -94.83410720454592 * tfactors.T9i + -38.1133 * tfactors.T913i + 1.77414 * tfactors.T913
                         + -1.48268 * tfactors.T9 + 0.121073 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  94.83410720454592 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -38.1133 * tfactors.T943i + (1.0/3.0) * 1.77414 * tfactors.T923i
                                  + -1.48268 + (5.0/3.0) * 0.121073 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Ni58_to_He4_Fe54_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni58 --> He4 + Fe54

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ni58_log_pf, dNi58_log_pf_dT9;

    // interpolating Ni58 partition function
    get_partition_function_cached(Ni58, tfactors.T9, pf_cache, Ni58_log_pf, dNi58_log_pf_dT9);

    amrex::Real Fe54_log_pf, dFe54_log_pf_dT9;

    // interpolating Fe54 partition function
    get_partition_function_cached(Fe54, tfactors.T9, pf_cache, Fe54_log_pf, dFe54_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Fe54_log_pf - Ni58_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dFe54_log_pf_dT9 - dNi58_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  85.23428031558652 + -74.25981652709197 * tfactors.T9i + -91.7628 * tfactors.T913i + 4.23027 * tfactors.T913
                         + -3.31305 * tfactors.T9 + 0.271293 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  74.25981652709197 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -91.7628 * tfactors.T943i + (1.0/3.0) * 4.23027 * tfactors.T923i
                                  + -3.31305 + (5.0/3.0) * 0.271293 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Cu59_to_p_Ni58_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Cu59 --> p + Ni58

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Cu59_log_pf, dCu59_log_pf_dT9;

    // interpolating Cu59 partition function
    get_partition_function_cached(Cu59, tfactors.T9, pf_cache, Cu59_log_pf, dCu59_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Ni58_log_pf, dNi58_log_pf_dT9;

    // interpolating Ni58 partition function
    get_partition_function_cached(Ni58, tfactors.T9, pf_cache, Ni58_log_pf, dNi58_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Ni58_log_pf - Cu59_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dNi58_log_pf_dT9 - dCu59_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from si13n
    ln_set_rate =  58.97928371193912 + -39.67086986195644 * tfactors.T9i + -39.0526 * tfactors.T913i + 1.0436 * tfactors.T913
                         + -2.10834 * tfactors.T9 + 0.239559 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  39.67086986195644 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -39.0526 * tfactors.T943i + (1.0/3.0) * 1.0436 * tfactors.T923i
                                  + -2.10834 + (5.0/3.0) * 0.239559 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Cu59_to_He4_Co55_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Cu59 --> He4 + Co55

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Cu59_log_pf, dCu59_log_pf_dT9;

    // interpolating Cu59 partition function
    get_partition_function_cached(Cu59, tfactors.T9, pf_cache, Cu59_log_pf, dCu59_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Co55_log_pf, dCo55_log_pf_dT9;

    // interpolating Co55 partition function
    get_partition_function_cached(Co55, tfactors.T9, pf_cache, Co55_log_pf, dCo55_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Co55_log_pf - Cu59_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCo55_log_pf_dT9 - dCu59_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  90.97354557942113 + -55.1611006026015 * tfactors.T9i + -94.1404 * tfactors.T913i + -2.62786 * tfactors.T913
                         + -2.12066 * tfactors.T9 + 0.237999 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  55.1611006026015 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -94.1404 * tfactors.T943i + (1.0/3.0) * -2.62786 * tfactors.T923i
                                  + -2.12066 + (5.0/3.0) * 0.237999 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Zn60_to_p_Cu59_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Zn60 --> p + Cu59

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Cu59_log_pf, dCu59_log_pf_dT9;

    // interpolating Cu59 partition function
    get_partition_function_cached(Cu59, tfactors.T9, pf_cache, Cu59_log_pf, dCu59_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Zn60_log_pf, dZn60_log_pf_dT9;

    // interpolating Zn60 partition function
    get_partition_function_cached(Zn60, tfactors.T9, pf_cache, Zn60_log_pf, dZn60_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Cu59_log_pf - Zn60_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dCu59_log_pf_dT9 - dZn60_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  62.31218638341413 + -59.240729222177784 * tfactors.T9i + -39.9808 * tfactors.T913i + 1.17311 * tfactors.T913
                         + -2.90486 * tfactors.T9 + 0.339644 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  59.240729222177784 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -39.9808 * tfactors.T943i + (1.0/3.0) * 1.17311 * tfactors.T923i
                                  + -2.90486 + (5.0/3.0) * 0.339644 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_Zn60_to_He4_Ni56_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Zn60 --> He4 + Ni56

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Ni56_log_pf, dNi56_log_pf_dT9;

    // interpolating Ni56 partition function
    get_partition_function_cached(Ni56, tfactors.T9, pf_cache, Ni56_log_pf, dNi56_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Zn60_log_pf, dZn60_log_pf_dT9;

    // interpolating Zn60 partition function
    get_partition_function_cached(Zn60, tfactors.T9, pf_cache, Zn60_log_pf, dZn60_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Ni56_log_pf - Zn60_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dNi56_log_pf_dT9 - dZn60_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  86.06332792090964 + -31.2372260432927 * tfactors.T9i + -96.4898 * tfactors.T913i + 6.47209 * tfactors.T913
                         + -5.2029 * tfactors.T9 + 0.533391 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  31.2372260432927 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -96.4898 * tfactors.T943i + (1.0/3.0) * 6.47209 * tfactors.T923i
                                  + -5.2029 + (5.0/3.0) * 0.533391 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Mn52_to_n_Fe52_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mn52 + p --> n + Fe52

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Mn52_log_pf, dMn52_log_pf_dT9;

    // interpolating Mn52 partition function
    get_partition_function_cached(Mn52, tfactors.T9, pf_cache, Mn52_log_pf, dMn52_log_pf_dT9);

    amrex::Real Fe52_log_pf, dFe52_log_pf_dT9;

    // interpolating Fe52 partition function
    get_partition_function_cached(Fe52, tfactors.T9, pf_cache, Fe52_log_pf, dFe52_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Fe52_log_pf - p_log_pf - Mn52_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dFe52_log_pf_dT9 - dp_log_pf_dT9 - dMn52_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  15.385756109542921 + -36.68927427806889 * tfactors.T9i + 0.569856 * tfactors.T913
                         + 0.390117 * tfactors.T9 + -0.0479043 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  36.68927427806889 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 0.569856 * tfactors.T923i
                                  + 0.390117 + (5.0/3.0) * -0.0479043 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Mn52_to_n_Co55_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mn52 + He4 --> n + Co55

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Mn52_log_pf, dMn52_log_pf_dT9;

    // interpolating Mn52 partition function
    get_partition_function_cached(Mn52, tfactors.T9, pf_cache, Mn52_log_pf, dMn52_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Co55_log_pf, dCo55_log_pf_dT9;

    // interpolating Co55 partition function
    get_partition_function_cached(Co55, tfactors.T9, pf_cache, Co55_log_pf, dCo55_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Co55_log_pf - He4_log_pf - Mn52_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dCo55_log_pf_dT9 - dHe4_log_pf_dT9 - dMn52_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  -14.314644699269245 + -27.01290711629877 * tfactors.T9i + 8.15876 * tfactors.T913
                         + 3.05978 * tfactors.T9 + -0.406535 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  27.01290711629877 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 8.15876 * tfactors.T923i
                                  + 3.05978 + (5.0/3.0) * -0.406535 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Fe53_to_He4_Mn50_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe53 + p --> He4 + Mn50

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Mn50_log_pf, dMn50_log_pf_dT9;

    // interpolating Mn50 partition function
    get_partition_function_cached(Mn50, tfactors.T9, pf_cache, Mn50_log_pf, dMn50_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Fe53_log_pf, dFe53_log_pf_dT9;

    // interpolating Fe53 partition function
    get_partition_function_cached(Fe53, tfactors.T9, pf_cache, Fe53_log_pf, dFe53_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = He4_log_pf + Mn50_log_pf - p_log_pf - Fe53_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dMn50_log_pf_dT9 - dp_log_pf_dT9 - dFe53_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  61.873981640134645 + -40.111652884372496 * tfactors.T9i + -89.2309 * tfactors.T913i + -1.4758 * tfactors.T913
                         + -0.437831 * tfactors.T9 + -0.00354594 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  40.111652884372496 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -89.2309 * tfactors.T943i + (1.0/3.0) * -1.4758 * tfactors.T923i
                                  + -0.437831 + (5.0/3.0) * -0.00354594 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Fe53_to_n_Ni56_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe53 + He4 --> n + Ni56

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Ni56_log_pf, dNi56_log_pf_dT9;

    // interpolating Ni56 partition function
    get_partition_function_cached(Ni56, tfactors.T9, pf_cache, Ni56_log_pf, dNi56_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Fe53_log_pf, dFe53_log_pf_dT9;

    // interpolating Fe53 partition function
    get_partition_function_cached(Fe53, tfactors.T9, pf_cache, Fe53_log_pf, dFe53_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Ni56_log_pf - He4_log_pf - Fe53_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dNi56_log_pf_dT9 - dHe4_log_pf_dT9 - dFe53_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  -11.207137423641244 + -31.17324290806502 * tfactors.T9i + 6.43259 * tfactors.T913
                         + 2.29591 * tfactors.T9 + -0.254724 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  31.17324290806502 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 6.43259 * tfactors.T923i
                                  + 2.29591 + (5.0/3.0) * -0.254724 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Fe54_to_He4_Mn51_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe54 + p --> He4 + Mn51

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Mn51_log_pf, dMn51_log_pf_dT9;

    // interpolating Mn51 partition function
    get_partition_function_cached(Mn51, tfactors.T9, pf_cache, Mn51_log_pf, dMn51_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Fe54_log_pf, dFe54_log_pf_dT9;

    // interpolating Fe54 partition function
    get_partition_function_cached(Fe54, tfactors.T9, pf_cache, Fe54_log_pf, dFe54_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Mn51_log_pf - p_log_pf - Fe54_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dMn51_log_pf_dT9 - dp_log_pf_dT9 - dFe54_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  65.35241368298804 + -36.52341983603332 * tfactors.T9i + -89.274 * tfactors.T913i + -0.862452 * tfactors.T913
                         + -0.635672 * tfactors.T9 + 0.0196464 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  36.52341983603332 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -89.274 * tfactors.T943i + (1.0/3.0) * -0.862452 * tfactors.T923i
                                  + -0.635672 + (5.0/3.0) * 0.0196464 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Fe54_to_n_Ni57_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe54 + He4 --> n + Ni57

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Ni57_log_pf, dNi57_log_pf_dT9;

    // interpolating Ni57 partition function
    get_partition_function_cached(Ni57, tfactors.T9, pf_cache, Ni57_log_pf, dNi57_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Fe54_log_pf, dFe54_log_pf_dT9;

    // interpolating Fe54 partition function
    get_partition_function_cached(Fe54, tfactors.T9, pf_cache, Fe54_log_pf, dFe54_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Ni57_log_pf - He4_log_pf - Fe54_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dNi57_log_pf_dT9 - dHe4_log_pf_dT9 - dFe54_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  14.633898659996703 + -67.50350779116859 * tfactors.T9i + 0.48562 * tfactors.T913
                         + -0.0249302 * tfactors.T9 + 0.0509605 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  67.50350779116859 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 0.48562 * tfactors.T923i
                                  + -0.0249302 + (5.0/3.0) * 0.0509605 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Fe54_to_p_Co57_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe54 + He4 --> p + Co57

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Co57_log_pf, dCo57_log_pf_dT9;

    // interpolating Co57 partition function
    get_partition_function_cached(Co57, tfactors.T9, pf_cache, Co57_log_pf, dCo57_log_pf_dT9);

    amrex::Real Fe54_log_pf, dFe54_log_pf_dT9;

    // interpolating Fe54 partition function
    get_partition_function_cached(Fe54, tfactors.T9, pf_cache, Fe54_log_pf, dFe54_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Co57_log_pf - He4_log_pf - Fe54_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dCo57_log_pf_dT9 - dHe4_log_pf_dT9 - dFe54_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  -1.404090444334669 + -20.574290677453938 * tfactors.T9i + -38.1133 * tfactors.T913i + 29.3541 * tfactors.T913
                         + -4.75966 * tfactors.T9 + 0.40418 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  20.574290677453938 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -38.1133 * tfactors.T943i + (1.0/3.0) * 29.3541 * tfactors.T923i
                                  + -4.75966 + (5.0/3.0) * 0.40418 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Fe55_to_n_Co55_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe55 + p --> n + Co55

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Co55_log_pf, dCo55_log_pf_dT9;

    // interpolating Co55 partition function
    get_partition_function_cached(Co55, tfactors.T9, pf_cache, Co55_log_pf, dCo55_log_pf_dT9);

    amrex::Real Fe55_log_pf, dFe55_log_pf_dT9;

    // interpolating Fe55 partition function
    get_partition_function_cached(Fe55, tfactors.T9, pf_cache, Fe55_log_pf, dFe55_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Co55_log_pf - p_log_pf - Fe55_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dCo55_log_pf_dT9 - dp_log_pf_dT9 - dFe55_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  21.4351983233521 + -49.13059420131897 * tfactors.T9i + -1.62382 * tfactors.T913
                         + 0.58115 * tfactors.T9 + -0.0537057 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  49.13059420131897 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.62382 * tfactors.T923i
                                  + 0.58115 + (5.0/3.0) * -0.0537057 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Fe55_to_He4_Mn52_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe55 + p --> He4 + Mn52

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Mn52_log_pf, dMn52_log_pf_dT9;

    // interpolating Mn52 partition function
    get_partition_function_cached(Mn52, tfactors.T9, pf_cache, Mn52_log_pf, dMn52_log_pf_dT9);

    amrex::Real Fe55_log_pf, dFe55_log_pf_dT9;

    // interpolating Fe55 partition function
    get_partition_function_cached(Fe55, tfactors.T9, pf_cache, Fe55_log_pf, dFe55_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Mn52_log_pf - p_log_pf - Fe55_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dMn52_log_pf_dT9 - dp_log_pf_dT9 - dFe55_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  59.132043022621346 + -22.1176870850202 * tfactors.T9i + -89.3157 * tfactors.T913i + 5.65639 * tfactors.T913
                         + -2.28952 * tfactors.T9 + 0.186017 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  22.1176870850202 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -89.3157 * tfactors.T943i + (1.0/3.0) * 5.65639 * tfactors.T923i
                                  + -2.28952 + (5.0/3.0) * 0.186017 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Fe55_to_n_Ni58_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe55 + He4 --> n + Ni58

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ni58_log_pf, dNi58_log_pf_dT9;

    // interpolating Ni58 partition function
    get_partition_function_cached(Ni58, tfactors.T9, pf_cache, Ni58_log_pf, dNi58_log_pf_dT9);

    amrex::Real Fe55_log_pf, dFe55_log_pf_dT9;

    // interpolating Fe55 partition function
    get_partition_function_cached(Fe55, tfactors.T9, pf_cache, Fe55_log_pf, dFe55_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Ni58_log_pf - He4_log_pf - Fe55_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dNi58_log_pf_dT9 - dHe4_log_pf_dT9 - dFe55_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  -4.711973544129908 + -33.64036346067392 * tfactors.T9i + 3.44996 * tfactors.T913
                         + 2.98226 * tfactors.T9 + -0.387699 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  33.64036346067392 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 3.44996 * tfactors.T923i
                                  + 2.98226 + (5.0/3.0) * -0.387699 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Fe56_to_n_Co56_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Fe56 + p --> n + Co56

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Fe56_log_pf, dFe56_log_pf_dT9;

    // interpolating Fe56 partition function
    get_partition_function_cached(Fe56, tfactors.T9, pf_cache, Fe56_log_pf, dFe56_log_pf_dT9);

    amrex::Real Co56_log_pf, dCo56_log_pf_dT9;

    // interpolating Co56 partition function
    get_partition_function_cached(Co56, tfactors.T9, pf_cache, Co56_log_pf, dCo56_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Co56_log_pf - p_log_pf - Fe56_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dCo56_log_pf_dT9 - dp_log_pf_dT9 - dFe56_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  22.653426224164285 + -62.07264908160129 * tfactors.T9i + -1.13331 * tfactors.T913
                         + 0.347185 * tfactors.T9 + -0.0328879 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  62.07264908160129 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.13331 * tfactors.T923i
                                  + 0.347185 + (5.0/3.0) * -0.0328879 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Co56_to_n_Ni56_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co56 + p --> n + Ni56

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Ni56_log_pf, dNi56_log_pf_dT9;

    // interpolating Ni56 partition function
    get_partition_function_cached(Ni56, tfactors.T9, pf_cache, Ni56_log_pf, dNi56_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Co56_log_pf, dCo56_log_pf_dT9;

    // interpolating Co56 partition function
    get_partition_function_cached(Co56, tfactors.T9, pf_cache, Co56_log_pf, dCo56_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Ni56_log_pf - p_log_pf - Co56_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dNi56_log_pf_dT9 - dp_log_pf_dT9 - dCo56_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  13.374260291785419 + -33.830037058152215 * tfactors.T9i + 1.76846 * tfactors.T913
                         + 0.197992 * tfactors.T9 + -0.017494 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  33.830037058152215 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 1.76846 * tfactors.T923i
                                  + 0.197992 + (5.0/3.0) * -0.017494 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Co56_to_He4_Fe53_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co56 + p --> He4 + Fe53

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Fe53_log_pf, dFe53_log_pf_dT9;

    // interpolating Fe53 partition function
    get_partition_function_cached(Fe53, tfactors.T9, pf_cache, Fe53_log_pf, dFe53_log_pf_dT9);

    amrex::Real Co56_log_pf, dCo56_log_pf_dT9;

    // interpolating Co56 partition function
    get_partition_function_cached(Co56, tfactors.T9, pf_cache, Co56_log_pf, dCo56_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Fe53_log_pf - p_log_pf - Co56_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dFe53_log_pf_dT9 - dp_log_pf_dT9 - dCo56_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  56.029567715426666 + -2.6567941500871917 * tfactors.T9i + -91.7231 * tfactors.T913i + 10.052 * tfactors.T913
                         + -3.86332 * tfactors.T9 + 0.391105 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  2.6567941500871917 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -91.7231 * tfactors.T943i + (1.0/3.0) * 10.052 * tfactors.T923i
                                  + -3.86332 + (5.0/3.0) * 0.391105 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Co56_to_n_Cu59_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co56 + He4 --> n + Cu59

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Cu59_log_pf, dCu59_log_pf_dT9;

    // interpolating Cu59 partition function
    get_partition_function_cached(Cu59, tfactors.T9, pf_cache, Cu59_log_pf, dCu59_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Co56_log_pf, dCo56_log_pf_dT9;

    // interpolating Co56 partition function
    get_partition_function_cached(Co56, tfactors.T9, pf_cache, Co56_log_pf, dCo56_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Cu59_log_pf - He4_log_pf - Co56_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dCu59_log_pf_dT9 - dHe4_log_pf_dT9 - dCo56_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  8.911918754289905 + -61.8335402370445 * tfactors.T9i + -1.00882 * tfactors.T913
                         + 1.88993 * tfactors.T9 + -0.202319 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  61.8335402370445 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.00882 * tfactors.T923i
                                  + 1.88993 + (5.0/3.0) * -0.202319 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Co57_to_n_Ni57_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Co57 + p --> n + Ni57

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Ni57_log_pf, dNi57_log_pf_dT9;

    // interpolating Ni57 partition function
    get_partition_function_cached(Ni57, tfactors.T9, pf_cache, Ni57_log_pf, dNi57_log_pf_dT9);

    amrex::Real Co57_log_pf, dCo57_log_pf_dT9;

    // interpolating Co57 partition function
    get_partition_function_cached(Co57, tfactors.T9, pf_cache, Co57_log_pf, dCo57_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = n_log_pf + Ni57_log_pf - p_log_pf - Co57_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dNi57_log_pf_dT9 - dp_log_pf_dT9 - dCo57_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  19.30298910433137 + -46.92921711371466 * tfactors.T9i + -1.13177 * tfactors.T913
                         + 0.56875 * tfactors.T9 + -0.0579913 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  46.92921711371466 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.13177 * tfactors.T923i
                                  + 0.56875 + (5.0/3.0) * -0.0579913 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Ni56_to_p_Cu59_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni56 + He4 --> p + Cu59

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Cu59_log_pf, dCu59_log_pf_dT9;

    // interpolating Cu59 partition function
    get_partition_function_cached(Cu59, tfactors.T9, pf_cache, Cu59_log_pf, dCu59_log_pf_dT9);

    amrex::Real Ni56_log_pf, dNi56_log_pf_dT9;

    // interpolating Ni56 partition function
    get_partition_function_cached(Ni56, tfactors.T9, pf_cache, Ni56_log_pf, dNi56_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = p_log_pf + Cu59_log_pf - He4_log_pf - Ni56_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dCu59_log_pf_dT9 - dHe4_log_pf_dT9 - dNi56_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  20.864958462504486 + -28.003503178913892 * tfactors.T9i + -39.9808 * tfactors.T913i + 13.6773 * tfactors.T913
                         + -3.76429 * tfactors.T9 + 0.438096 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  28.003503178913892 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -39.9808 * tfactors.T943i + (1.0/3.0) * 13.6773 * tfactors.T923i
                                  + -3.76429 + (5.0/3.0) * 0.438096 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_He4_Ni57_to_n_Zn60_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni57 + He4 --> n + Zn60

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real n_log_pf, dn_log_pf_dT9;

    // setting n log(partition function) to 0.0 by default, independent of T
    n_log_pf = 0.0_rt;
    dn_log_pf_dT9 = 0.0_rt;

    amrex::Real Ni57_log_pf, dNi57_log_pf_dT9;

    // interpolating Ni57 partition function
    get_partition_function_cached(Ni57, tfactors.T9, pf_cache, Ni57_log_pf, dNi57_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Zn60_log_pf, dZn60_log_pf_dT9;

    // interpolating Zn60 partition function
    get_partition_function_cached(Zn60, tfactors.T9, pf_cache, Zn60_log_pf, dZn60_log_pf_dT9);

    amrex::Real net_log_pf = n_log_pf + Zn60_log_pf - He4_log_pf - Ni57_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dZn60_log_pf_dT9 - dHe4_log_pf_dT9 - dNi57_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  9.748611363087353 + -87.68260435279066 * tfactors.T9i + 5.43552 * tfactors.T913
                         + -1.25169 * tfactors.T9 + 0.148383 * tfactors.T953;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  87.68260435279066 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 5.43552 * tfactors.T923i
                                  + -1.25169 + (5.0/3.0) * 0.148383 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}

template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void rate_p_Ni58_to_He4_Co55_derived(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ni58 + p --> He4 + Co55

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Co55_log_pf, dCo55_log_pf_dT9;

    // interpolating Co55 partition function
    get_partition_function_cached(Co55, tfactors.T9, pf_cache, Co55_log_pf, dCo55_log_pf_dT9);

    amrex::Real Ni58_log_pf, dNi58_log_pf_dT9;

    // interpolating Ni58 partition function
    get_partition_function_cached(Ni58, tfactors.T9, pf_cache, Ni58_log_pf, dNi58_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Co55_log_pf - p_log_pf - Ni58_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCo55_log_pf_dT9 - dp_log_pf_dT9 - dNi58_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  63.596461867482006 + -15.490230740645051 * tfactors.T9i + -94.1404 * tfactors.T913i + 3.39179 * tfactors.T913
                         + -1.71062 * tfactors.T9 + 0.133003 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf + log_scor;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  15.490230740645051 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -94.1404 * tfactors.T943i + (1.0/3.0) * 3.39179 * tfactors.T923i
                                  + -1.71062 + (5.0/3.0) * 0.133003 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9 + dlog_scor_dT * 1.0e9_rt;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

}


template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
fill_derived_rates(const tf_t& tfactors, T& rate_eval)
{

    amrex::ignore_unused(tfactors);
    amrex::ignore_unused(rate_eval);

    [[maybe_unused]] amrex::Real rate;
    [[maybe_unused]] amrex::Real drate_dT;

    part_fun::pf_cache_t pf_cache{};

    pf_cache.index_temp_array_1 = interp_net::find_index(tfactors.T9, part_fun::temp_array_1);

    {
        // N13_to_p_C12_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_N13_to_p_C12_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_N13_to_p_C12_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_N13_to_p_C12_derived) = drate_dT;
        }
    }

    {
        // N14_to_p_C13_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_N14_to_p_C13_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_N14_to_p_C13_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_N14_to_p_C13_derived) = drate_dT;
        }
    }

    {
        // N15_to_p_C14_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_N15_to_p_C14_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_N15_to_p_C14_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_N15_to_p_C14_derived) = drate_dT;
        }
    }

    {
        // O16_to_p_N15_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_O16_to_p_N15_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_O16_to_p_N15_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_O16_to_p_N15_derived) = drate_dT;
        }
    }

    {
        // O16_to_He4_C12_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_O16_to_He4_C12_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_O16_to_He4_C12_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_O16_to_He4_C12_derived) = drate_dT;
        }
    }

    {
        // O18_to_He4_C14_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_O18_to_He4_C14_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_O18_to_He4_C14_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_O18_to_He4_C14_derived) = drate_dT;
        }
    }

    {
        // F18_to_p_O17_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_F18_to_p_O17_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_F18_to_p_O17_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_F18_to_p_O17_derived) = drate_dT;
        }
    }

    {
        // F18_to_He4_N14_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_F18_to_He4_N14_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_F18_to_He4_N14_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_F18_to_He4_N14_derived) = drate_dT;
        }
    }

    {
        // F19_to_p_O18_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_F19_to_p_O18_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_F19_to_p_O18_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_F19_to_p_O18_derived) = drate_dT;
        }
    }

    {
        // F19_to_He4_N15_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_F19_to_He4_N15_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_F19_to_He4_N15_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_F19_to_He4_N15_derived) = drate_dT;
        }
    }

    {
        // Ne20_to_p_F19_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ne20_to_p_F19_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ne20_to_p_F19_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne20_to_p_F19_derived) = drate_dT;
        }
    }

    {
        // Ne20_to_He4_O16_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ne20_to_He4_O16_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ne20_to_He4_O16_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne20_to_He4_O16_derived) = drate_dT;
        }
    }

    {
        // Ne21_to_He4_O17_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ne21_to_He4_O17_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ne21_to_He4_O17_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne21_to_He4_O17_derived) = drate_dT;
        }
    }

    {
        // Ne22_to_He4_O18_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ne22_to_He4_O18_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ne22_to_He4_O18_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne22_to_He4_O18_derived) = drate_dT;
        }
    }

    {
        // Na22_to_p_Ne21_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Na22_to_p_Ne21_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Na22_to_p_Ne21_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Na22_to_p_Ne21_derived) = drate_dT;
        }
    }

    {
        // Na22_to_He4_F18_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Na22_to_He4_F18_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Na22_to_He4_F18_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Na22_to_He4_F18_derived) = drate_dT;
        }
    }

    {
        // Na23_to_p_Ne22_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Na23_to_p_Ne22_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Na23_to_p_Ne22_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Na23_to_p_Ne22_derived) = drate_dT;
        }
    }

    {
        // Na23_to_He4_F19_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Na23_to_He4_F19_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Na23_to_He4_F19_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Na23_to_He4_F19_derived) = drate_dT;
        }
    }

    {
        // Mg24_to_p_Na23_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Mg24_to_p_Na23_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Mg24_to_p_Na23_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg24_to_p_Na23_derived) = drate_dT;
        }
    }

    {
        // Mg24_to_He4_Ne20_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Mg24_to_He4_Ne20_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Mg24_to_He4_Ne20_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg24_to_He4_Ne20_derived) = drate_dT;
        }
    }

    {
        // Mg25_to_He4_Ne21_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Mg25_to_He4_Ne21_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Mg25_to_He4_Ne21_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg25_to_He4_Ne21_derived) = drate_dT;
        }
    }

    {
        // Al26_to_p_Mg25_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Al26_to_p_Mg25_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Al26_to_p_Mg25_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Al26_to_p_Mg25_derived) = drate_dT;
        }
    }

    {
        // Al26_to_He4_Na22_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Al26_to_He4_Na22_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Al26_to_He4_Na22_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Al26_to_He4_Na22_derived) = drate_dT;
        }
    }

    {
        // Al27_to_He4_Na23_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Al27_to_He4_Na23_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Al27_to_He4_Na23_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Al27_to_He4_Na23_derived) = drate_dT;
        }
    }

    {
        // Si28_to_p_Al27_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Si28_to_p_Al27_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Si28_to_p_Al27_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si28_to_p_Al27_derived) = drate_dT;
        }
    }

    {
        // Si28_to_He4_Mg24_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Si28_to_He4_Mg24_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Si28_to_He4_Mg24_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si28_to_He4_Mg24_derived) = drate_dT;
        }
    }

    {
        // Si29_to_He4_Mg25_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Si29_to_He4_Mg25_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Si29_to_He4_Mg25_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si29_to_He4_Mg25_derived) = drate_dT;
        }
    }

    {
        // P30_to_p_Si29_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_P30_to_p_Si29_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_P30_to_p_Si29_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P30_to_p_Si29_derived) = drate_dT;
        }
    }

    {
        // P30_to_He4_Al26_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_P30_to_He4_Al26_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_P30_to_He4_Al26_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P30_to_He4_Al26_derived) = drate_dT;
        }
    }

    {
        // P31_to_He4_Al27_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_P31_to_He4_Al27_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_P31_to_He4_Al27_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P31_to_He4_Al27_derived) = drate_dT;
        }
    }

    {
        // S32_to_p_P31_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_S32_to_p_P31_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_S32_to_p_P31_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S32_to_p_P31_derived) = drate_dT;
        }
    }

    {
        // S32_to_He4_Si28_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_S32_to_He4_Si28_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_S32_to_He4_Si28_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S32_to_He4_Si28_derived) = drate_dT;
        }
    }

    {
        // S33_to_He4_Si29_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_S33_to_He4_Si29_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_S33_to_He4_Si29_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S33_to_He4_Si29_derived) = drate_dT;
        }
    }

    {
        // Cl34_to_p_S33_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Cl34_to_p_S33_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Cl34_to_p_S33_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl34_to_p_S33_derived) = drate_dT;
        }
    }

    {
        // Cl34_to_He4_P30_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Cl34_to_He4_P30_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Cl34_to_He4_P30_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl34_to_He4_P30_derived) = drate_dT;
        }
    }

    {
        // Cl35_to_He4_P31_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Cl35_to_He4_P31_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Cl35_to_He4_P31_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl35_to_He4_P31_derived) = drate_dT;
        }
    }

    {
        // Ar36_to_p_Cl35_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ar36_to_p_Cl35_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ar36_to_p_Cl35_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar36_to_p_Cl35_derived) = drate_dT;
        }
    }

    {
        // Ar36_to_He4_S32_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ar36_to_He4_S32_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ar36_to_He4_S32_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar36_to_He4_S32_derived) = drate_dT;
        }
    }

    {
        // Ar37_to_He4_S33_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ar37_to_He4_S33_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ar37_to_He4_S33_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar37_to_He4_S33_derived) = drate_dT;
        }
    }

    {
        // K38_to_p_Ar37_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_K38_to_p_Ar37_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_K38_to_p_Ar37_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K38_to_p_Ar37_derived) = drate_dT;
        }
    }

    {
        // K38_to_He4_Cl34_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_K38_to_He4_Cl34_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_K38_to_He4_Cl34_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K38_to_He4_Cl34_derived) = drate_dT;
        }
    }

    {
        // K39_to_He4_Cl35_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_K39_to_He4_Cl35_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_K39_to_He4_Cl35_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K39_to_He4_Cl35_derived) = drate_dT;
        }
    }

    {
        // Ca40_to_p_K39_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ca40_to_p_K39_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ca40_to_p_K39_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca40_to_p_K39_derived) = drate_dT;
        }
    }

    {
        // Ca40_to_He4_Ar36_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ca40_to_He4_Ar36_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ca40_to_He4_Ar36_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca40_to_He4_Ar36_derived) = drate_dT;
        }
    }

    {
        // Ca41_to_He4_Ar37_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ca41_to_He4_Ar37_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ca41_to_He4_Ar37_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca41_to_He4_Ar37_derived) = drate_dT;
        }
    }

    {
        // Sc42_to_p_Ca41_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Sc42_to_p_Ca41_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Sc42_to_p_Ca41_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc42_to_p_Ca41_derived) = drate_dT;
        }
    }

    {
        // Sc42_to_He4_K38_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Sc42_to_He4_K38_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Sc42_to_He4_K38_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc42_to_He4_K38_derived) = drate_dT;
        }
    }

    {
        // Sc43_to_He4_K39_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Sc43_to_He4_K39_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Sc43_to_He4_K39_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc43_to_He4_K39_derived) = drate_dT;
        }
    }

    {
        // Ti44_to_p_Sc43_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ti44_to_p_Sc43_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ti44_to_p_Sc43_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti44_to_p_Sc43_derived) = drate_dT;
        }
    }

    {
        // Ti44_to_He4_Ca40_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ti44_to_He4_Ca40_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ti44_to_He4_Ca40_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti44_to_He4_Ca40_derived) = drate_dT;
        }
    }

    {
        // Ti45_to_He4_Ca41_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ti45_to_He4_Ca41_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ti45_to_He4_Ca41_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti45_to_He4_Ca41_derived) = drate_dT;
        }
    }

    {
        // V46_to_p_Ti45_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_V46_to_p_Ti45_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_V46_to_p_Ti45_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V46_to_p_Ti45_derived) = drate_dT;
        }
    }

    {
        // V46_to_He4_Sc42_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_V46_to_He4_Sc42_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_V46_to_He4_Sc42_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V46_to_He4_Sc42_derived) = drate_dT;
        }
    }

    {
        // V47_to_He4_Sc43_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_V47_to_He4_Sc43_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_V47_to_He4_Sc43_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V47_to_He4_Sc43_derived) = drate_dT;
        }
    }

    {
        // Cr48_to_p_V47_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Cr48_to_p_V47_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Cr48_to_p_V47_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr48_to_p_V47_derived) = drate_dT;
        }
    }

    {
        // Cr48_to_He4_Ti44_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Cr48_to_He4_Ti44_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Cr48_to_He4_Ti44_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr48_to_He4_Ti44_derived) = drate_dT;
        }
    }

    {
        // Cr49_to_He4_Ti45_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Cr49_to_He4_Ti45_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Cr49_to_He4_Ti45_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr49_to_He4_Ti45_derived) = drate_dT;
        }
    }

    {
        // Mn50_to_p_Cr49_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Mn50_to_p_Cr49_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Mn50_to_p_Cr49_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn50_to_p_Cr49_derived) = drate_dT;
        }
    }

    {
        // Mn50_to_He4_V46_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Mn50_to_He4_V46_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Mn50_to_He4_V46_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn50_to_He4_V46_derived) = drate_dT;
        }
    }

    {
        // Mn51_to_He4_V47_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Mn51_to_He4_V47_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Mn51_to_He4_V47_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn51_to_He4_V47_derived) = drate_dT;
        }
    }

    {
        // Fe52_to_p_Mn51_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Fe52_to_p_Mn51_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Fe52_to_p_Mn51_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe52_to_p_Mn51_derived) = drate_dT;
        }
    }

    {
        // Fe52_to_He4_Cr48_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Fe52_to_He4_Cr48_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Fe52_to_He4_Cr48_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe52_to_He4_Cr48_derived) = drate_dT;
        }
    }

    {
        // Co55_to_He4_Mn51_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Co55_to_He4_Mn51_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Co55_to_He4_Mn51_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co55_to_He4_Mn51_derived) = drate_dT;
        }
    }

    {
        // Ni56_to_p_Co55_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ni56_to_p_Co55_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ni56_to_p_Co55_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni56_to_p_Co55_derived) = drate_dT;
        }
    }

    {
        // Ni56_to_He4_Fe52_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ni56_to_He4_Fe52_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ni56_to_He4_Fe52_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni56_to_He4_Fe52_derived) = drate_dT;
        }
    }

    {
        // C12_to_He4_He4_He4_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_C12_to_He4_He4_He4_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_C12_to_He4_He4_He4_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_C12_to_He4_He4_He4_derived) = drate_dT;
        }
    }

    {
        // He4_C12_to_p_N15_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_C12);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_C12);
        }
#endif
        rate_He4_C12_to_p_N15_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_C12_to_p_N15_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_C12_to_p_N15_derived) = drate_dT;
        }
    }

    {
        // He4_N14_to_p_O17_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_N14);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_N14);
        }
#endif
        rate_He4_N14_to_p_O17_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_N14_to_p_O17_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_N14_to_p_O17_derived) = drate_dT;
        }
    }

    {
        // He4_N15_to_p_O18_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_N15);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_N15);
        }
#endif
        rate_He4_N15_to_p_O18_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_N15_to_p_O18_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_N15_to_p_O18_derived) = drate_dT;
        }
    }

    {
        // p_O16_to_He4_N13_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_O16);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_O16);
        }
#endif
        rate_p_O16_to_He4_N13_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_O16_to_He4_N13_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_O16_to_He4_N13_derived) = drate_dT;
        }
    }

    {
        // He4_O16_to_p_F19_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_O16);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_O16);
        }
#endif
        rate_He4_O16_to_p_F19_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_O16_to_p_F19_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_O16_to_p_F19_derived) = drate_dT;
        }
    }

    {
        // He4_Ne20_to_p_Na23_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Ne20);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ne20);
        }
#endif
        rate_He4_Ne20_to_p_Na23_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Ne20_to_p_Na23_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Ne20_to_p_Na23_derived) = drate_dT;
        }
    }

    {
        // He4_Ne20_to_C12_C12_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Ne20);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ne20);
        }
#endif
        rate_He4_Ne20_to_C12_C12_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Ne20_to_C12_C12_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Ne20_to_C12_C12_derived) = drate_dT;
        }
    }

    {
        // p_Ne21_to_He4_F18_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Ne21);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ne21);
        }
#endif
        rate_p_Ne21_to_He4_F18_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Ne21_to_He4_F18_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Ne21_to_He4_F18_derived) = drate_dT;
        }
    }

    {
        // p_Ne22_to_He4_F19_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Ne22);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ne22);
        }
#endif
        rate_p_Ne22_to_He4_F19_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Ne22_to_He4_F19_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Ne22_to_He4_F19_derived) = drate_dT;
        }
    }

    {
        // p_Na23_to_C12_C12_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Na23);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Na23);
        }
#endif
        rate_p_Na23_to_C12_C12_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Na23_to_C12_C12_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Na23_to_C12_C12_derived) = drate_dT;
        }
    }

    {
        // He4_Mg24_to_p_Al27_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Mg24);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mg24);
        }
#endif
        rate_He4_Mg24_to_p_Al27_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Mg24_to_p_Al27_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Mg24_to_p_Al27_derived) = drate_dT;
        }
    }

    {
        // He4_Mg24_to_C12_O16_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Mg24);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mg24);
        }
#endif
        rate_He4_Mg24_to_C12_O16_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Mg24_to_C12_O16_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Mg24_to_C12_O16_derived) = drate_dT;
        }
    }

    {
        // p_Mg25_to_He4_Na22_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Mg25);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mg25);
        }
#endif
        rate_p_Mg25_to_He4_Na22_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Mg25_to_He4_Na22_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Mg25_to_He4_Na22_derived) = drate_dT;
        }
    }

    {
        // p_Al27_to_C12_O16_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Al27);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Al27);
        }
#endif
        rate_p_Al27_to_C12_O16_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Al27_to_C12_O16_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Al27_to_C12_O16_derived) = drate_dT;
        }
    }

    {
        // He4_Si28_to_p_P31_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Si28);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si28);
        }
#endif
        rate_He4_Si28_to_p_P31_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Si28_to_p_P31_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Si28_to_p_P31_derived) = drate_dT;
        }
    }

    {
        // He4_Si28_to_C12_Ne20_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Si28);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si28);
        }
#endif
        rate_He4_Si28_to_C12_Ne20_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Si28_to_C12_Ne20_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Si28_to_C12_Ne20_derived) = drate_dT;
        }
    }

    {
        // He4_Si28_to_O16_O16_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Si28);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si28);
        }
#endif
        rate_He4_Si28_to_O16_O16_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Si28_to_O16_O16_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Si28_to_O16_O16_derived) = drate_dT;
        }
    }

    {
        // p_Si29_to_He4_Al26_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Si29);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Si29);
        }
#endif
        rate_p_Si29_to_He4_Al26_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Si29_to_He4_Al26_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Si29_to_He4_Al26_derived) = drate_dT;
        }
    }

    {
        // p_P31_to_C12_Ne20_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_P31);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_P31);
        }
#endif
        rate_p_P31_to_C12_Ne20_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_P31_to_C12_Ne20_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_P31_to_C12_Ne20_derived) = drate_dT;
        }
    }

    {
        // p_P31_to_O16_O16_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_P31);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_P31);
        }
#endif
        rate_p_P31_to_O16_O16_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_P31_to_O16_O16_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_P31_to_O16_O16_derived) = drate_dT;
        }
    }

    {
        // He4_S32_to_p_Cl35_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_S32);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_S32);
        }
#endif
        rate_He4_S32_to_p_Cl35_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_S32_to_p_Cl35_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_S32_to_p_Cl35_derived) = drate_dT;
        }
    }

    {
        // p_S33_to_He4_P30_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_S33);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_S33);
        }
#endif
        rate_p_S33_to_He4_P30_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_S33_to_He4_P30_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_S33_to_He4_P30_derived) = drate_dT;
        }
    }

    {
        // He4_Ar36_to_p_K39_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Ar36);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ar36);
        }
#endif
        rate_He4_Ar36_to_p_K39_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Ar36_to_p_K39_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Ar36_to_p_K39_derived) = drate_dT;
        }
    }

    {
        // p_Ar37_to_He4_Cl34_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Ar37);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ar37);
        }
#endif
        rate_p_Ar37_to_He4_Cl34_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Ar37_to_He4_Cl34_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Ar37_to_He4_Cl34_derived) = drate_dT;
        }
    }

    {
        // He4_Ca40_to_p_Sc43_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Ca40);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca40);
        }
#endif
        rate_He4_Ca40_to_p_Sc43_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Ca40_to_p_Sc43_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Ca40_to_p_Sc43_derived) = drate_dT;
        }
    }

    {
        // p_Ca41_to_He4_K38_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Ca41);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca41);
        }
#endif
        rate_p_Ca41_to_He4_K38_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Ca41_to_He4_K38_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Ca41_to_He4_K38_derived) = drate_dT;
        }
    }

    {
        // p_Ti45_to_He4_Sc42_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Ti45);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti45);
        }
#endif
        rate_p_Ti45_to_He4_Sc42_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Ti45_to_He4_Sc42_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Ti45_to_He4_Sc42_derived) = drate_dT;
        }
    }

    {
        // p_V47_to_He4_Ti44_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_V47);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V47);
        }
#endif
        rate_p_V47_to_He4_Ti44_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_V47_to_He4_Ti44_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_V47_to_He4_Ti44_derived) = drate_dT;
        }
    }

    {
        // p_Cr49_to_He4_V46_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Cr49);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr49);
        }
#endif
        rate_p_Cr49_to_He4_V46_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Cr49_to_He4_V46_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Cr49_to_He4_V46_derived) = drate_dT;
        }
    }

    {
        // p_Mn51_to_He4_Cr48_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Mn51);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn51);
        }
#endif
        rate_p_Mn51_to_He4_Cr48_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Mn51_to_He4_Cr48_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Mn51_to_He4_Cr48_derived) = drate_dT;
        }
    }

    {
        // p_Mn52_to_He4_Cr49_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Mn52);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn52);
        }
#endif
        rate_p_Mn52_to_He4_Cr49_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Mn52_to_He4_Cr49_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Mn52_to_He4_Cr49_derived) = drate_dT;
        }
    }

    {
        // p_Co55_to_He4_Fe52_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Co55);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co55);
        }
#endif
        rate_p_Co55_to_He4_Fe52_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Co55_to_He4_Fe52_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Co55_to_He4_Fe52_derived) = drate_dT;
        }
    }

    {
        // Mn51_to_n_Mn50_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Mn51_to_n_Mn50_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Mn51_to_n_Mn50_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn51_to_n_Mn50_derived) = drate_dT;
        }
    }

    {
        // Mn52_to_n_Mn51_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Mn52_to_n_Mn51_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Mn52_to_n_Mn51_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn52_to_n_Mn51_derived) = drate_dT;
        }
    }

    {
        // Fe53_to_n_Fe52_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Fe53_to_n_Fe52_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Fe53_to_n_Fe52_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe53_to_n_Fe52_derived) = drate_dT;
        }
    }

    {
        // Fe53_to_p_Mn52_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Fe53_to_p_Mn52_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Fe53_to_p_Mn52_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe53_to_p_Mn52_derived) = drate_dT;
        }
    }

    {
        // Fe54_to_n_Fe53_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Fe54_to_n_Fe53_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Fe54_to_n_Fe53_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe54_to_n_Fe53_derived) = drate_dT;
        }
    }

    {
        // Fe55_to_n_Fe54_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Fe55_to_n_Fe54_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Fe55_to_n_Fe54_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe55_to_n_Fe54_derived) = drate_dT;
        }
    }

    {
        // Fe56_to_n_Fe55_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Fe56_to_n_Fe55_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Fe56_to_n_Fe55_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe56_to_n_Fe55_derived) = drate_dT;
        }
    }

    {
        // Co55_to_p_Fe54_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Co55_to_p_Fe54_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Co55_to_p_Fe54_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co55_to_p_Fe54_derived) = drate_dT;
        }
    }

    {
        // Co56_to_n_Co55_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Co56_to_n_Co55_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Co56_to_n_Co55_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co56_to_n_Co55_derived) = drate_dT;
        }
    }

    {
        // Co56_to_p_Fe55_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Co56_to_p_Fe55_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Co56_to_p_Fe55_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co56_to_p_Fe55_derived) = drate_dT;
        }
    }

    {
        // Co56_to_He4_Mn52_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Co56_to_He4_Mn52_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Co56_to_He4_Mn52_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co56_to_He4_Mn52_derived) = drate_dT;
        }
    }

    {
        // Co57_to_n_Co56_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Co57_to_n_Co56_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Co57_to_n_Co56_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co57_to_n_Co56_derived) = drate_dT;
        }
    }

    {
        // Co57_to_p_Fe56_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Co57_to_p_Fe56_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Co57_to_p_Fe56_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co57_to_p_Fe56_derived) = drate_dT;
        }
    }

    {
        // Ni57_to_n_Ni56_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ni57_to_n_Ni56_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ni57_to_n_Ni56_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni57_to_n_Ni56_derived) = drate_dT;
        }
    }

    {
        // Ni57_to_p_Co56_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ni57_to_p_Co56_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ni57_to_p_Co56_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni57_to_p_Co56_derived) = drate_dT;
        }
    }

    {
        // Ni57_to_He4_Fe53_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ni57_to_He4_Fe53_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ni57_to_He4_Fe53_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni57_to_He4_Fe53_derived) = drate_dT;
        }
    }

    {
        // Ni58_to_n_Ni57_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ni58_to_n_Ni57_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ni58_to_n_Ni57_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni58_to_n_Ni57_derived) = drate_dT;
        }
    }

    {
        // Ni58_to_p_Co57_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ni58_to_p_Co57_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ni58_to_p_Co57_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni58_to_p_Co57_derived) = drate_dT;
        }
    }

    {
        // Ni58_to_He4_Fe54_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Ni58_to_He4_Fe54_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Ni58_to_He4_Fe54_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni58_to_He4_Fe54_derived) = drate_dT;
        }
    }

    {
        // Cu59_to_p_Ni58_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Cu59_to_p_Ni58_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Cu59_to_p_Ni58_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu59_to_p_Ni58_derived) = drate_dT;
        }
    }

    {
        // Cu59_to_He4_Co55_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Cu59_to_He4_Co55_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Cu59_to_He4_Co55_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu59_to_He4_Co55_derived) = drate_dT;
        }
    }

    {
        // Zn60_to_p_Cu59_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Zn60_to_p_Cu59_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Zn60_to_p_Cu59_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn60_to_p_Cu59_derived) = drate_dT;
        }
    }

    {
        // Zn60_to_He4_Ni56_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
        rate_Zn60_to_He4_Ni56_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_Zn60_to_He4_Ni56_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn60_to_He4_Ni56_derived) = drate_dT;
        }
    }

    {
        // p_Mn52_to_n_Fe52_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Mn52);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn52);
        }
#endif
        rate_p_Mn52_to_n_Fe52_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Mn52_to_n_Fe52_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Mn52_to_n_Fe52_derived) = drate_dT;
        }
    }

    {
        // He4_Mn52_to_n_Co55_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Mn52);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn52);
        }
#endif
        rate_He4_Mn52_to_n_Co55_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Mn52_to_n_Co55_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Mn52_to_n_Co55_derived) = drate_dT;
        }
    }

    {
        // p_Fe53_to_He4_Mn50_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Fe53);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe53);
        }
#endif
        rate_p_Fe53_to_He4_Mn50_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Fe53_to_He4_Mn50_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Fe53_to_He4_Mn50_derived) = drate_dT;
        }
    }

    {
        // He4_Fe53_to_n_Ni56_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Fe53);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe53);
        }
#endif
        rate_He4_Fe53_to_n_Ni56_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Fe53_to_n_Ni56_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Fe53_to_n_Ni56_derived) = drate_dT;
        }
    }

    {
        // p_Fe54_to_He4_Mn51_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Fe54);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe54);
        }
#endif
        rate_p_Fe54_to_He4_Mn51_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Fe54_to_He4_Mn51_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Fe54_to_He4_Mn51_derived) = drate_dT;
        }
    }

    {
        // He4_Fe54_to_n_Ni57_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Fe54);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe54);
        }
#endif
        rate_He4_Fe54_to_n_Ni57_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Fe54_to_n_Ni57_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Fe54_to_n_Ni57_derived) = drate_dT;
        }
    }

    {
        // He4_Fe54_to_p_Co57_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Fe54);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe54);
        }
#endif
        rate_He4_Fe54_to_p_Co57_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Fe54_to_p_Co57_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Fe54_to_p_Co57_derived) = drate_dT;
        }
    }

    {
        // p_Fe55_to_n_Co55_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Fe55);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe55);
        }
#endif
        rate_p_Fe55_to_n_Co55_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Fe55_to_n_Co55_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Fe55_to_n_Co55_derived) = drate_dT;
        }
    }

    {
        // p_Fe55_to_He4_Mn52_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Fe55);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe55);
        }
#endif
        rate_p_Fe55_to_He4_Mn52_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Fe55_to_He4_Mn52_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Fe55_to_He4_Mn52_derived) = drate_dT;
        }
    }

    {
        // He4_Fe55_to_n_Ni58_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Fe55);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe55);
        }
#endif
        rate_He4_Fe55_to_n_Ni58_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Fe55_to_n_Ni58_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Fe55_to_n_Ni58_derived) = drate_dT;
        }
    }

    {
        // p_Fe56_to_n_Co56_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Fe56);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe56);
        }
#endif
        rate_p_Fe56_to_n_Co56_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Fe56_to_n_Co56_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Fe56_to_n_Co56_derived) = drate_dT;
        }
    }

    {
        // p_Co56_to_n_Ni56_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Co56);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co56);
        }
#endif
        rate_p_Co56_to_n_Ni56_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Co56_to_n_Ni56_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Co56_to_n_Ni56_derived) = drate_dT;
        }
    }

    {
        // p_Co56_to_He4_Fe53_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Co56);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co56);
        }
#endif
        rate_p_Co56_to_He4_Fe53_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Co56_to_He4_Fe53_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Co56_to_He4_Fe53_derived) = drate_dT;
        }
    }

    {
        // He4_Co56_to_n_Cu59_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Co56);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co56);
        }
#endif
        rate_He4_Co56_to_n_Cu59_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Co56_to_n_Cu59_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Co56_to_n_Cu59_derived) = drate_dT;
        }
    }

    {
        // p_Co57_to_n_Ni57_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Co57);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co57);
        }
#endif
        rate_p_Co57_to_n_Ni57_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Co57_to_n_Ni57_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Co57_to_n_Ni57_derived) = drate_dT;
        }
    }

    {
        // He4_Ni56_to_p_Cu59_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Ni56);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni56);
        }
#endif
        rate_He4_Ni56_to_p_Cu59_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Ni56_to_p_Cu59_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Ni56_to_p_Cu59_derived) = drate_dT;
        }
    }

    {
        // He4_Ni57_to_n_Zn60_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_He4_Ni57);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni57);
        }
#endif
        rate_He4_Ni57_to_n_Zn60_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_He4_Ni57_to_n_Zn60_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_Ni57_to_n_Zn60_derived) = drate_dT;
        }
    }

    {
        // p_Ni58_to_He4_Co55_derived

        amrex::Real log_scor {0.0_rt};
        amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
        log_scor = rate_eval.log_screen(k_p_Ni58);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni58);
        }
#endif
        rate_p_Ni58_to_He4_Co55_derived<do_T_derivatives, T>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT, rate_eval, pf_cache);
        rate_eval.screened_rates(k_p_Ni58_to_He4_Co55_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_Ni58_to_He4_Co55_derived) = drate_dT;
        }
    }


}
#endif
