#ifndef DERIVED_RATES_H
#define DERIVED_RATES_H

#include <tfactors.H>
#include <actual_network.H>
#include <partition_functions.H>
#include <temperature_table_rates.H>
#include <interp_tools.H>

using namespace Rates;

template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_N13_to_p_C12_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  22.553277271250543 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -13.692 * tfactors.T943i + (1.0/3.0) * -0.230881 * tfactors.T923i
                                  + 4.44362 + (5.0/3.0) * -3.15898 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  26.331767271250545 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -5.10735 * tfactors.T943i + (1.0/3.0) * -2.24111 * tfactors.T923i
                                  + 0.148883;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_N14_to_p_C13_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  87.62065170634826 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -13.72 * tfactors.T943i + (1.0/3.0) * -0.450018 * tfactors.T923i
                                  + 3.70823 + (5.0/3.0) * -1.70545 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  93.40212170634825 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.196703 * tfactors.T923i
                                  + 0.142126 + (5.0/3.0) * -0.0238912 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  101.17495170634825 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_O14_to_p_N13_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // O14 --> p + N13

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real O14_log_pf, dO14_log_pf_dT9;

    // interpolating O14 partition function
    get_partition_function_cached(O14, tfactors.T9, pf_cache, O14_log_pf, dO14_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real N13_log_pf, dN13_log_pf_dT9;

    // setting N13 log(partition function) to 0.0 by default, independent of T
    N13_log_pf = 0.0_rt;
    dN13_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = p_log_pf + N13_log_pf - O14_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dN13_log_pf_dT9 - dO14_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from lg06r
    ln_set_rate =  35.3038971632548 + -59.81629660012578 * tfactors.T9i + 1.57122 * tfactors.T913i;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  59.81629660012578 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 1.57122 * tfactors.T943i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from lg06n
    ln_set_rate =  42.44239716325481 + -53.690276600125785 * tfactors.T9i + -15.1676 * tfactors.T913i + 0.0955166 * tfactors.T913
                         + 3.0659 * tfactors.T9 + -0.507339 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  53.690276600125785 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -15.1676 * tfactors.T943i + (1.0/3.0) * 0.0955166 * tfactors.T923i
                                  + 3.0659 + (5.0/3.0) * -0.507339 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_O15_to_p_N14_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // O15 --> p + N14

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real O15_log_pf, dO15_log_pf_dT9;

    // interpolating O15 partition function
    get_partition_function_cached(O15, tfactors.T9, pf_cache, O15_log_pf, dO15_log_pf_dT9);

    amrex::Real N14_log_pf, dN14_log_pf_dT9;

    // setting N14 log(partition function) to 0.0 by default, independent of T
    N14_log_pf = 0.0_rt;
    dN14_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = p_log_pf + N14_log_pf - O15_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dN14_log_pf_dT9 - dO15_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from im05r
    ln_set_rate =  30.76303704754867 + -89.56670699689953 * tfactors.T9i
                         + 1.5682 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  89.56670699689953 * tfactors.T9i * tfactors.T9i
                                  + 1.5682 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from im05r
    ln_set_rate =  31.68169704754867 + -87.67370699689953 * tfactors.T9i;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  87.67370699689953 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from im05n
    ln_set_rate =  44.144157047548674 + -84.67570699689952 * tfactors.T9i + -15.193 * tfactors.T913i + -4.63975 * tfactors.T913
                         + 9.73458 * tfactors.T9 + -9.55051 * tfactors.T953 + 1.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  84.67570699689952 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -15.193 * tfactors.T943i + (1.0/3.0) * -4.63975 * tfactors.T923i
                                  + 9.73458 + (5.0/3.0) * -9.55051 * tfactors.T923 + 1.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from im05n
    ln_set_rate =  41.03725704754867 + -84.67570699689952 * tfactors.T9i + -15.193 * tfactors.T913i + -0.161954 * tfactors.T913
                         + -7.52123 * tfactors.T9 + -0.987565 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  84.67570699689952 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -15.193 * tfactors.T943i + (1.0/3.0) * -0.161954 * tfactors.T923i
                                  + -7.52123 + (5.0/3.0) * -0.987565 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_O16_to_p_N15_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  150.96226378057287 * tfactors.T9i * tfactors.T9i
                                  + 0.0459037;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  143.65591378057286 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  140.73276378057287 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -15.24 * tfactors.T943i + (1.0/3.0) * 0.334926 * tfactors.T923i
                                  + 4.59088 + (5.0/3.0) * -4.78468 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_O16_to_He4_C12_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  84.95157686791683 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 103.411 * tfactors.T943i + (1.0/3.0) * -420.567 * tfactors.T923i
                                  + 64.0874 + (5.0/3.0) * -12.4624 * tfactors.T923 + 138.803 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  84.50314686791683 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 58.9128 * tfactors.T943i + (1.0/3.0) * -148.273 * tfactors.T923i
                                  + 9.08324 + (5.0/3.0) * -0.541041 * tfactors.T923 + 71.8554 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_F17_to_p_O16_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // F17 --> p + O16

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real F17_log_pf, dF17_log_pf_dT9;

    // interpolating F17 partition function
    get_partition_function_cached(F17, tfactors.T9, pf_cache, F17_log_pf, dF17_log_pf_dT9);

    amrex::Real O16_log_pf, dO16_log_pf_dT9;

    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors.T9, pf_cache, O16_log_pf, dO16_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + O16_log_pf - F17_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dO16_log_pf_dT9 - dF17_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ia08n
    ln_set_rate =  40.93184403787935 + -6.965832070525503 * tfactors.T9i + -16.696 * tfactors.T913i + -1.16252 * tfactors.T913
                         + 0.267703 * tfactors.T9 + -0.0338411 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  6.965832070525503 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -16.696 * tfactors.T943i + (1.0/3.0) * -1.16252 * tfactors.T923i
                                  + 0.267703 + (5.0/3.0) * -0.0338411 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_F18_to_p_O17_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  71.29605321275191 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 2.31435 * tfactors.T923i
                                  + -0.302835 + (5.0/3.0) * 0.020133 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  65.81406921275192 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  65.06777321275192 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -16.4035 * tfactors.T943i + (1.0/3.0) * 4.31885 * tfactors.T923i
                                  + -0.709921 + (5.0/3.0) * -2.0 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_F18_to_He4_N14_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  62.20224752987261 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -5.6227 * tfactors.T943i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  56.396987529872604 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  51.236647529872606 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -36.2504 * tfactors.T943i
                                  + (5.0/3.0) * -5.0 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_F19_to_p_O18_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  92.7618744782197 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -16.7246 * tfactors.T943i
                                  + (5.0/3.0) * -3.0 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  99.4871744782197 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 3.99059 * tfactors.T923i
                                  + -0.593127 + (5.0/3.0) * 0.0877534 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  94.41868447821969 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  93.0066174782197 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_F19_to_He4_N15_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  50.77806404408662 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 35.4292 * tfactors.T923i
                                  + -5.5767 + (5.0/3.0) * 0.441293 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  53.56282404408662 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  50.75615404408662 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  46.57820404408662 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -36.2324 * tfactors.T943i
                                  + (5.0/3.0) * -2.0 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Ne18_to_p_F17_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ne18 --> p + F17

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real F17_log_pf, dF17_log_pf_dT9;

    // interpolating F17 partition function
    get_partition_function_cached(F17, tfactors.T9, pf_cache, F17_log_pf, dF17_log_pf_dT9);

    amrex::Real Ne18_log_pf, dNe18_log_pf_dT9;

    // interpolating Ne18 partition function
    get_partition_function_cached(Ne18, tfactors.T9, pf_cache, Ne18_log_pf, dNe18_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + F17_log_pf - Ne18_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dF17_log_pf_dT9 - dNe18_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from cb09
    ln_set_rate =  17.58370046698957 + -45.55769965436449 * tfactors.T9i + -14.2191 * tfactors.T913i + 34.0647 * tfactors.T913
                         + -16.5698 * tfactors.T9 + 2.48116 * tfactors.T953 + -0.6337600000000001 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  45.55769965436449 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -14.2191 * tfactors.T943i + (1.0/3.0) * 34.0647 * tfactors.T923i
                                  + -16.5698 + (5.0/3.0) * 2.48116 * tfactors.T923 + -0.6337600000000001 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from cb09
    ln_set_rate =  53.00858046698957 + -50.485039254364494 * tfactors.T9i + -21.3249 * tfactors.T913i + -0.230774 * tfactors.T913
                         + 0.917931 * tfactors.T9 + -0.0440377 * tfactors.T953 + -5.86014 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  50.485039254364494 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -21.3249 * tfactors.T943i + (1.0/3.0) * -0.230774 * tfactors.T923i
                                  + 0.917931 + (5.0/3.0) * -0.0440377 * tfactors.T923 + -5.86014 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Ne18_to_He4_O14_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ne18 --> He4 + O14

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real O14_log_pf, dO14_log_pf_dT9;

    // interpolating O14 partition function
    get_partition_function_cached(O14, tfactors.T9, pf_cache, O14_log_pf, dO14_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ne18_log_pf, dNe18_log_pf_dT9;

    // interpolating Ne18 partition function
    get_partition_function_cached(Ne18, tfactors.T9, pf_cache, Ne18_log_pf, dNe18_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + O14_log_pf - Ne18_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dO14_log_pf_dT9 - dNe18_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from wh87n
    ln_set_rate =  51.160347955038674 + -59.35823432141371 * tfactors.T9i + -39.38 * tfactors.T913i + -0.0772187 * tfactors.T913
                         + -0.635361 * tfactors.T9 + 0.106236 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  59.35823432141371 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -39.38 * tfactors.T943i + (1.0/3.0) * -0.0772187 * tfactors.T923i
                                  + -0.635361 + (5.0/3.0) * 0.106236 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from wh87r
    ln_set_rate =  20.01796795503867 + -71.51723432141371 * tfactors.T9i
                         + 6.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  71.51723432141371 * tfactors.T9i * tfactors.T9i
                                  + 6.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from wh87r
    ln_set_rate =  28.243807955038672 + -81.96823432141372 * tfactors.T9i;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  81.96823432141372 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from wh87r
    ln_set_rate =  22.56327795503867 + -71.08823432141371 * tfactors.T9i;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  71.08823432141371 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Ne19_to_p_F18_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ne19 --> p + F18

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

    amrex::Real Ne19_log_pf, dNe19_log_pf_dT9;

    // interpolating Ne19 partition function
    get_partition_function_cached(Ne19, tfactors.T9, pf_cache, Ne19_log_pf, dNe19_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + F18_log_pf - Ne19_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dF18_log_pf_dT9 - dNe19_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from il10r
    ln_set_rate =  18.192220240787112 + -77.27667559674641 * tfactors.T9i + 13.1683 * tfactors.T913
                         + -1.92023 * tfactors.T9 + 0.16901 * tfactors.T953;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  77.27667559674641 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 13.1683 * tfactors.T923i
                                  + -1.92023 + (5.0/3.0) * 0.16901 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  -5.39950975921289 + -74.78415559674642 * tfactors.T9i + 22.4903 * tfactors.T913
                         + 0.307872 * tfactors.T9 + -0.296226 * tfactors.T953;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  74.78415559674642 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 22.4903 * tfactors.T923i
                                  + 0.307872 + (5.0/3.0) * -0.296226 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10n
    ln_set_rate =  81.45789024078711 + -74.38520559674642 * tfactors.T9i + -21.4023 * tfactors.T913i + -93.766 * tfactors.T913
                         + 179.258 * tfactors.T9 + -202.561 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  74.38520559674642 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -21.4023 * tfactors.T943i + (1.0/3.0) * -93.766 * tfactors.T923i
                                  + 179.258 + (5.0/3.0) * -202.561 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Ne19_to_He4_O15_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Ne19 --> He4 + O15

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real O15_log_pf, dO15_log_pf_dT9;

    // interpolating O15 partition function
    get_partition_function_cached(O15, tfactors.T9, pf_cache, O15_log_pf, dO15_log_pf_dT9);

    amrex::Real Ne19_log_pf, dNe19_log_pf_dT9;

    // interpolating Ne19 partition function
    get_partition_function_cached(Ne19, tfactors.T9, pf_cache, Ne19_log_pf, dNe19_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + O15_log_pf - Ne19_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dO15_log_pf_dT9 - dNe19_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from dc11r
    ln_set_rate =  -7.510242070092687 + -45.15053612970033 * tfactors.T9i + -3.24609 * tfactors.T913i + 44.4647 * tfactors.T913
                         + -9.79962 * tfactors.T9 + 0.841782 * tfactors.T953;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  45.15053612970033 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -3.24609 * tfactors.T943i + (1.0/3.0) * 44.4647 * tfactors.T923i
                                  + -9.79962 + (5.0/3.0) * 0.841782 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from dc11r
    ln_set_rate =  24.694111429907313 + -46.83053612970032 * tfactors.T9i;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  46.83053612970032 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from dc11n
    ln_set_rate =  51.03075792990731 + -40.946146129700324 * tfactors.T9i + -39.578 * tfactors.T913i
                         + -3.0 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  40.946146129700324 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -39.578 * tfactors.T943i
                                  + (5.0/3.0) * -3.0 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Ne20_to_p_F19_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  156.78628300750532 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 31.6442 * tfactors.T943i + (1.0/3.0) * -58.6563 * tfactors.T923i
                                  + 67.7365 + (5.0/3.0) * -22.9721 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  150.75597300750533 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -11.3832 * tfactors.T943i + (1.0/3.0) * 5.47872 * tfactors.T923i
                                  + -1.07203 + (5.0/3.0) * 0.11196 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  149.04214300750533 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -18.116 * tfactors.T943i + (1.0/3.0) * -1.4622 * tfactors.T923i
                                  + 6.95113 + (5.0/3.0) * -2.90366 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Ne20_to_He4_O16_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  65.24608327101186 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  54.887583271011856 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -39.7262 * tfactors.T943i + (1.0/3.0) * -0.210799 * tfactors.T923i
                                  + 0.442879 + (5.0/3.0) * -0.0797753 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  67.65188327101185 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -3.65925 * tfactors.T923i
                                  + 0.714224 + (5.0/3.0) * -0.00107508 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Ne21_to_He4_O17_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  90.7784662204353 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 123.363 * tfactors.T943i + (1.0/3.0) * -87.4351 * tfactors.T923i
                                  + -3.40974e-06 + (5.0/3.0) * -57.0469 * tfactors.T923 + 83.7218 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  98.9451062204353 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 3.31162e-08 * tfactors.T943i + (1.0/3.0) * 130.258 * tfactors.T923i
                                  + -7.92551e-05 + (5.0/3.0) * -4.13772 * tfactors.T923 + -41.2753 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  91.26872622043531 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 2.87641 * tfactors.T943i + (1.0/3.0) * -3.54489 * tfactors.T923i
                                  + -2.11222e-08 + (5.0/3.0) * -3.90649e-09 * tfactors.T923 + 6.25778 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Ne22_to_He4_O18_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  114.18179645973349 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  112.85484845973349 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  143.2255364597335 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  113.76407645973349 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -44.3823 * tfactors.T943i + (1.0/3.0) * -46.6617 * tfactors.T923i
                                  + 7.88059 + (5.0/3.0) * -0.590829 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Na21_to_p_Ne20_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Na21 --> p + Ne20

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Na21_log_pf, dNa21_log_pf_dT9;

    // interpolating Na21 partition function
    get_partition_function_cached(Na21, tfactors.T9, pf_cache, Na21_log_pf, dNa21_log_pf_dT9);

    amrex::Real Ne20_log_pf, dNe20_log_pf_dT9;

    // interpolating Ne20 partition function
    get_partition_function_cached(Ne20, tfactors.T9, pf_cache, Ne20_log_pf, dNe20_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Ne20_log_pf - Na21_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dNe20_log_pf_dT9 - dNa21_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ly18
    ln_set_rate =  195319.2650977437 + -89.3610144370849 * tfactors.T9i + 21894.7 * tfactors.T913i + -319153.0 * tfactors.T913
                         + 224369.0 * tfactors.T9 + -188049.0 * tfactors.T953 + 48704.9 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  89.3610144370849 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 21894.7 * tfactors.T943i + (1.0/3.0) * -319153.0 * tfactors.T923i
                                  + 224369.0 + (5.0/3.0) * -188049.0 * tfactors.T923 + 48704.9 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from ly18
    ln_set_rate =  230.14209774368004 + -28.373725437084907 * tfactors.T9i + 15.325 * tfactors.T913i + -294.859 * tfactors.T913
                         + 107.692 * tfactors.T9 + -46.2072 * tfactors.T953 + 59.3398 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  28.373725437084907 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 15.325 * tfactors.T943i + (1.0/3.0) * -294.859 * tfactors.T923i
                                  + 107.692 + (5.0/3.0) * -46.2072 * tfactors.T923 + 59.3398 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from ly18
    ln_set_rate =  28.096127743680043 + -37.059014437084905 * tfactors.T9i + 20.5893 * tfactors.T913i + -17.5841 * tfactors.T913
                         + 0.243226 * tfactors.T9 + -0.000231418 * tfactors.T953 + 14.3398 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  37.059014437084905 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 20.5893 * tfactors.T943i + (1.0/3.0) * -17.5841 * tfactors.T923i
                                  + 0.243226 + (5.0/3.0) * -0.000231418 * tfactors.T923 + 14.3398 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from ly18
    ln_set_rate =  252.28409774368004 + -32.674594437084906 * tfactors.T9i + 258.57 * tfactors.T913i + -506.387 * tfactors.T913
                         + 22.1576 * tfactors.T9 + -0.721182 * tfactors.T953 + 231.788 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  32.674594437084906 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 258.57 * tfactors.T943i + (1.0/3.0) * -506.387 * tfactors.T923i
                                  + 22.1576 + (5.0/3.0) * -0.721182 * tfactors.T923 + 231.788 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Na21_to_He4_F17_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Na21 --> He4 + F17

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real F17_log_pf, dF17_log_pf_dT9;

    // interpolating F17 partition function
    get_partition_function_cached(F17, tfactors.T9, pf_cache, F17_log_pf, dF17_log_pf_dT9);

    amrex::Real Na21_log_pf, dNa21_log_pf_dT9;

    // interpolating Na21 partition function
    get_partition_function_cached(Na21, tfactors.T9, pf_cache, Na21_log_pf, dNa21_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + F17_log_pf - Na21_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dF17_log_pf_dT9 - dNa21_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from rpsmr
    ln_set_rate =  66.33566650730282 + -77.87093563755207 * tfactors.T9i + 15.559 * tfactors.T913i + -68.3231 * tfactors.T913
                         + 2.54275 * tfactors.T9 + -0.0989207 * tfactors.T953 + 38.3877 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  77.87093563755207 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 15.559 * tfactors.T943i + (1.0/3.0) * -68.3231 * tfactors.T923i
                                  + 2.54275 + (5.0/3.0) * -0.0989207 * tfactors.T923 + 38.3877 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Na22_to_p_Ne21_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  78.19798607071328 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -19.2096 * tfactors.T943i
                                  + (5.0/3.0) * -1.0 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  82.41183607071328 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 21.1176 * tfactors.T943i + (1.0/3.0) * 34.0411 * tfactors.T923i
                                  + -4.45593 + (5.0/3.0) * 0.328613 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  79.59755607071328 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  78.39416607071328 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Na22_to_He4_F18_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  100.22898907841348 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 18.8956 * tfactors.T943i + (1.0/3.0) * -65.6134 * tfactors.T923i
                                  + 1.71114 + (5.0/3.0) * -0.0260999 * tfactors.T923 + 39.3396 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Na23_to_p_Ne22_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  103.77146219384292 * tfactors.T9i * tfactors.T9i
                                  + 1.189235 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  102.46235419384291 * tfactors.T9i * tfactors.T9i
                                  + 0.009810000000000096 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  102.85114819384292 * tfactors.T9i * tfactors.T9i
                                  + 4.735580000000001 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  110.77516219384292 * tfactors.T9i * tfactors.T9i
                                  + 0.732533 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  106.65075219384292 * tfactors.T9i * tfactors.T9i
                                  + 1.656226 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  104.66929219384292 * tfactors.T9i * tfactors.T9i
                                  + -2.79964 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Na23_to_He4_F19_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  123.58233417534949 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 39.7219 * tfactors.T943i + (1.0/3.0) * -100.401 * tfactors.T923i
                                  + 3.15808 + (5.0/3.0) * -0.0629822 * tfactors.T923 + 55.9823 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Mg24_to_p_Na23_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  138.96838756220697 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.360588 * tfactors.T923i
                                  + 1.4187 + (5.0/3.0) * -0.184061 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  137.30028756220696 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  135.68809756220696 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -20.6428 * tfactors.T943i + (1.0/3.0) * 1.52954 * tfactors.T923i
                                  + 2.7487 + (5.0/3.0) * -1.0 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Mg24_to_He4_Ne20_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  120.89510873006793 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 16.9229 * tfactors.T923i
                                  + -2.57325 + (5.0/3.0) * 0.208997 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  117.33446873006793 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  110.62025873006793 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  108.11420873006793 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -46.2525 * tfactors.T943i + (1.0/3.0) * 5.58901 * tfactors.T923i
                                  + 7.61843 + (5.0/3.0) * -3.683 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Al27_to_He4_Na23_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  117.11189602392653 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -50.2042 * tfactors.T943i + (1.0/3.0) * -1.64239 * tfactors.T923i
                                  + -1.59995 + (5.0/3.0) * 0.184933 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Si28_to_p_Al27_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  136.3412367256539 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 23.8634 * tfactors.T923i
                                  + -3.70135 + (5.0/3.0) * 0.28964 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  134.8245897256539 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -26.8327 * tfactors.T943i + (1.0/3.0) * -116.137 * tfactors.T923i
                                  + 0.00950567 + (5.0/3.0) * 0.00999755 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  134.4372767256539 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -23.2205 * tfactors.T943i
                                  + (5.0/3.0) * -2.0 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Si28_to_He4_Mg24_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  131.4900751873663 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  128.6942751873663 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 21.3721 * tfactors.T943i + (1.0/3.0) * 37.7649 * tfactors.T923i
                                  + -4.10635 + (5.0/3.0) * 0.249618 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_P31_to_He4_Al27_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  112.19944704325766 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -56.5351 * tfactors.T943i + (1.0/3.0) * -0.896208 * tfactors.T923i
                                  + -1.72024 + (5.0/3.0) * 0.185409 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_S32_to_p_P31_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  106.63906812813134 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 8.09341 * tfactors.T923i
                                  + -0.615971 + (5.0/3.0) * 0.031159 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  105.12160812813134 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  102.86202812813134 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -25.3278 * tfactors.T943i + (1.0/3.0) * 6.4931 * tfactors.T923i
                                  + -9.27513 + (5.0/3.0) * -0.610439 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_S32_to_He4_Si28_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  80.62419844573513 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -59.4896 * tfactors.T943i + (1.0/3.0) * 4.47205 * tfactors.T923i
                                  + -4.78989 + (5.0/3.0) * 0.557201 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_C12_to_He4_He4_He4_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  85.44440046993657 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -23.57 * tfactors.T943i + (1.0/3.0) * 20.4886 * tfactors.T923i
                                  + -12.9882 + (5.0/3.0) * -20.0 * tfactors.T923 + 0.8333300000000001 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  84.41994046993656 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -37.06 * tfactors.T943i + (1.0/3.0) * 29.3493 * tfactors.T923i
                                  + -115.507 + (5.0/3.0) * -10.0 * tfactors.T923 + 1.66667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  88.54650046993656 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -13.49 * tfactors.T943i + (1.0/3.0) * 21.4259 * tfactors.T923i
                                  + -1.34769 + (5.0/3.0) * 0.0879816 * tfactors.T923 + -10.1653 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_He4_C12_to_p_N15_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  57.62215691264642 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -15.253 * tfactors.T943i + (1.0/3.0) * 1.59318 * tfactors.T923i
                                  + 2.4479 + (5.0/3.0) * -2.19708 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  58.78595691264642 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 22.7105 * tfactors.T923i
                                  + -2.90707 + (5.0/3.0) * 0.205754 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  65.02815691264642 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  59.643326912646415 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 30.8497 * tfactors.T923i
                                  + -8.50433 + (5.0/3.0) * -1.54426 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_He4_N14_to_p_O17_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  14.584520682879308 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  13.831125682879309 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -16.9078 * tfactors.T943i
                                  + (5.0/3.0) * -2.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  18.78977568287931 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 5.10182 * tfactors.T923i
                                  + 0.379373 + (5.0/3.0) * -0.0672515 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  15.945895682879309 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_He4_N15_to_p_O18_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  46.42955443414268 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  46.183670434142684 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -16.6979 * tfactors.T943i
                                  + (5.0/3.0) * -3.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  51.50702043414268 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 11.6568 * tfactors.T923i
                                  + -2.16303 + (5.0/3.0) * 0.209965 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  47.84667043414268 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_He4_O15_to_p_F18_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // O15 + He4 --> p + F18

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real O15_log_pf, dO15_log_pf_dT9;

    // interpolating O15 partition function
    get_partition_function_cached(O15, tfactors.T9, pf_cache, O15_log_pf, dO15_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real F18_log_pf, dF18_log_pf_dT9;

    // interpolating F18 partition function
    get_partition_function_cached(F18, tfactors.T9, pf_cache, F18_log_pf, dF18_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = p_log_pf + F18_log_pf - He4_log_pf - O15_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dF18_log_pf_dT9 - dHe4_log_pf_dT9 - dO15_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from il10r
    ln_set_rate =  1.0671723108797964 + -36.45580946704611 * tfactors.T9i + 13.3223 * tfactors.T913
                         + -1.36696 * tfactors.T9 + 0.0757363 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  36.45580946704611 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 13.3223 * tfactors.T923i
                                  + -1.36696 + (5.0/3.0) * 0.0757363 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10r
    ln_set_rate =  -32.4286676891202 + -33.81549146704611 * tfactors.T9i + 61.738 * tfactors.T913
                         + -108.29 * tfactors.T9 + -34.2365 * tfactors.T953 + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  33.81549146704611 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 61.738 * tfactors.T923i
                                  + -108.29 + (5.0/3.0) * -34.2365 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from il10n
    ln_set_rate =  61.315932310879795 + -33.43905946704611 * tfactors.T9i + -21.4023 * tfactors.T913i + -80.8891 * tfactors.T913
                         + 134.6 * tfactors.T9 + -126.504 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  33.43905946704611 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -21.4023 * tfactors.T943i + (1.0/3.0) * -80.8891 * tfactors.T923i
                                  + 134.6 + (5.0/3.0) * -126.504 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_p_O16_to_He4_N13_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  60.557329596678294 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -35.829 * tfactors.T943i + (1.0/3.0) * -0.530275 * tfactors.T923i
                                  + -0.982462 + (5.0/3.0) * 0.0808059 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_He4_O16_to_p_F19_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  94.28220973651266 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  94.15455973651265 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -18.116 * tfactors.T943i
                                  + 1.86674 + (5.0/3.0) * -7.5666 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  97.44055973651265 * tfactors.T9i * tfactors.T9i
                                  + -0.21103 + 2.87702 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  97.90640973651266 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  96.62283973651265 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_p_F17_to_He4_O14_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // F17 + p --> He4 + O14

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real O14_log_pf, dO14_log_pf_dT9;

    // interpolating O14 partition function
    get_partition_function_cached(O14, tfactors.T9, pf_cache, O14_log_pf, dO14_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real F17_log_pf, dF17_log_pf_dT9;

    // interpolating F17 partition function
    get_partition_function_cached(F17, tfactors.T9, pf_cache, F17_log_pf, dF17_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + O14_log_pf - p_log_pf - F17_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dO14_log_pf_dT9 - dp_log_pf_dT9 - dF17_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from Ha96r
    ln_set_rate =  11.415567488049104 + -25.85518506704921 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  25.85518506704921 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from Ha96r
    ln_set_rate =  17.938467488049106 + -39.83288506704921 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  39.83288506704921 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from Ha96n
    ln_set_rate =  40.122467488049104 + -13.83288506704921 * tfactors.T9i + -39.388 * tfactors.T913i + -17.4673 * tfactors.T913
                         + 35.3029 * tfactors.T9 + -24.8162 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  13.83288506704921 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -39.388 * tfactors.T943i + (1.0/3.0) * -17.4673 * tfactors.T923i
                                  + 35.3029 + (5.0/3.0) * -24.8162 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from Ha96r
    ln_set_rate =  15.595367488049106 + -36.34288506704921 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  36.34288506704921 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from Ha96r
    ln_set_rate =  10.405067488049104 + -27.43288506704921 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  27.43288506704921 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from Ha96r
    ln_set_rate =  -106.80433251195089 + -14.285921067049209 * tfactors.T9i
                         + -1.5 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  14.285921067049209 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_He4_F17_to_p_Ne20_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // F17 + He4 --> p + Ne20

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real F17_log_pf, dF17_log_pf_dT9;

    // interpolating F17 partition function
    get_partition_function_cached(F17, tfactors.T9, pf_cache, F17_log_pf, dF17_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ne20_log_pf, dNe20_log_pf_dT9;

    // interpolating Ne20 partition function
    get_partition_function_cached(Ne20, tfactors.T9, pf_cache, Ne20_log_pf, dNe20_log_pf_dT9);

    amrex::Real net_log_pf = p_log_pf + Ne20_log_pf - He4_log_pf - F17_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dNe20_log_pf_dT9 - dHe4_log_pf_dT9 - dF17_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from nacr
    ln_set_rate =  38.645331236377224 + -0.004848799532844339 * tfactors.T9i + -43.18 * tfactors.T913i + 4.46827 * tfactors.T913
                         + -1.63915 * tfactors.T9 + 0.123483 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  0.004848799532844339 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -43.18 * tfactors.T943i + (1.0/3.0) * 4.46827 * tfactors.T923i
                                  + -1.63915 + (5.0/3.0) * 0.123483 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_He4_Ne20_to_p_Na23_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  29.889658832139027 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 19.7297 * tfactors.T923i
                                  + -2.20987 + (5.0/3.0) * 0.153374 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  29.434918832139026 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  27.573888832139026 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -20.0024 * tfactors.T943i + (1.0/3.0) * 11.5988 * tfactors.T923i
                                  + -1.37398 + (5.0/3.0) * -1.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_p_Ne21_to_He4_F18_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  22.04702300768339 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 21.4461 * tfactors.T943i + (1.0/3.0) * -73.252 * tfactors.T923i
                                  + 2.42329 + (5.0/3.0) * -0.077278 * tfactors.T923 + 40.7604 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_p_Ne22_to_He4_F19_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  22.26126198150658 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -38.7722 * tfactors.T943i + (1.0/3.0) * -13.3654 * tfactors.T923i
                                  + 0.863648 + (5.0/3.0) * -0.0451491 * tfactors.T923 + 1.33333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  65.19766198150657 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -34.5008 * tfactors.T943i + (1.0/3.0) * 56.9316 * tfactors.T923i
                                  + 2.09613 + (5.0/3.0) * -32.496 * tfactors.T923 + 0.333333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  152.44286198150658 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 12625.1 * tfactors.T943i + (1.0/3.0) * -49107.1 * tfactors.T923i
                                  + 9227.53 + (5.0/3.0) * -2086.65 * tfactors.T923 + 14520.2 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_p_Na21_to_He4_Ne18_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Na21 + p --> He4 + Ne18

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Na21_log_pf, dNa21_log_pf_dT9;

    // interpolating Na21 partition function
    get_partition_function_cached(Na21, tfactors.T9, pf_cache, Na21_log_pf, dNa21_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ne18_log_pf, dNe18_log_pf_dT9;

    // interpolating Ne18 partition function
    get_partition_function_cached(Ne18, tfactors.T9, pf_cache, Ne18_log_pf, dNe18_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Ne18_log_pf - p_log_pf - Na21_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dNe18_log_pf_dT9 - dp_log_pf_dT9 - dNa21_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from mo14r
    ln_set_rate =  -23.82911395968675 + -36.38615638318757 * tfactors.T9i + 8.5632e-08 * tfactors.T913i + 24.8579 * tfactors.T913
                         + 0.0823845 * tfactors.T9 + -0.365374 * tfactors.T953 + -2.21415e-06 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  36.38615638318757 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 8.5632e-08 * tfactors.T943i + (1.0/3.0) * 24.8579 * tfactors.T923i
                                  + 0.0823845 + (5.0/3.0) * -0.365374 * tfactors.T923 + -2.21415e-06 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from mo14r
    ln_set_rate =  -22.686913959686752 + -40.87721638318757 * tfactors.T9i + 4.73034e-07 * tfactors.T913i + 36.29 * tfactors.T913
                         + -6.56565 * tfactors.T9 + -2.96287e-06 * tfactors.T953 + -9.00373e-07 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  40.87721638318757 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 4.73034e-07 * tfactors.T943i + (1.0/3.0) * 36.29 * tfactors.T923i
                                  + -6.56565 + (5.0/3.0) * -2.96287e-06 * tfactors.T923 + -9.00373e-07 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
    }

    // avoid underflows by zeroing rates in [0.0, 1.e-100]
    ln_set_rate = std::max(ln_set_rate, -230.0);
    set_rate = std::exp(ln_set_rate);
    rate += set_rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
    }

    // ReacLib set derived from mo14r
    ln_set_rate =  -37.64591395968675 + -34.48551638318757 * tfactors.T9i + 7.09521e-06 * tfactors.T913i + 34.1789 * tfactors.T913
                         + -1.72974 * tfactors.T9 + -0.0395081 * tfactors.T953 + -7.82759e-07 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  34.48551638318757 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 7.09521e-06 * tfactors.T943i + (1.0/3.0) * 34.1789 * tfactors.T923i
                                  + -1.72974 + (5.0/3.0) * -0.0395081 * tfactors.T923 + -7.82759e-07 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_He4_Na21_to_p_Mg24_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Na21 + He4 --> p + Mg24

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Na21_log_pf, dNa21_log_pf_dT9;

    // interpolating Na21 partition function
    get_partition_function_cached(Na21, tfactors.T9, pf_cache, Na21_log_pf, dNa21_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Mg24_log_pf, dMg24_log_pf_dT9;

    // interpolating Mg24 partition function
    get_partition_function_cached(Mg24, tfactors.T9, pf_cache, Mg24_log_pf, dMg24_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real net_log_pf = p_log_pf + Mg24_log_pf - He4_log_pf - Na21_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dMg24_log_pf_dT9 - dHe4_log_pf_dT9 - dNa21_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from nacr
    ln_set_rate =  39.83098890373066 + -0.00380570701699412 * tfactors.T9i + -49.9621 * tfactors.T913i + 5.90498 * tfactors.T913
                         + -1.6598 * tfactors.T9 + 0.117817 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  0.00380570701699412 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -49.9621 * tfactors.T943i + (1.0/3.0) * 5.90498 * tfactors.T923i
                                  + -1.6598 + (5.0/3.0) * 0.117817 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_p_Na22_to_He4_Ne19_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Na22 + p --> He4 + Ne19

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real Na22_log_pf, dNa22_log_pf_dT9;

    // interpolating Na22 partition function
    get_partition_function_cached(Na22, tfactors.T9, pf_cache, Na22_log_pf, dNa22_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ne19_log_pf, dNe19_log_pf_dT9;

    // interpolating Ne19 partition function
    get_partition_function_cached(Ne19, tfactors.T9, pf_cache, Ne19_log_pf, dNe19_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Ne19_log_pf - p_log_pf - Na22_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dNe19_log_pf_dT9 - dp_log_pf_dT9 - dNa22_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  43.08404589973112 + -24.01421348166705 * tfactors.T9i + -46.6346 * tfactors.T913i + 0.866532 * tfactors.T913
                         + -0.893541 * tfactors.T9 + 0.0747971 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  24.01421348166705 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -46.6346 * tfactors.T943i + (1.0/3.0) * 0.866532 * tfactors.T923i
                                  + -0.893541 + (5.0/3.0) * 0.0747971 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_He4_Mg24_to_p_Al27_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  22.818701538287602 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 18.0416 * tfactors.T923i
                                  + -1.54137 + (5.0/3.0) * 0.0847506 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  19.5392135382876 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 5.18642 * tfactors.T943i + (1.0/3.0) * -34.7936 * tfactors.T923i
                                  + 168.225 + (5.0/3.0) * -115.825 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  18.576201538287602 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -26.4162 * tfactors.T943i
                                  + (5.0/3.0) * -2.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_He4_Si28_to_p_P31_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  25.663579682396225 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 21.521 * tfactors.T923i
                                  + -1.90355 + (5.0/3.0) * 0.092724 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  24.114989682396224 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  22.237829682396224 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -31.932 * tfactors.T943i + (1.0/3.0) * -77.0334 * tfactors.T923i
                                  + -43.6847 + (5.0/3.0) * -4.28955 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_He4_S32_to_p_Cl35_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  25.589703706656966 * tfactors.T9i * tfactors.T9i
                                  + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  22.187684706656967 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 25.5338 * tfactors.T923i
                                  + 6.45824 + (5.0/3.0) * -0.950294 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  21.654753706656965 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -30.9147 * tfactors.T943i + (1.0/3.0) * -1.2345 * tfactors.T923i
                                  + 22.5118 + (5.0/3.0) * -33.0589 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  27.664513706656965 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 5.33756 * tfactors.T923i
                                  + 1.64418 + (5.0/3.0) * -0.246167 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Ar36_to_He4_S32_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  77.06468541928012 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -65.3709 * tfactors.T943i + (1.0/3.0) * 5.68294 * tfactors.T923i
                                  + -5.00388 + (5.0/3.0) * 0.571407 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Ar36_to_p_Cl35_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  99.28409012589869 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  98.71943912589869 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -27.8971 * tfactors.T943i + (1.0/3.0) * -16.2304 * tfactors.T923i
                                  + 35.255 + (5.0/3.0) * -25.8411 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  102.37035912589869 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 18.0179 * tfactors.T923i
                                  + -2.86304 + (5.0/3.0) * 0.250854 * tfactors.T923;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  100.7293991258987 * tfactors.T9i * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_He4_Ar36_to_p_K39_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  14.951286608188834 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -30.0732 * tfactors.T943i + (1.0/3.0) * 7.03263 * tfactors.T923i
                                  + -1.10085 + (5.0/3.0) * 0.133768 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Ca40_to_He4_Ar36_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  81.69320665503751 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -71.0046 * tfactors.T943i + (1.0/3.0) * 4.0656 * tfactors.T923i
                                  + -5.26509 + (5.0/3.0) * 0.683546 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Ca40_to_p_K39_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  101.86683326322635 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 802.18 * tfactors.T943i + (1.0/3.0) * -4010.27 * tfactors.T923i
                                  + 1137.69 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  109.20919326322635 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 641.844 * tfactors.T943i + (1.0/3.0) * -1248.49 * tfactors.T923i
                                  + 566.426 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  98.30957326322635 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 41.1723 * tfactors.T943i + (1.0/3.0) * -149.299 * tfactors.T923i
                                  + 10.5229 + (5.0/3.0) * -0.68208 * tfactors.T923 + 60.7367 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_He4_Ca40_to_p_Sc43_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  40.87525788838128 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -32.1734 * tfactors.T943i + (1.0/3.0) * 0.0296879 * tfactors.T923i
                                  + -0.95232 + (5.0/3.0) * 0.129022 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Ti44_to_He4_Ca40_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  59.49768581557457 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -76.4273 * tfactors.T943i + (1.0/3.0) * 3.87451 * tfactors.T923i
                                  + -3.61477 + (5.0/3.0) * 0.367451 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Ti44_to_p_Sc43_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  100.37294370395585 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -32.1734 * tfactors.T943i + (1.0/3.0) * -1.77078 * tfactors.T923i
                                  + -2.21706 + (5.0/3.0) * 0.298499 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Cr48_to_He4_Ti44_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  89.32364150067336 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -81.667 * tfactors.T943i + (1.0/3.0) * -10.6333 * tfactors.T923i
                                  + -0.672613 + (5.0/3.0) * 0.161209 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Cr48_to_p_V47_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  100.08956223422913 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -34.0548 * tfactors.T943i + (1.0/3.0) * -3.41973 * tfactors.T923i
                                  + 1.16501 + (5.0/3.0) * -0.105543 * tfactors.T923 + -6.20886 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  99.32517223422914 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 317.171 * tfactors.T943i + (1.0/3.0) * -911.679 * tfactors.T923i
                                  + 94.4245 + (5.0/3.0) * -10.1973 * tfactors.T923 + 332.227 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  93.78359723422915 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -45.9868 * tfactors.T943i + (1.0/3.0) * 13.6822 * tfactors.T923i
                                  + -0.376902 + (5.0/3.0) * -0.0194875 * tfactors.T923 + -6.9232499999999995 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  94.54467623422914 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -110.655 * tfactors.T943i + (1.0/3.0) * 83.0232 * tfactors.T923i
                                  + -19.7762 + (5.0/3.0) * 3.03961 * tfactors.T923 + -47.9742 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_p_V47_to_He4_Ti44_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  6.086479266444236 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 70.2835 * tfactors.T923i
                                  + -7.99061 + (5.0/3.0) * 0.486213 * tfactors.T923 + -1.5 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Fe52_to_He4_Cr48_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  92.09363997638769 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -86.7459 * tfactors.T943i + (1.0/3.0) * -9.79373 * tfactors.T923i
                                  + -0.772169 + (5.0/3.0) * 0.155883 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Fe52_to_p_Mn51_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  85.61663846070292 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -36.1825 * tfactors.T943i + (1.0/3.0) * 0.873042 * tfactors.T923i
                                  + -2.89731 + (5.0/3.0) * 0.364394 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_p_Mn51_to_He4_Cr48_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  6.47700151570639 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -86.7459 * tfactors.T943i + (1.0/3.0) * 1.05653 * tfactors.T923i
                                  + -1.15757 + (5.0/3.0) * 0.0877546 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Ni56_to_He4_Fe52_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  92.84097094326391 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -91.6819 * tfactors.T943i + (1.0/3.0) * -9.51885 * tfactors.T923i
                                  + -0.533014 + (5.0/3.0) * 0.0892607 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_Ni56_to_p_Co55_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  83.16460378149378 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -38.1053 * tfactors.T943i + (1.0/3.0) * -0.210947 * tfactors.T923i
                                  + -2.68377 + (5.0/3.0) * 0.355814 * tfactors.T923 + 0.833333 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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
void rate_p_Co55_to_He4_Fe52_removed(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  9.676367161770125 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -91.6819 * tfactors.T943i + (1.0/3.0) * -0.329235 * tfactors.T923i
                                  + -0.780924 + (5.0/3.0) * 0.0425179 * tfactors.T923 + -0.666667 * tfactors.T9i;

        dln_set_rate_dT9 += net_dlog_pf_dT9;
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

    rate_N13_to_p_C12_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_N13_to_p_C12_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_N13_to_p_C12_derived) = drate_dT;

    }
    rate_N14_to_p_C13_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_N14_to_p_C13_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_N14_to_p_C13_derived) = drate_dT;

    }
    rate_O14_to_p_N13_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_O14_to_p_N13_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_O14_to_p_N13_derived) = drate_dT;

    }
    rate_O15_to_p_N14_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_O15_to_p_N14_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_O15_to_p_N14_derived) = drate_dT;

    }
    rate_O16_to_p_N15_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_O16_to_p_N15_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_O16_to_p_N15_derived) = drate_dT;

    }
    rate_O16_to_He4_C12_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_O16_to_He4_C12_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_O16_to_He4_C12_derived) = drate_dT;

    }
    rate_F17_to_p_O16_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_F17_to_p_O16_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_F17_to_p_O16_derived) = drate_dT;

    }
    rate_F18_to_p_O17_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_F18_to_p_O17_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_F18_to_p_O17_derived) = drate_dT;

    }
    rate_F18_to_He4_N14_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_F18_to_He4_N14_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_F18_to_He4_N14_derived) = drate_dT;

    }
    rate_F19_to_p_O18_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_F19_to_p_O18_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_F19_to_p_O18_derived) = drate_dT;

    }
    rate_F19_to_He4_N15_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_F19_to_He4_N15_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_F19_to_He4_N15_derived) = drate_dT;

    }
    rate_Ne18_to_p_F17_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Ne18_to_p_F17_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ne18_to_p_F17_derived) = drate_dT;

    }
    rate_Ne18_to_He4_O14_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Ne18_to_He4_O14_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ne18_to_He4_O14_derived) = drate_dT;

    }
    rate_Ne19_to_p_F18_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Ne19_to_p_F18_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ne19_to_p_F18_derived) = drate_dT;

    }
    rate_Ne19_to_He4_O15_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Ne19_to_He4_O15_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ne19_to_He4_O15_derived) = drate_dT;

    }
    rate_Ne20_to_p_F19_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Ne20_to_p_F19_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ne20_to_p_F19_derived) = drate_dT;

    }
    rate_Ne20_to_He4_O16_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Ne20_to_He4_O16_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ne20_to_He4_O16_derived) = drate_dT;

    }
    rate_Ne21_to_He4_O17_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Ne21_to_He4_O17_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ne21_to_He4_O17_derived) = drate_dT;

    }
    rate_Ne22_to_He4_O18_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Ne22_to_He4_O18_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ne22_to_He4_O18_derived) = drate_dT;

    }
    rate_Na21_to_p_Ne20_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Na21_to_p_Ne20_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Na21_to_p_Ne20_derived) = drate_dT;

    }
    rate_Na21_to_He4_F17_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Na21_to_He4_F17_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Na21_to_He4_F17_derived) = drate_dT;

    }
    rate_Na22_to_p_Ne21_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Na22_to_p_Ne21_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Na22_to_p_Ne21_derived) = drate_dT;

    }
    rate_Na22_to_He4_F18_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Na22_to_He4_F18_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Na22_to_He4_F18_derived) = drate_dT;

    }
    rate_Na23_to_p_Ne22_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Na23_to_p_Ne22_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Na23_to_p_Ne22_derived) = drate_dT;

    }
    rate_Na23_to_He4_F19_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Na23_to_He4_F19_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Na23_to_He4_F19_derived) = drate_dT;

    }
    rate_Mg24_to_p_Na23_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Mg24_to_p_Na23_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Mg24_to_p_Na23_derived) = drate_dT;

    }
    rate_Mg24_to_He4_Ne20_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Mg24_to_He4_Ne20_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Mg24_to_He4_Ne20_derived) = drate_dT;

    }
    rate_Al27_to_He4_Na23_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Al27_to_He4_Na23_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Al27_to_He4_Na23_derived) = drate_dT;

    }
    rate_Si28_to_p_Al27_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Si28_to_p_Al27_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Si28_to_p_Al27_derived) = drate_dT;

    }
    rate_Si28_to_He4_Mg24_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Si28_to_He4_Mg24_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Si28_to_He4_Mg24_derived) = drate_dT;

    }
    rate_P31_to_He4_Al27_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_P31_to_He4_Al27_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_P31_to_He4_Al27_derived) = drate_dT;

    }
    rate_S32_to_p_P31_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_S32_to_p_P31_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_S32_to_p_P31_derived) = drate_dT;

    }
    rate_S32_to_He4_Si28_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_S32_to_He4_Si28_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_S32_to_He4_Si28_derived) = drate_dT;

    }
    rate_C12_to_He4_He4_He4_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_C12_to_He4_He4_He4_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_C12_to_He4_He4_He4_derived) = drate_dT;

    }
    rate_He4_C12_to_p_N15_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_He4_C12_to_p_N15_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_C12_to_p_N15_derived) = drate_dT;

    }
    rate_He4_N14_to_p_O17_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_He4_N14_to_p_O17_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_N14_to_p_O17_derived) = drate_dT;

    }
    rate_He4_N15_to_p_O18_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_He4_N15_to_p_O18_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_N15_to_p_O18_derived) = drate_dT;

    }
    rate_He4_O15_to_p_F18_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_He4_O15_to_p_F18_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_O15_to_p_F18_derived) = drate_dT;

    }
    rate_p_O16_to_He4_N13_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_p_O16_to_He4_N13_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_O16_to_He4_N13_derived) = drate_dT;

    }
    rate_He4_O16_to_p_F19_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_He4_O16_to_p_F19_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_O16_to_p_F19_derived) = drate_dT;

    }
    rate_p_F17_to_He4_O14_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_p_F17_to_He4_O14_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_F17_to_He4_O14_derived) = drate_dT;

    }
    rate_He4_F17_to_p_Ne20_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_He4_F17_to_p_Ne20_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_F17_to_p_Ne20_derived) = drate_dT;

    }
    rate_He4_Ne20_to_p_Na23_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_He4_Ne20_to_p_Na23_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Ne20_to_p_Na23_derived) = drate_dT;

    }
    rate_p_Ne21_to_He4_F18_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_p_Ne21_to_He4_F18_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Ne21_to_He4_F18_derived) = drate_dT;

    }
    rate_p_Ne22_to_He4_F19_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_p_Ne22_to_He4_F19_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Ne22_to_He4_F19_derived) = drate_dT;

    }
    rate_p_Na21_to_He4_Ne18_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_p_Na21_to_He4_Ne18_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Na21_to_He4_Ne18_derived) = drate_dT;

    }
    rate_He4_Na21_to_p_Mg24_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_He4_Na21_to_p_Mg24_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Na21_to_p_Mg24_derived) = drate_dT;

    }
    rate_p_Na22_to_He4_Ne19_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_p_Na22_to_He4_Ne19_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Na22_to_He4_Ne19_derived) = drate_dT;

    }
    rate_He4_Mg24_to_p_Al27_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_He4_Mg24_to_p_Al27_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Mg24_to_p_Al27_derived) = drate_dT;

    }
    rate_He4_Si28_to_p_P31_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_He4_Si28_to_p_P31_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Si28_to_p_P31_derived) = drate_dT;

    }
    rate_He4_S32_to_p_Cl35_removed<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_He4_S32_to_p_Cl35_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_S32_to_p_Cl35_removed) = drate_dT;

    }
    rate_Ar36_to_He4_S32_removed<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Ar36_to_He4_S32_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ar36_to_He4_S32_removed) = drate_dT;

    }
    rate_Ar36_to_p_Cl35_removed<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Ar36_to_p_Cl35_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ar36_to_p_Cl35_removed) = drate_dT;

    }
    rate_He4_Ar36_to_p_K39_removed<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_He4_Ar36_to_p_K39_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Ar36_to_p_K39_removed) = drate_dT;

    }
    rate_Ca40_to_He4_Ar36_removed<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Ca40_to_He4_Ar36_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ca40_to_He4_Ar36_removed) = drate_dT;

    }
    rate_Ca40_to_p_K39_removed<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Ca40_to_p_K39_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ca40_to_p_K39_removed) = drate_dT;

    }
    rate_He4_Ca40_to_p_Sc43_removed<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_He4_Ca40_to_p_Sc43_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Ca40_to_p_Sc43_removed) = drate_dT;

    }
    rate_Ti44_to_He4_Ca40_removed<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Ti44_to_He4_Ca40_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ti44_to_He4_Ca40_removed) = drate_dT;

    }
    rate_Ti44_to_p_Sc43_removed<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Ti44_to_p_Sc43_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ti44_to_p_Sc43_removed) = drate_dT;

    }
    rate_Cr48_to_He4_Ti44_removed<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Cr48_to_He4_Ti44_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Cr48_to_He4_Ti44_removed) = drate_dT;

    }
    rate_Cr48_to_p_V47_removed<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Cr48_to_p_V47_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Cr48_to_p_V47_removed) = drate_dT;

    }
    rate_p_V47_to_He4_Ti44_removed<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_p_V47_to_He4_Ti44_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_V47_to_He4_Ti44_removed) = drate_dT;

    }
    rate_Fe52_to_He4_Cr48_removed<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Fe52_to_He4_Cr48_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe52_to_He4_Cr48_removed) = drate_dT;

    }
    rate_Fe52_to_p_Mn51_removed<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Fe52_to_p_Mn51_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe52_to_p_Mn51_removed) = drate_dT;

    }
    rate_p_Mn51_to_He4_Cr48_removed<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_p_Mn51_to_He4_Cr48_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Mn51_to_He4_Cr48_removed) = drate_dT;

    }
    rate_Ni56_to_He4_Fe52_removed<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Ni56_to_He4_Fe52_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni56_to_He4_Fe52_removed) = drate_dT;

    }
    rate_Ni56_to_p_Co55_removed<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Ni56_to_p_Co55_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni56_to_p_Co55_removed) = drate_dT;

    }
    rate_p_Co55_to_He4_Fe52_removed<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_p_Co55_to_He4_Fe52_removed) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_p_Co55_to_He4_Fe52_removed) = drate_dT;

    }

}
#endif
