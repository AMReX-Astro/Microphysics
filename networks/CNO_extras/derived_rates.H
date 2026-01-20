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
    get_partition_function_cached(O14, tfactors, pf_cache, O14_log_pf, dO14_log_pf_dT9);

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
    get_partition_function_cached(O15, tfactors, pf_cache, O15_log_pf, dO15_log_pf_dT9);

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
    get_partition_function_cached(O16, tfactors, pf_cache, O16_log_pf, dO16_log_pf_dT9);

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
    get_partition_function_cached(O16, tfactors, pf_cache, O16_log_pf, dO16_log_pf_dT9);

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
    get_partition_function_cached(F17, tfactors, pf_cache, F17_log_pf, dF17_log_pf_dT9);

    amrex::Real O16_log_pf, dO16_log_pf_dT9;

    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors, pf_cache, O16_log_pf, dO16_log_pf_dT9);

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
    get_partition_function_cached(F18, tfactors, pf_cache, F18_log_pf, dF18_log_pf_dT9);

    amrex::Real O17_log_pf, dO17_log_pf_dT9;

    // interpolating O17 partition function
    get_partition_function_cached(O17, tfactors, pf_cache, O17_log_pf, dO17_log_pf_dT9);

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
    get_partition_function_cached(F18, tfactors, pf_cache, F18_log_pf, dF18_log_pf_dT9);

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
    get_partition_function_cached(O18, tfactors, pf_cache, O18_log_pf, dO18_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real F19_log_pf, dF19_log_pf_dT9;

    // interpolating F19 partition function
    get_partition_function_cached(F19, tfactors, pf_cache, F19_log_pf, dF19_log_pf_dT9);

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
    get_partition_function_cached(F19, tfactors, pf_cache, F19_log_pf, dF19_log_pf_dT9);

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
    get_partition_function_cached(F17, tfactors, pf_cache, F17_log_pf, dF17_log_pf_dT9);

    amrex::Real Ne18_log_pf, dNe18_log_pf_dT9;

    // interpolating Ne18 partition function
    get_partition_function_cached(Ne18, tfactors, pf_cache, Ne18_log_pf, dNe18_log_pf_dT9);

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
    get_partition_function_cached(O14, tfactors, pf_cache, O14_log_pf, dO14_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ne18_log_pf, dNe18_log_pf_dT9;

    // interpolating Ne18 partition function
    get_partition_function_cached(Ne18, tfactors, pf_cache, Ne18_log_pf, dNe18_log_pf_dT9);

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
    get_partition_function_cached(F18, tfactors, pf_cache, F18_log_pf, dF18_log_pf_dT9);

    amrex::Real Ne19_log_pf, dNe19_log_pf_dT9;

    // interpolating Ne19 partition function
    get_partition_function_cached(Ne19, tfactors, pf_cache, Ne19_log_pf, dNe19_log_pf_dT9);

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
    get_partition_function_cached(O15, tfactors, pf_cache, O15_log_pf, dO15_log_pf_dT9);

    amrex::Real Ne19_log_pf, dNe19_log_pf_dT9;

    // interpolating Ne19 partition function
    get_partition_function_cached(Ne19, tfactors, pf_cache, Ne19_log_pf, dNe19_log_pf_dT9);

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
    get_partition_function_cached(F19, tfactors, pf_cache, F19_log_pf, dF19_log_pf_dT9);

    amrex::Real Ne20_log_pf, dNe20_log_pf_dT9;

    // interpolating Ne20 partition function
    get_partition_function_cached(Ne20, tfactors, pf_cache, Ne20_log_pf, dNe20_log_pf_dT9);

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
    get_partition_function_cached(O16, tfactors, pf_cache, O16_log_pf, dO16_log_pf_dT9);

    amrex::Real Ne20_log_pf, dNe20_log_pf_dT9;

    // interpolating Ne20 partition function
    get_partition_function_cached(Ne20, tfactors, pf_cache, Ne20_log_pf, dNe20_log_pf_dT9);

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
void rate_Mg22_to_He4_Ne18_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

    // Mg22 --> He4 + Ne18

    rate = 0.0;
    drate_dT = 0.0;

    // Evaluate partition function terms

    amrex::Real Mg22_log_pf, dMg22_log_pf_dT9;

    // interpolating Mg22 partition function
    get_partition_function_cached(Mg22, tfactors, pf_cache, Mg22_log_pf, dMg22_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ne18_log_pf, dNe18_log_pf_dT9;

    // interpolating Ne18 partition function
    get_partition_function_cached(Ne18, tfactors, pf_cache, Ne18_log_pf, dNe18_log_pf_dT9);

    amrex::Real net_log_pf = He4_log_pf + Ne18_log_pf - Mg22_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dNe18_log_pf_dT9 - dMg22_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from ths8r
    ln_set_rate =  57.680013188545345 + -94.4898569232192 * tfactors.T9i + -46.4859 * tfactors.T913i + 0.956741 * tfactors.T913
                         + -0.914402 * tfactors.T9 + 0.0722478 * tfactors.T953 + 0.833333 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  94.4898569232192 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -46.4859 * tfactors.T943i + (1.0/3.0) * 0.956741 * tfactors.T923i
                                  + -0.914402 + (5.0/3.0) * 0.0722478 * tfactors.T923 + 0.833333 * tfactors.T9i;

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
    get_partition_function_cached(Mg24, tfactors, pf_cache, Mg24_log_pf, dMg24_log_pf_dT9);

    amrex::Real Ne20_log_pf, dNe20_log_pf_dT9;

    // interpolating Ne20 partition function
    get_partition_function_cached(Ne20, tfactors, pf_cache, Ne20_log_pf, dNe20_log_pf_dT9);

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
    get_partition_function_cached(O17, tfactors, pf_cache, O17_log_pf, dO17_log_pf_dT9);

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
    get_partition_function_cached(O18, tfactors, pf_cache, O18_log_pf, dO18_log_pf_dT9);

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
    get_partition_function_cached(O15, tfactors, pf_cache, O15_log_pf, dO15_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real F18_log_pf, dF18_log_pf_dT9;

    // interpolating F18 partition function
    get_partition_function_cached(F18, tfactors, pf_cache, F18_log_pf, dF18_log_pf_dT9);

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
    get_partition_function_cached(O16, tfactors, pf_cache, O16_log_pf, dO16_log_pf_dT9);

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
    get_partition_function_cached(O16, tfactors, pf_cache, O16_log_pf, dO16_log_pf_dT9);

    amrex::Real p_log_pf, dp_log_pf_dT9;

    // setting p log(partition function) to 0.0 by default, independent of T
    p_log_pf = 0.0_rt;
    dp_log_pf_dT9 = 0.0_rt;

    amrex::Real F19_log_pf, dF19_log_pf_dT9;

    // interpolating F19 partition function
    get_partition_function_cached(F19, tfactors, pf_cache, F19_log_pf, dF19_log_pf_dT9);

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
    get_partition_function_cached(O14, tfactors, pf_cache, O14_log_pf, dO14_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real F17_log_pf, dF17_log_pf_dT9;

    // interpolating F17 partition function
    get_partition_function_cached(F17, tfactors, pf_cache, F17_log_pf, dF17_log_pf_dT9);

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
    get_partition_function_cached(F17, tfactors, pf_cache, F17_log_pf, dF17_log_pf_dT9);

    amrex::Real He4_log_pf, dHe4_log_pf_dT9;

    // setting He4 log(partition function) to 0.0 by default, independent of T
    He4_log_pf = 0.0_rt;
    dHe4_log_pf_dT9 = 0.0_rt;

    amrex::Real Ne20_log_pf, dNe20_log_pf_dT9;

    // interpolating Ne20 partition function
    get_partition_function_cached(Ne20, tfactors, pf_cache, Ne20_log_pf, dNe20_log_pf_dT9);

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
void rate_He4_Ne20_to_C12_C12_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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
    get_partition_function_cached(Ne20, tfactors, pf_cache, Ne20_log_pf, dNe20_log_pf_dT9);

    amrex::Real net_log_pf = C12_log_pf + C12_log_pf - He4_log_pf - Ne20_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dC12_log_pf_dT9 + dC12_log_pf_dT9 - dHe4_log_pf_dT9 - dNe20_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from cf88r
    ln_set_rate =  61.474151468919175 + -53.57824966896331 * tfactors.T9i + -84.165 * tfactors.T913i + -1.56627 * tfactors.T913
                         + -0.0736084 * tfactors.T9 + -0.072797 * tfactors.T953 + -0.666667 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  53.57824966896331 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -84.165 * tfactors.T943i + (1.0/3.0) * -1.56627 * tfactors.T923i
                                  + -0.0736084 + (5.0/3.0) * -0.072797 * tfactors.T923 + -0.666667 * tfactors.T9i;

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
void rate_He4_Mg24_to_C12_O16_derived(const tf_t& tfactors, amrex::Real& rate, amrex::Real& drate_dT, [[maybe_unused]] const T& rate_eval, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

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
    get_partition_function_cached(Mg24, tfactors, pf_cache, Mg24_log_pf, dMg24_log_pf_dT9);

    amrex::Real O16_log_pf, dO16_log_pf_dT9;

    // interpolating O16 partition function
    get_partition_function_cached(O16, tfactors, pf_cache, O16_log_pf, dO16_log_pf_dT9);

    amrex::Real net_log_pf = C12_log_pf + O16_log_pf - He4_log_pf - Mg24_log_pf;
    [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dC12_log_pf_dT9 + dO16_log_pf_dT9 - dHe4_log_pf_dT9 - dMg24_log_pf_dT9;

    amrex::Real ln_set_rate{0.0};
    amrex::Real dln_set_rate_dT9{0.0};
    amrex::Real set_rate{0.0};

    // ReacLib set derived from cf88r
    ln_set_rate =  49.57295819231048 + -78.2098115311048 * tfactors.T9i + -133.413 * tfactors.T913i + 50.1572 * tfactors.T913
                         + -3.15987 * tfactors.T9 + 0.0178251 * tfactors.T953 + -23.7027 * tfactors.lnT9;

    ln_set_rate += net_log_pf;

    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        dln_set_rate_dT9 =  78.2098115311048 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -133.413 * tfactors.T943i + (1.0/3.0) * 50.1572 * tfactors.T923i
                                  + -3.15987 + (5.0/3.0) * 0.0178251 * tfactors.T923 + -23.7027 * tfactors.T9i;

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
    rate_Mg22_to_He4_Ne18_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Mg22_to_He4_Ne18_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Mg22_to_He4_Ne18_derived) = drate_dT;

    }
    rate_Mg24_to_He4_Ne20_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_Mg24_to_He4_Ne20_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Mg24_to_He4_Ne20_derived) = drate_dT;

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
    rate_He4_Ne20_to_C12_C12_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_He4_Ne20_to_C12_C12_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Ne20_to_C12_C12_derived) = drate_dT;

    }
    rate_He4_Mg24_to_C12_O16_derived<do_T_derivatives, T>(tfactors, rate, drate_dT, rate_eval, pf_cache);
    rate_eval.screened_rates(k_He4_Mg24_to_C12_O16_derived) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_He4_Mg24_to_C12_O16_derived) = drate_dT;

    }

}
#endif
