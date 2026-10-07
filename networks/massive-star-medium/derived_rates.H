#ifndef DERIVED_RATES_H
#define DERIVED_RATES_H

#include <tfactors.H>
#include <actual_network.H>
#include <partition_functions.H>
#include <temperature_table_rates.H>
#include <interp_tools.H>

namespace derived_rates {

    using namespace Rates;

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_O16_to_He4_C12_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // O16 --> He4 + C12

        amrex::Real rate{}, drate_dT{};
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_O16_to_He4_C12_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_O16_to_He4_C12_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ne20_to_He4_O16_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ne20 --> He4 + O16

        amrex::Real rate{}, drate_dT{};
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ne20_to_He4_O16_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne20_to_He4_O16_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mg24_to_n_Mg23_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Mg24 --> n + Mg23

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Mg23_log_pf, dMg23_log_pf_dT9;

        // interpolating Mg23 partition function
        get_partition_function_cached(Mg23, tfactors.T9, pf_cache, Mg23_log_pf, dMg23_log_pf_dT9);

        amrex::Real Mg24_log_pf, dMg24_log_pf_dT9;

        // interpolating Mg24 partition function
        get_partition_function_cached(Mg24, tfactors.T9, pf_cache, Mg24_log_pf, dMg24_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Mg23_log_pf - Mg24_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dMg23_log_pf_dT9 - dMg24_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  32.05683058133977 + -191.83691284889673 * tfactors.T9i + 2.66964 * tfactors.T913
                             + -0.448904 * tfactors.T9 + 0.0326505 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  191.83691284889673 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 2.66964 * tfactors.T923i
                                      + -0.448904 + (5.0/3.0) * 0.0326505 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Mg24_to_n_Mg23_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg24_to_n_Mg23_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si27_to_He4_Mg23_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Si27 --> He4 + Mg23

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Mg23_log_pf, dMg23_log_pf_dT9;

        // interpolating Mg23 partition function
        get_partition_function_cached(Mg23, tfactors.T9, pf_cache, Mg23_log_pf, dMg23_log_pf_dT9);

        amrex::Real Si27_log_pf, dSi27_log_pf_dT9;

        // interpolating Si27 partition function
        get_partition_function_cached(Si27, tfactors.T9, pf_cache, Si27_log_pf, dSi27_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Mg23_log_pf - Si27_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dMg23_log_pf_dT9 - dSi27_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  68.14909115483175 + -108.33694817182082 * tfactors.T9i + -53.203 * tfactors.T913i + -4.6318 * tfactors.T913
                             + -0.130951 * tfactors.T9 + 0.014691 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  108.33694817182082 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -53.203 * tfactors.T943i + (1.0/3.0) * -4.6318 * tfactors.T923i
                                      + -0.130951 + (5.0/3.0) * 0.014691 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Si27_to_He4_Mg23_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si27_to_He4_Mg23_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si28_to_n_Si27_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Si28 --> n + Si27

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Si27_log_pf, dSi27_log_pf_dT9;

        // interpolating Si27 partition function
        get_partition_function_cached(Si27, tfactors.T9, pf_cache, Si27_log_pf, dSi27_log_pf_dT9);

        amrex::Real Si28_log_pf, dSi28_log_pf_dT9;

        // interpolating Si28 partition function
        get_partition_function_cached(Si28, tfactors.T9, pf_cache, Si28_log_pf, dSi28_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Si27_log_pf - Si28_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dSi27_log_pf_dT9 - dSi28_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  34.60590343256573 + -199.36103986444223 * tfactors.T9i + 0.559804 * tfactors.T913
                             + -0.0233808 * tfactors.T9 + -0.00121313 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  199.36103986444223 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 0.559804 * tfactors.T923i
                                      + -0.0233808 + (5.0/3.0) * -0.00121313 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Si28_to_n_Si27_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si28_to_n_Si27_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si28_to_p_Al27_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Si28 --> p + Al27

        amrex::Real rate{}, drate_dT{};
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Si28_to_p_Al27_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si28_to_p_Al27_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si28_to_He4_Mg24_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Si28 --> He4 + Mg24

        amrex::Real rate{}, drate_dT{};
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Si28_to_He4_Mg24_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si28_to_He4_Mg24_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P29_to_p_Si28_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // P29 --> p + Si28

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real P29_log_pf, dP29_log_pf_dT9;

        // interpolating P29 partition function
        get_partition_function_cached(P29, tfactors.T9, pf_cache, P29_log_pf, dP29_log_pf_dT9);

        amrex::Real Si28_log_pf, dSi28_log_pf_dT9;

        // interpolating Si28 partition function
        get_partition_function_cached(Si28, tfactors.T9, pf_cache, Si28_log_pf, dSi28_log_pf_dT9);

        amrex::Real net_log_pf = p_log_pf + Si28_log_pf - P29_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dSi28_log_pf_dT9 - dP29_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from il10n
        ln_set_rate =  39.156867542086886 + -31.90051818093998 * tfactors.T9i + -23.8173 * tfactors.T913i + 7.08203 * tfactors.T913
                             + -1.44753 * tfactors.T9 + 0.0804296 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  31.90051818093998 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -23.8173 * tfactors.T943i + (1.0/3.0) * 7.08203 * tfactors.T923i
                                      + -1.44753 + (5.0/3.0) * 0.0804296 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from il10r
        ln_set_rate =  28.71871754208689 + -36.04283818093998 * tfactors.T9i;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  36.04283818093998 * tfactors.T9i * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_P29_to_p_Si28_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P29_to_p_Si28_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P31_to_p_Si30_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // P31 --> p + Si30

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real Si30_log_pf, dSi30_log_pf_dT9;

        // interpolating Si30 partition function
        get_partition_function_cached(Si30, tfactors.T9, pf_cache, Si30_log_pf, dSi30_log_pf_dT9);

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real P31_log_pf, dP31_log_pf_dT9;

        // interpolating P31 partition function
        get_partition_function_cached(P31, tfactors.T9, pf_cache, P31_log_pf, dP31_log_pf_dT9);

        amrex::Real net_log_pf = p_log_pf + Si30_log_pf - P31_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dSi30_log_pf_dT9 - dP31_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from de20r
        ln_set_rate =  35.352341431948254 + -91.31403689784659 * tfactors.T9i
                             + 0.3809 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  91.31403689784659 * tfactors.T9i * tfactors.T9i
                                      + 0.3809 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from de20r
        ln_set_rate =  31.780501431948256 + -89.85529689784659 * tfactors.T9i
                             + 2.7883 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  89.85529689784659 * tfactors.T9i * tfactors.T9i
                                      + 2.7883 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from de20r
        ln_set_rate =  -311.28315856805176 + -85.80625689784658 * tfactors.T9i
                             + -77.047 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  85.80625689784658 * tfactors.T9i * tfactors.T9i
                                      + -77.047 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from de20r
        ln_set_rate =  18.740761431948258 + -85.92472689784658 * tfactors.T9i
                             + 6.49034 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  85.92472689784658 * tfactors.T9i * tfactors.T9i
                                      + 6.49034 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from de20r
        ln_set_rate =  4.063041431948257 + -85.61824789784659 * tfactors.T9i
                             + 2.8033099999999997 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  85.61824789784659 * tfactors.T9i * tfactors.T9i
                                      + 2.8033099999999997 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from de20r
        ln_set_rate =  -1115.3571585680518 + -180.5498868978466 * tfactors.T9i
                             + -895.258 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  180.5498868978466 * tfactors.T9i * tfactors.T9i
                                      + -895.258 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from de20r
        ln_set_rate =  32.94828143194825 + -90.26261689784658 * tfactors.T9i
                             + -0.0708200000000001 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  90.26261689784658 * tfactors.T9i * tfactors.T9i
                                      + -0.0708200000000001 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from de20r
        ln_set_rate =  -11.876558568051745 + -85.26592089784658 * tfactors.T9i
                             + -0.12839 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  85.26592089784658 * tfactors.T9i * tfactors.T9i
                                      + -0.12839 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_P31_to_p_Si30_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P31_to_p_Si30_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P31_to_He4_Al27_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // P31 --> He4 + Al27

        amrex::Real rate{}, drate_dT{};
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_P31_to_He4_Al27_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P31_to_He4_Al27_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S30_to_p_P29_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // S30 --> p + P29

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real S30_log_pf, dS30_log_pf_dT9;

        // interpolating S30 partition function
        get_partition_function_cached(S30, tfactors.T9, pf_cache, S30_log_pf, dS30_log_pf_dT9);

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real P29_log_pf, dP29_log_pf_dT9;

        // interpolating P29 partition function
        get_partition_function_cached(P29, tfactors.T9, pf_cache, P29_log_pf, dP29_log_pf_dT9);

        amrex::Real net_log_pf = p_log_pf + P29_log_pf - S30_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dP29_log_pf_dT9 - dS30_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from il10r
        ln_set_rate =  -1.6327732225774767 + -54.72612338903181 * tfactors.T9i + 8.48834 * tfactors.T913i + 25.65 * tfactors.T913
                             + -3.79773 * tfactors.T9 + 0.320391 * tfactors.T953;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  54.72612338903181 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 8.48834 * tfactors.T943i + (1.0/3.0) * 25.65 * tfactors.T923i
                                      + -3.79773 + (5.0/3.0) * 0.320391 * tfactors.T923;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from il10n
        ln_set_rate =  45.34032677742252 + -51.00674338903181 * tfactors.T9i + -25.6007 * tfactors.T913i
                             + -2.0 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  51.00674338903181 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -25.6007 * tfactors.T943i
                                      + (5.0/3.0) * -2.0 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_S30_to_p_P29_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S30_to_p_P29_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S32_to_p_P31_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // S32 --> p + P31

        amrex::Real rate{}, drate_dT{};
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_S32_to_p_P31_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S32_to_p_P31_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S32_to_He4_Si28_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // S32 --> He4 + Si28

        amrex::Real rate{}, drate_dT{};
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_S32_to_He4_Si28_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S32_to_He4_Si28_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S34_to_He4_Si30_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // S34 --> He4 + Si30

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real Si30_log_pf, dSi30_log_pf_dT9;

        // interpolating Si30 partition function
        get_partition_function_cached(Si30, tfactors.T9, pf_cache, Si30_log_pf, dSi30_log_pf_dT9);

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real S34_log_pf, dS34_log_pf_dT9;

        // interpolating S34 partition function
        get_partition_function_cached(S34, tfactors.T9, pf_cache, S34_log_pf, dS34_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Si30_log_pf - S34_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dSi30_log_pf_dT9 - dS34_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  72.40691003335057 + -91.95006887812762 * tfactors.T9i + -59.6559 * tfactors.T913i + 3.70141 * tfactors.T913
                             + -3.12537 * tfactors.T9 + 0.307626 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  91.95006887812762 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -59.6559 * tfactors.T943i + (1.0/3.0) * 3.70141 * tfactors.T923i
                                      + -3.12537 + (5.0/3.0) * 0.307626 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_S34_to_He4_Si30_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S34_to_He4_Si30_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl33_to_p_S32_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Cl33 --> p + S32

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real Cl33_log_pf, dCl33_log_pf_dT9;

        // interpolating Cl33 partition function
        get_partition_function_cached(Cl33, tfactors.T9, pf_cache, Cl33_log_pf, dCl33_log_pf_dT9);

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real S32_log_pf, dS32_log_pf_dT9;

        // interpolating S32 partition function
        get_partition_function_cached(S32, tfactors.T9, pf_cache, S32_log_pf, dS32_log_pf_dT9);

        amrex::Real net_log_pf = p_log_pf + S32_log_pf - Cl33_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dS32_log_pf_dT9 - dCl33_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from il10n
        ln_set_rate =  74.7623172827099 + -26.42040054322257 * tfactors.T9i + -29.7741 * tfactors.T913i + -87.4473 * tfactors.T913
                             + 182.189 * tfactors.T9 + -128.625 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  26.42040054322257 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -29.7741 * tfactors.T943i + (1.0/3.0) * -87.4473 * tfactors.T923i
                                      + 182.189 + (5.0/3.0) * -128.625 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from il10r
        ln_set_rate =  91.6269172827099 + -29.424110543222568 * tfactors.T9i + -33.7204 * tfactors.T913i + -32.7355 * tfactors.T913
                             + 3.92526 * tfactors.T9 + -0.250479 * tfactors.T953;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  29.424110543222568 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -33.7204 * tfactors.T943i + (1.0/3.0) * -32.7355 * tfactors.T923i
                                      + 3.92526 + (5.0/3.0) * -0.250479 * tfactors.T923;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from il10r
        ln_set_rate =  -4.945882717290097 + -27.29450754322257 * tfactors.T9i;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  27.29450754322257 * tfactors.T9i * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Cl33_to_p_S32_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl33_to_p_S32_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl33_to_He4_P29_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Cl33 --> He4 + P29

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real Cl33_log_pf, dCl33_log_pf_dT9;

        // interpolating Cl33 partition function
        get_partition_function_cached(Cl33, tfactors.T9, pf_cache, Cl33_log_pf, dCl33_log_pf_dT9);

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real P29_log_pf, dP29_log_pf_dT9;

        // interpolating P29 partition function
        get_partition_function_cached(P29, tfactors.T9, pf_cache, P29_log_pf, dP29_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + P29_log_pf - Cl33_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dP29_log_pf_dT9 - dCl33_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  66.20475225259291 + -75.14408080801772 * tfactors.T9i + -62.3802 * tfactors.T913i + 0.593062 * tfactors.T913
                             + -1.14226 * tfactors.T9 + 0.0934776 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  75.14408080801772 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -62.3802 * tfactors.T943i + (1.0/3.0) * 0.593062 * tfactors.T923i
                                      + -1.14226 + (5.0/3.0) * 0.0934776 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Cl33_to_He4_P29_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl33_to_He4_P29_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl35_to_p_S34_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Cl35 --> p + S34

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real Cl35_log_pf, dCl35_log_pf_dT9;

        // interpolating Cl35 partition function
        get_partition_function_cached(Cl35, tfactors.T9, pf_cache, Cl35_log_pf, dCl35_log_pf_dT9);

        amrex::Real S34_log_pf, dS34_log_pf_dT9;

        // interpolating S34 partition function
        get_partition_function_cached(S34, tfactors.T9, pf_cache, S34_log_pf, dS34_log_pf_dT9);

        amrex::Real net_log_pf = p_log_pf + S34_log_pf - Cl35_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dS34_log_pf_dT9 - dCl35_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from se19r
        ln_set_rate =  14.594869555370497 + -75.95634244116214 * tfactors.T9i;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  75.95634244116214 * tfactors.T9i * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from se19r
        ln_set_rate =  30.549759555370496 + -79.68602244116214 * tfactors.T9i;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  79.68602244116214 * tfactors.T9i * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from se19r
        ln_set_rate =  7.417229555370497 + -75.33433244116213 * tfactors.T9i;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  75.33433244116213 * tfactors.T9i * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from se19r
        ln_set_rate =  29.296749555370496 + -79.03617244116214 * tfactors.T9i;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  79.03617244116214 * tfactors.T9i * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from se19r
        ln_set_rate =  4.358729555370498 + -75.10224244116213 * tfactors.T9i;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  75.10224244116213 * tfactors.T9i * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from se19r
        ln_set_rate =  28.101369555370496 + -78.57199244116214 * tfactors.T9i;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  78.57199244116214 * tfactors.T9i * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from se19r
        ln_set_rate =  22.729819555370497 + -77.09822244116214 * tfactors.T9i;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  77.09822244116214 * tfactors.T9i * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from se19r
        ln_set_rate =  -37.72937044462951 + -101.61919244116214 * tfactors.T9i + 454.53 * tfactors.T913i + -375.557 * tfactors.T913
                             + 12.6533 * tfactors.T9 + -0.408677 * tfactors.T953 + 246.618 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  101.61919244116214 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 454.53 * tfactors.T943i + (1.0/3.0) * -375.557 * tfactors.T923i
                                      + 12.6533 + (5.0/3.0) * -0.408677 * tfactors.T923 + 246.618 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Cl35_to_p_S34_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl35_to_p_S34_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl35_to_He4_P31_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Cl35 --> He4 + P31

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  81.20727442147438 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -62.5433 * tfactors.T943i + (1.0/3.0) * -2.95026 * tfactors.T923i
                                      + -0.89652 + (5.0/3.0) * 0.0774126 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Cl35_to_He4_P31_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl35_to_He4_P31_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar34_to_p_Cl33_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ar34 --> p + Cl33

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real Cl33_log_pf, dCl33_log_pf_dT9;

        // interpolating Cl33 partition function
        get_partition_function_cached(Cl33, tfactors.T9, pf_cache, Cl33_log_pf, dCl33_log_pf_dT9);

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real Ar34_log_pf, dAr34_log_pf_dT9;

        // interpolating Ar34 partition function
        get_partition_function_cached(Ar34, tfactors.T9, pf_cache, Ar34_log_pf, dAr34_log_pf_dT9);

        amrex::Real net_log_pf = p_log_pf + Cl33_log_pf - Ar34_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dCl33_log_pf_dT9 - dAr34_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  60.19620056329993 + -54.12302068534123 * tfactors.T9i + -27.8815 * tfactors.T913i + -3.18731 * tfactors.T913
                             + -1.76254 * tfactors.T9 + 0.264735 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  54.12302068534123 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -27.8815 * tfactors.T943i + (1.0/3.0) * -3.18731 * tfactors.T923i
                                      + -1.76254 + (5.0/3.0) * 0.264735 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ar34_to_p_Cl33_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar34_to_p_Cl33_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar34_to_He4_S30_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ar34 --> He4 + S30

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real S30_log_pf, dS30_log_pf_dT9;

        // interpolating S30 partition function
        get_partition_function_cached(S30, tfactors.T9, pf_cache, S30_log_pf, dS30_log_pf_dT9);

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Ar34_log_pf, dAr34_log_pf_dT9;

        // interpolating Ar34 partition function
        get_partition_function_cached(Ar34, tfactors.T9, pf_cache, Ar34_log_pf, dAr34_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + S30_log_pf - Ar34_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dS30_log_pf_dT9 - dAr34_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  68.01882603847031 + -78.26035810433434 * tfactors.T9i + -65.211 * tfactors.T913i + -1.41447 * tfactors.T913
                             + -0.542976 * tfactors.T9 + 0.0211165 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  78.26035810433434 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -65.211 * tfactors.T943i + (1.0/3.0) * -1.41447 * tfactors.T923i
                                      + -0.542976 + (5.0/3.0) * 0.0211165 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ar34_to_He4_S30_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar34_to_He4_S30_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar36_to_p_Cl35_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ar36 --> p + Cl35

        amrex::Real rate{}, drate_dT{};
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ar36_to_p_Cl35_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar36_to_p_Cl35_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar36_to_He4_S32_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ar36 --> He4 + S32

        amrex::Real rate{}, drate_dT{};
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ar36_to_He4_S32_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar36_to_He4_S32_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar38_to_He4_S34_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ar38 --> He4 + S34

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real Ar38_log_pf, dAr38_log_pf_dT9;

        // interpolating Ar38 partition function
        get_partition_function_cached(Ar38, tfactors.T9, pf_cache, Ar38_log_pf, dAr38_log_pf_dT9);

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real S34_log_pf, dS34_log_pf_dT9;

        // interpolating S34 partition function
        get_partition_function_cached(S34, tfactors.T9, pf_cache, S34_log_pf, dS34_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + S34_log_pf - Ar38_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dS34_log_pf_dT9 - dAr38_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from cn18
        ln_set_rate =  253.65592843778612 + -101.61351496454262 * tfactors.T9i + 635.332 * tfactors.T913i + -902.05 * tfactors.T913
                             + 42.5623 * tfactors.T9 + -2.08434 * tfactors.T953 + 478.85 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  101.61351496454262 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 635.332 * tfactors.T943i + (1.0/3.0) * -902.05 * tfactors.T923i
                                      + 42.5623 + (5.0/3.0) * -2.08434 * tfactors.T923 + 478.85 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ar38_to_He4_S34_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar38_to_He4_S34_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K37_to_p_Ar36_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // K37 --> p + Ar36

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real K37_log_pf, dK37_log_pf_dT9;

        // interpolating K37 partition function
        get_partition_function_cached(K37, tfactors.T9, pf_cache, K37_log_pf, dK37_log_pf_dT9);

        amrex::Real Ar36_log_pf, dAr36_log_pf_dT9;

        // interpolating Ar36 partition function
        get_partition_function_cached(Ar36, tfactors.T9, pf_cache, Ar36_log_pf, dAr36_log_pf_dT9);

        amrex::Real net_log_pf = p_log_pf + Ar36_log_pf - K37_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dAr36_log_pf_dT9 - dK37_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from il10r
        ln_set_rate =  25.48017518360185 + -30.537200136219777 * tfactors.T9i + 7.50962 * tfactors.T913
                             + -1.68675 * tfactors.T9 + 0.189725 * tfactors.T953;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  30.537200136219777 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 7.50962 * tfactors.T923i
                                      + -1.68675 + (5.0/3.0) * 0.189725 * tfactors.T923;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from il10r
        ln_set_rate =  26.99508518360185 + -25.181990136219778 * tfactors.T9i;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  25.181990136219778 * tfactors.T9i * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from il10n
        ln_set_rate =  43.78969518360185 + -21.556890136219778 * tfactors.T9i + -28.9682 * tfactors.T913i
                             + -1.0 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  21.556890136219778 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -28.9682 * tfactors.T943i
                                      + (5.0/3.0) * -1.0 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_K37_to_p_Ar36_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K37_to_p_Ar36_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K37_to_He4_Cl33_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // K37 --> He4 + Cl33

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real Cl33_log_pf, dCl33_log_pf_dT9;

        // interpolating Cl33 partition function
        get_partition_function_cached(Cl33, tfactors.T9, pf_cache, Cl33_log_pf, dCl33_log_pf_dT9);

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real K37_log_pf, dK37_log_pf_dT9;

        // interpolating K37 partition function
        get_partition_function_cached(K37, tfactors.T9, pf_cache, K37_log_pf, dK37_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Cl33_log_pf - K37_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCl33_log_pf_dT9 - dK37_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  71.60805297248307 + -72.20117501230614 * tfactors.T9i + -68.1442 * tfactors.T913i + -2.72746 * tfactors.T913
                             + -0.317537 * tfactors.T9 + 0.0162782 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  72.20117501230614 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -68.1442 * tfactors.T943i + (1.0/3.0) * -2.72746 * tfactors.T923i
                                      + -0.317537 + (5.0/3.0) * 0.0162782 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_K37_to_He4_Cl33_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K37_to_He4_Cl33_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K39_to_p_Ar38_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // K39 --> p + Ar38

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real Ar38_log_pf, dAr38_log_pf_dT9;

        // interpolating Ar38 partition function
        get_partition_function_cached(Ar38, tfactors.T9, pf_cache, Ar38_log_pf, dAr38_log_pf_dT9);

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real K39_log_pf, dK39_log_pf_dT9;

        // interpolating K39 partition function
        get_partition_function_cached(K39, tfactors.T9, pf_cache, K39_log_pf, dK39_log_pf_dT9);

        amrex::Real net_log_pf = p_log_pf + Ar38_log_pf - K39_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dAr38_log_pf_dT9 - dK39_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  57.583260221976474 + -74.05232999435097 * tfactors.T9i + -29.0021 * tfactors.T913i + -0.525968 * tfactors.T913
                             + -1.94216 * tfactors.T9 + 0.267346 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  74.05232999435097 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -29.0021 * tfactors.T943i + (1.0/3.0) * -0.525968 * tfactors.T923i
                                      + -1.94216 + (5.0/3.0) * 0.267346 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_K39_to_p_Ar38_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K39_to_p_Ar38_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K39_to_He4_Cl35_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // K39 --> He4 + Cl35

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  83.76815251770985 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -68.2848 * tfactors.T943i + (1.0/3.0) * 0.0178545 * tfactors.T923i
                                      + -2.06783 + (5.0/3.0) * 0.199659 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_K39_to_He4_Cl35_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K39_to_He4_Cl35_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca38_to_p_K37_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ca38 --> p + K37

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real K37_log_pf, dK37_log_pf_dT9;

        // interpolating K37 partition function
        get_partition_function_cached(K37, tfactors.T9, pf_cache, K37_log_pf, dK37_log_pf_dT9);

        amrex::Real Ca38_log_pf, dCa38_log_pf_dT9;

        // interpolating Ca38 partition function
        get_partition_function_cached(Ca38, tfactors.T9, pf_cache, Ca38_log_pf, dCa38_log_pf_dT9);

        amrex::Real net_log_pf = p_log_pf + K37_log_pf - Ca38_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dK37_log_pf_dT9 - dCa38_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  61.28349410873902 + -52.76888946578814 * tfactors.T9i + -30.0596 * tfactors.T913i + -2.56332 * tfactors.T913
                             + -2.026 * tfactors.T9 + 0.289027 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  52.76888946578814 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -30.0596 * tfactors.T943i + (1.0/3.0) * -2.56332 * tfactors.T923i
                                      + -2.026 + (5.0/3.0) * 0.289027 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ca38_to_p_K37_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca38_to_p_K37_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca38_to_He4_Ar34_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ca38 --> He4 + Ar34

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Ar34_log_pf, dAr34_log_pf_dT9;

        // interpolating Ar34 partition function
        get_partition_function_cached(Ar34, tfactors.T9, pf_cache, Ar34_log_pf, dAr34_log_pf_dT9);

        amrex::Real Ca38_log_pf, dCa38_log_pf_dT9;

        // interpolating Ca38 partition function
        get_partition_function_cached(Ca38, tfactors.T9, pf_cache, Ca38_log_pf, dCa38_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Ar34_log_pf - Ca38_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dAr34_log_pf_dT9 - dCa38_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  73.80824651792216 + -70.84704379274585 * tfactors.T9i + -70.8666 * tfactors.T913i + -5.98563 * tfactors.T913
                             + 0.582779 * tfactors.T9 + -0.093842 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  70.84704379274585 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -70.8666 * tfactors.T943i + (1.0/3.0) * -5.98563 * tfactors.T923i
                                      + 0.582779 + (5.0/3.0) * -0.093842 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ca38_to_He4_Ar34_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca38_to_He4_Ar34_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca40_to_p_K39_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ca40 --> p + K39

        amrex::Real rate{}, drate_dT{};
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ca40_to_p_K39_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca40_to_p_K39_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca40_to_He4_Ar36_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ca40 --> He4 + Ar36

        amrex::Real rate{}, drate_dT{};
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ca40_to_He4_Ar36_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca40_to_He4_Ar36_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca42_to_He4_Ar38_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ca42 --> He4 + Ar38

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real Ar38_log_pf, dAr38_log_pf_dT9;

        // interpolating Ar38 partition function
        get_partition_function_cached(Ar38, tfactors.T9, pf_cache, Ar38_log_pf, dAr38_log_pf_dT9);

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Ca42_log_pf, dCa42_log_pf_dT9;

        // interpolating Ca42 partition function
        get_partition_function_cached(Ca42, tfactors.T9, pf_cache, Ca42_log_pf, dCa42_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Ar38_log_pf - Ca42_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dAr38_log_pf_dT9 - dCa42_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  77.58296390367897 + -72.6138316767612 * tfactors.T9i + -71.1296 * tfactors.T913i + 2.75299 * tfactors.T913
                             + -2.43406 * tfactors.T9 + 0.199511 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  72.6138316767612 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -71.1296 * tfactors.T943i + (1.0/3.0) * 2.75299 * tfactors.T923i
                                      + -2.43406 + (5.0/3.0) * 0.199511 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ca42_to_He4_Ar38_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca42_to_He4_Ar38_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc41_to_p_Ca40_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Sc41 --> p + Ca40

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real Ca40_log_pf, dCa40_log_pf_dT9;

        // interpolating Ca40 partition function
        get_partition_function_cached(Ca40, tfactors.T9, pf_cache, Ca40_log_pf, dCa40_log_pf_dT9);

        amrex::Real Sc41_log_pf, dSc41_log_pf_dT9;

        // interpolating Sc41 partition function
        get_partition_function_cached(Sc41, tfactors.T9, pf_cache, Sc41_log_pf, dSc41_log_pf_dT9);

        amrex::Real net_log_pf = p_log_pf + Ca40_log_pf - Sc41_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dCa40_log_pf_dT9 - dSc41_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from il10n
        ln_set_rate =  41.615999266911004 + -12.590078983707985 * tfactors.T9i + -31.1711 * tfactors.T913i
                             + -1.0 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  12.590078983707985 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -31.1711 * tfactors.T943i
                                      + (5.0/3.0) * -1.0 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from il10r
        ln_set_rate =  72.295299266911 + -46.102378983707986 * tfactors.T9i + -33.7287 * tfactors.T913
                             + 5.22449 * tfactors.T9 + -0.357104 * tfactors.T953;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  46.102378983707986 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -33.7287 * tfactors.T923i
                                      + 5.22449 + (5.0/3.0) * -0.357104 * tfactors.T923;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from il10r
        ln_set_rate =  27.237989266911 + -19.912438983707986 * tfactors.T9i;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  19.912438983707986 * tfactors.T9i * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Sc41_to_p_Ca40_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc41_to_p_Ca40_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc41_to_He4_K37_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Sc41 --> He4 + K37

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real K37_log_pf, dK37_log_pf_dT9;

        // interpolating K37 partition function
        get_partition_function_cached(K37, tfactors.T9, pf_cache, K37_log_pf, dK37_log_pf_dT9);

        amrex::Real Sc41_log_pf, dSc41_log_pf_dT9;

        // interpolating Sc41 partition function
        get_partition_function_cached(Sc41, tfactors.T9, pf_cache, Sc41_log_pf, dSc41_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + K37_log_pf - Sc41_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dK37_log_pf_dT9 - dSc41_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  74.8476253358184 + -72.72639550252572 * tfactors.T9i + -73.6776 * tfactors.T913i + -5.3733 * tfactors.T913
                             + 0.685776 * tfactors.T9 + -0.123735 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  72.72639550252572 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -73.6776 * tfactors.T943i + (1.0/3.0) * -5.3733 * tfactors.T923i
                                      + 0.685776 + (5.0/3.0) * -0.123735 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Sc41_to_He4_K37_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc41_to_He4_K37_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc43_to_p_Ca42_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Sc43 --> p + Ca42

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real Ca42_log_pf, dCa42_log_pf_dT9;

        // interpolating Ca42 partition function
        get_partition_function_cached(Ca42, tfactors.T9, pf_cache, Ca42_log_pf, dCa42_log_pf_dT9);

        amrex::Real Sc43_log_pf, dSc43_log_pf_dT9;

        // interpolating Sc43 partition function
        get_partition_function_cached(Sc43, tfactors.T9, pf_cache, Sc43_log_pf, dSc43_log_pf_dT9);

        amrex::Real net_log_pf = p_log_pf + Ca42_log_pf - Sc43_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dCa42_log_pf_dT9 - dSc43_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  57.429501245652034 + -57.20773369244204 * tfactors.T9i + -31.1381 * tfactors.T913i + -0.151667 * tfactors.T913
                             + -1.76926 * tfactors.T9 + 0.219569 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  57.20773369244204 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -31.1381 * tfactors.T943i + (1.0/3.0) * -0.151667 * tfactors.T923i
                                      + -1.76926 + (5.0/3.0) * 0.219569 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Sc43_to_p_Ca42_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc43_to_p_Ca42_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc43_to_He4_K39_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Sc43 --> He4 + K39

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  55.76923537484507 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -73.8006 * tfactors.T943i + (1.0/3.0) * 1.87885 * tfactors.T923i
                                      + -2.75862 + (5.0/3.0) * 0.279678 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Sc43_to_He4_K39_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc43_to_He4_K39_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti42_to_p_Sc41_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ti42 --> p + Sc41

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real Ti42_log_pf, dTi42_log_pf_dT9;

        // interpolating Ti42 partition function
        get_partition_function_cached(Ti42, tfactors.T9, pf_cache, Ti42_log_pf, dTi42_log_pf_dT9);

        amrex::Real Sc41_log_pf, dSc41_log_pf_dT9;

        // interpolating Sc41 partition function
        get_partition_function_cached(Sc41, tfactors.T9, pf_cache, Sc41_log_pf, dSc41_log_pf_dT9);

        amrex::Real net_log_pf = p_log_pf + Sc41_log_pf - Ti42_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dSc41_log_pf_dT9 - dTi42_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  64.12066756154039 + -43.52809564043815 * tfactors.T9i + -32.1615 * tfactors.T913i + -3.68811 * tfactors.T913
                             + -2.01234 * tfactors.T9 + 0.3232 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  43.52809564043815 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -32.1615 * tfactors.T943i + (1.0/3.0) * -3.68811 * tfactors.T923i
                                      + -2.01234 + (5.0/3.0) * 0.3232 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ti42_to_p_Sc41_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti42_to_p_Sc41_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti42_to_He4_Ca38_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ti42 --> He4 + Ca38

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Ti42_log_pf, dTi42_log_pf_dT9;

        // interpolating Ti42 partition function
        get_partition_function_cached(Ti42, tfactors.T9, pf_cache, Ti42_log_pf, dTi42_log_pf_dT9);

        amrex::Real Ca38_log_pf, dCa38_log_pf_dT9;

        // interpolating Ca38 partition function
        get_partition_function_cached(Ca38, tfactors.T9, pf_cache, Ca38_log_pf, dCa38_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Ca38_log_pf - Ti42_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCa38_log_pf_dT9 - dTi42_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  83.43049878861977 + -63.48560167717572 * tfactors.T9i + -76.3063 * tfactors.T913i + -10.7213 * tfactors.T913
                             + 0.926161 * tfactors.T9 + -0.0855699 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  63.48560167717572 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -76.3063 * tfactors.T943i + (1.0/3.0) * -10.7213 * tfactors.T923i
                                      + 0.926161 + (5.0/3.0) * -0.0855699 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ti42_to_He4_Ca38_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti42_to_He4_Ca38_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti44_to_p_Sc43_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ti44 --> p + Sc43

        amrex::Real rate{}, drate_dT{};
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ti44_to_p_Sc43_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti44_to_p_Sc43_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti44_to_He4_Ca40_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ti44 --> He4 + Ca40

        amrex::Real rate{}, drate_dT{};
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ti44_to_He4_Ca40_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti44_to_He4_Ca40_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti46_to_He4_Ca42_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ti46 --> He4 + Ca42

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real Ti46_log_pf, dTi46_log_pf_dT9;

        // interpolating Ti46 partition function
        get_partition_function_cached(Ti46, tfactors.T9, pf_cache, Ti46_log_pf, dTi46_log_pf_dT9);

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Ca42_log_pf, dCa42_log_pf_dT9;

        // interpolating Ca42 partition function
        get_partition_function_cached(Ca42, tfactors.T9, pf_cache, Ca42_log_pf, dCa42_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Ca42_log_pf - Ti46_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCa42_log_pf_dT9 - dTi46_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  79.79033156127323 + -92.90456370262497 * tfactors.T9i + -76.5379 * tfactors.T913i + 2.80099 * tfactors.T913
                             + -2.34788 * tfactors.T9 + 0.150014 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  92.90456370262497 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -76.5379 * tfactors.T943i + (1.0/3.0) * 2.80099 * tfactors.T923i
                                      + -2.34788 + (5.0/3.0) * 0.150014 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ti46_to_He4_Ca42_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti46_to_He4_Ca42_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V45_to_p_Ti44_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // V45 --> p + Ti44

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real V45_log_pf, dV45_log_pf_dT9;

        // interpolating V45 partition function
        get_partition_function_cached(V45, tfactors.T9, pf_cache, V45_log_pf, dV45_log_pf_dT9);

        amrex::Real Ti44_log_pf, dTi44_log_pf_dT9;

        // interpolating Ti44 partition function
        get_partition_function_cached(Ti44, tfactors.T9, pf_cache, Ti44_log_pf, dTi44_log_pf_dT9);

        amrex::Real net_log_pf = p_log_pf + Ti44_log_pf - V45_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dTi44_log_pf_dT9 - dV45_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from nfisn
        ln_set_rate =  35.89240351152864 + -20.283684291835424 * tfactors.T9i + -22.1209 * tfactors.T913i + 2.87067 * tfactors.T913
                             + -0.0200651 * tfactors.T9 + -0.0194061 * tfactors.T953 + -5.01065 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  20.283684291835424 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -22.1209 * tfactors.T943i + (1.0/3.0) * 2.87067 * tfactors.T923i
                                      + -0.0200651 + (5.0/3.0) * -0.0194061 * tfactors.T923 + -5.01065 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from nfisn
        ln_set_rate =  40.33000351152864 + -18.718495291835424 * tfactors.T9i + -41.059 * tfactors.T913i + 12.0361 * tfactors.T913
                             + -0.649888 * tfactors.T9 + 0.0175707 * tfactors.T953 + -4.96732 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  18.718495291835424 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -41.059 * tfactors.T943i + (1.0/3.0) * 12.0361 * tfactors.T923i
                                      + -0.649888 + (5.0/3.0) * 0.0175707 * tfactors.T923 + -4.96732 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from nfisn
        ln_set_rate =  43.39120351152864 + -20.763494291835425 * tfactors.T9i + -155.515 * tfactors.T913i + 135.903 * tfactors.T913
                             + -3.34812 * tfactors.T9 + -0.0139974 * tfactors.T953 + -85.9481 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  20.763494291835425 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -155.515 * tfactors.T943i + (1.0/3.0) * 135.903 * tfactors.T923i
                                      + -3.34812 + (5.0/3.0) * -0.0139974 * tfactors.T923 + -85.9481 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from nfisn
        ln_set_rate =  -96.49529648847135 + -53.72779429183542 * tfactors.T9i + 84.8583 * tfactors.T913i + 64.855 * tfactors.T913
                             + -6.59381 * tfactors.T9 + 0.413146 * tfactors.T953 + 0.07498000000000005 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  53.72779429183542 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 84.8583 * tfactors.T943i + (1.0/3.0) * 64.855 * tfactors.T923i
                                      + -6.59381 + (5.0/3.0) * 0.413146 * tfactors.T923 + 0.07498000000000005 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from nfisn
        ln_set_rate =  24.537293511528638 + -21.334824291835425 * tfactors.T9i + -27.0951 * tfactors.T913i + 22.7374 * tfactors.T913
                             + -2.18213 * tfactors.T9 + 0.223757 * tfactors.T953 + -13.8241 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  21.334824291835425 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -27.0951 * tfactors.T943i + (1.0/3.0) * 22.7374 * tfactors.T923i
                                      + -2.18213 + (5.0/3.0) * 0.223757 * tfactors.T923 + -13.8241 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_V45_to_p_Ti44_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V45_to_p_Ti44_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V45_to_He4_Sc41_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // V45 --> He4 + Sc41

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real V45_log_pf, dV45_log_pf_dT9;

        // interpolating V45 partition function
        get_partition_function_cached(V45, tfactors.T9, pf_cache, V45_log_pf, dV45_log_pf_dT9);

        amrex::Real Sc41_log_pf, dSc41_log_pf_dT9;

        // interpolating Sc41 partition function
        get_partition_function_cached(Sc41, tfactors.T9, pf_cache, Sc41_log_pf, dSc41_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Sc41_log_pf - V45_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dSc41_log_pf_dT9 - dV45_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  78.20656889296623 + -65.785501123702 * tfactors.T9i + -79.013 * tfactors.T913i + -6.64623 * tfactors.T913
                             + 1.10793 * tfactors.T9 + -0.169427 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  65.785501123702 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -79.013 * tfactors.T943i + (1.0/3.0) * -6.64623 * tfactors.T923i
                                      + 1.10793 + (5.0/3.0) * -0.169427 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_V45_to_He4_Sc41_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V45_to_He4_Sc41_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V47_to_p_Ti46_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // V47 --> p + Ti46

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real Ti46_log_pf, dTi46_log_pf_dT9;

        // interpolating Ti46 partition function
        get_partition_function_cached(Ti46, tfactors.T9, pf_cache, Ti46_log_pf, dTi46_log_pf_dT9);

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real V47_log_pf, dV47_log_pf_dT9;

        // interpolating V47 partition function
        get_partition_function_cached(V47, tfactors.T9, pf_cache, V47_log_pf, dV47_log_pf_dT9);

        amrex::Real net_log_pf = p_log_pf + Ti46_log_pf - V47_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dTi46_log_pf_dT9 - dV47_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  58.07481350908684 + -59.96949296023878 * tfactors.T9i + -33.2034 * tfactors.T913i + 0.505619 * tfactors.T913
                             + -1.73654 * tfactors.T9 + 0.207342 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  59.96949296023878 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -33.2034 * tfactors.T943i + (1.0/3.0) * 0.505619 * tfactors.T923i
                                      + -1.73654 + (5.0/3.0) * 0.207342 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_V47_to_p_Ti46_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V47_to_p_Ti46_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V47_to_He4_Sc43_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // V47 --> He4 + Sc43

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  95.66632297040007 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -79.122 * tfactors.T943i + (1.0/3.0) * -7.07006 * tfactors.T923i
                                      + 0.424183 + (5.0/3.0) * -0.0665231 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_V47_to_He4_Sc43_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V47_to_He4_Sc43_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr46_to_p_V45_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Cr46 --> p + V45

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real V45_log_pf, dV45_log_pf_dT9;

        // interpolating V45 partition function
        get_partition_function_cached(V45, tfactors.T9, pf_cache, V45_log_pf, dV45_log_pf_dT9);

        amrex::Real Cr46_log_pf, dCr46_log_pf_dT9;

        // interpolating Cr46 partition function
        get_partition_function_cached(Cr46, tfactors.T9, pf_cache, Cr46_log_pf, dCr46_log_pf_dT9);

        amrex::Real net_log_pf = p_log_pf + V45_log_pf - Cr46_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dV45_log_pf_dT9 - dCr46_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  61.69853410036488 + -56.56704824688144 * tfactors.T9i + -34.1967 * tfactors.T913i + -1.70253 * tfactors.T913
                             + -2.38375 * tfactors.T9 + 0.33421 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  56.56704824688144 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -34.1967 * tfactors.T943i + (1.0/3.0) * -1.70253 * tfactors.T923i
                                      + -2.38375 + (5.0/3.0) * 0.33421 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Cr46_to_p_V45_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr46_to_p_V45_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr46_to_He4_Ti42_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Cr46 --> He4 + Ti42

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Ti42_log_pf, dTi42_log_pf_dT9;

        // interpolating Ti42 partition function
        get_partition_function_cached(Ti42, tfactors.T9, pf_cache, Ti42_log_pf, dTi42_log_pf_dT9);

        amrex::Real Cr46_log_pf, dCr46_log_pf_dT9;

        // interpolating Cr46 partition function
        get_partition_function_cached(Cr46, tfactors.T9, pf_cache, Cr46_log_pf, dCr46_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Ti42_log_pf - Cr46_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dTi42_log_pf_dT9 - dCr46_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  79.53613543179074 + -78.82445373014531 * tfactors.T9i + -81.5598 * tfactors.T913i + -7.76887 * tfactors.T913
                             + 1.52445 * tfactors.T9 + -0.287125 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  78.82445373014531 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -81.5598 * tfactors.T943i + (1.0/3.0) * -7.76887 * tfactors.T923i
                                      + 1.52445 + (5.0/3.0) * -0.287125 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Cr46_to_He4_Ti42_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr46_to_He4_Ti42_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr48_to_p_V47_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Cr48 --> p + V47

        amrex::Real rate{}, drate_dT{};
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Cr48_to_p_V47_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr48_to_p_V47_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr48_to_He4_Ti44_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Cr48 --> He4 + Ti44

        amrex::Real rate{}, drate_dT{};
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Cr48_to_He4_Ti44_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr48_to_He4_Ti44_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr50_to_p_V49_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Cr50 --> p + V49

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real V49_log_pf, dV49_log_pf_dT9;

        // interpolating V49 partition function
        get_partition_function_cached(V49, tfactors.T9, pf_cache, V49_log_pf, dV49_log_pf_dT9);

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real Cr50_log_pf, dCr50_log_pf_dT9;

        // interpolating Cr50 partition function
        get_partition_function_cached(Cr50, tfactors.T9, pf_cache, Cr50_log_pf, dCr50_log_pf_dT9);

        amrex::Real net_log_pf = p_log_pf + V49_log_pf - Cr50_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dV49_log_pf_dT9 - dCr50_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  61.00795168922196 + -111.26564068398237 * tfactors.T9i + -34.217 * tfactors.T913i + 1.28258 * tfactors.T913
                             + -1.5098 * tfactors.T9 + 0.142011 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  111.26564068398237 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -34.217 * tfactors.T943i + (1.0/3.0) * 1.28258 * tfactors.T923i
                                      + -1.5098 + (5.0/3.0) * 0.142011 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Cr50_to_p_V49_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr50_to_p_V49_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr50_to_He4_Ti46_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Cr50 --> He4 + Ti46

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real Ti46_log_pf, dTi46_log_pf_dT9;

        // interpolating Ti46 partition function
        get_partition_function_cached(Ti46, tfactors.T9, pf_cache, Ti46_log_pf, dTi46_log_pf_dT9);

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Cr50_log_pf, dCr50_log_pf_dT9;

        // interpolating Cr50 partition function
        get_partition_function_cached(Cr50, tfactors.T9, pf_cache, Cr50_log_pf, dCr50_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Ti46_log_pf - Cr50_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dTi46_log_pf_dT9 - dCr50_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  82.00592864021027 + -99.31153420274927 * tfactors.T9i + -81.7658 * tfactors.T913i + 3.48637 * tfactors.T913
                             + -3.9188 * tfactors.T9 + 0.440356 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  99.31153420274927 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -81.7658 * tfactors.T943i + (1.0/3.0) * 3.48637 * tfactors.T923i
                                      + -3.9188 + (5.0/3.0) * 0.440356 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Cr50_to_He4_Ti46_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr50_to_He4_Ti46_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn49_to_p_Cr48_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Mn49 --> p + Cr48

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real Cr48_log_pf, dCr48_log_pf_dT9;

        // interpolating Cr48 partition function
        get_partition_function_cached(Cr48, tfactors.T9, pf_cache, Cr48_log_pf, dCr48_log_pf_dT9);

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real Mn49_log_pf, dMn49_log_pf_dT9;

        // interpolating Mn49 partition function
        get_partition_function_cached(Mn49, tfactors.T9, pf_cache, Mn49_log_pf, dMn49_log_pf_dT9);

        amrex::Real net_log_pf = p_log_pf + Cr48_log_pf - Mn49_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dCr48_log_pf_dT9 - dMn49_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  60.43582082218232 + -24.22873759767555 * tfactors.T9i + -35.1969 * tfactors.T913i + -3.36393 * tfactors.T913
                             + -1.55965 * tfactors.T9 + 0.237119 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  24.22873759767555 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -35.1969 * tfactors.T943i + (1.0/3.0) * -3.36393 * tfactors.T923i
                                      + -1.55965 + (5.0/3.0) * 0.237119 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Mn49_to_p_Cr48_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn49_to_p_Cr48_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn49_to_He4_V45_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Mn49 --> He4 + V45

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real V45_log_pf, dV45_log_pf_dT9;

        // interpolating V45 partition function
        get_partition_function_cached(V45, tfactors.T9, pf_cache, V45_log_pf, dV45_log_pf_dT9);

        amrex::Real Mn49_log_pf, dMn49_log_pf_dT9;

        // interpolating Mn49 partition function
        get_partition_function_cached(Mn49, tfactors.T9, pf_cache, Mn49_log_pf, dMn49_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + V45_log_pf - Mn49_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dV45_log_pf_dT9 - dMn49_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  77.07158196897863 + -94.67448480649189 * tfactors.T9i + -84.1758 * tfactors.T913i + -2.09835 * tfactors.T913
                             + -0.662318 * tfactors.T9 + 0.0293536 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  94.67448480649189 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -84.1758 * tfactors.T943i + (1.0/3.0) * -2.09835 * tfactors.T923i
                                      + -0.662318 + (5.0/3.0) * 0.0293536 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Mn49_to_He4_V45_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn49_to_He4_V45_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn51_to_p_Cr50_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Mn51 --> p + Cr50

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real Mn51_log_pf, dMn51_log_pf_dT9;

        // interpolating Mn51 partition function
        get_partition_function_cached(Mn51, tfactors.T9, pf_cache, Mn51_log_pf, dMn51_log_pf_dT9);

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real Cr50_log_pf, dCr50_log_pf_dT9;

        // interpolating Cr50 partition function
        get_partition_function_cached(Cr50, tfactors.T9, pf_cache, Cr50_log_pf, dCr50_log_pf_dT9);

        amrex::Real net_log_pf = p_log_pf + Cr50_log_pf - Mn51_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dCr50_log_pf_dT9 - dMn51_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  57.756820578960856 + -61.16522250740343 * tfactors.T9i + -35.2065 * tfactors.T913i + 0.884621 * tfactors.T913
                             + -1.76421 * tfactors.T9 + 0.202494 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  61.16522250740343 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -35.2065 * tfactors.T943i + (1.0/3.0) * 0.884621 * tfactors.T923i
                                      + -1.76421 + (5.0/3.0) * 0.202494 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Mn51_to_p_Cr50_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn51_to_p_Cr50_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn51_to_He4_V47_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Mn51 --> He4 + V47

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  100.50726374990671 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -84.2732 * tfactors.T943i + (1.0/3.0) * -2.98061 * tfactors.T923i
                                      + -0.531361 + (5.0/3.0) * 0.023612 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Mn51_to_He4_V47_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn51_to_He4_V47_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn53_to_p_Cr52_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Mn53 --> p + Cr52

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real Mn53_log_pf, dMn53_log_pf_dT9;

        // interpolating Mn53 partition function
        get_partition_function_cached(Mn53, tfactors.T9, pf_cache, Mn53_log_pf, dMn53_log_pf_dT9);

        amrex::Real Cr52_log_pf, dCr52_log_pf_dT9;

        // interpolating Cr52 partition function
        get_partition_function_cached(Cr52, tfactors.T9, pf_cache, Cr52_log_pf, dCr52_log_pf_dT9);

        amrex::Real net_log_pf = p_log_pf + Cr52_log_pf - Mn53_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dCr52_log_pf_dT9 - dMn53_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  57.31253842897789 + -76.122866140133 * tfactors.T9i + -35.2154 * tfactors.T913i + 0.73019 * tfactors.T913
                             + -1.26018 * tfactors.T9 + 0.121141 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  76.122866140133 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -35.2154 * tfactors.T943i + (1.0/3.0) * 0.73019 * tfactors.T923i
                                      + -1.26018 + (5.0/3.0) * 0.121141 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Mn53_to_p_Cr52_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn53_to_p_Cr52_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn53_to_He4_V49_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Mn53 --> He4 + V49

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real V49_log_pf, dV49_log_pf_dT9;

        // interpolating V49 partition function
        get_partition_function_cached(V49, tfactors.T9, pf_cache, V49_log_pf, dV49_log_pf_dT9);

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Mn53_log_pf, dMn53_log_pf_dT9;

        // interpolating Mn53 partition function
        get_partition_function_cached(Mn53, tfactors.T9, pf_cache, Mn53_log_pf, dMn53_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + V49_log_pf - Mn53_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dV49_log_pf_dT9 - dMn53_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  87.43161285745258 + -106.21633853020451 * tfactors.T9i + -84.3633 * tfactors.T913i + -5.28933 * tfactors.T913
                             + -1.3103 * tfactors.T9 + 0.187088 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  106.21633853020451 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -84.3633 * tfactors.T943i + (1.0/3.0) * -5.28933 * tfactors.T923i
                                      + -1.3103 + (5.0/3.0) * 0.187088 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Mn53_to_He4_V49_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn53_to_He4_V49_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe50_to_p_Mn49_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Fe50 --> p + Mn49

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real Fe50_log_pf, dFe50_log_pf_dT9;

        // interpolating Fe50 partition function
        get_partition_function_cached(Fe50, tfactors.T9, pf_cache, Fe50_log_pf, dFe50_log_pf_dT9);

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real Mn49_log_pf, dMn49_log_pf_dT9;

        // interpolating Mn49 partition function
        get_partition_function_cached(Mn49, tfactors.T9, pf_cache, Mn49_log_pf, dMn49_log_pf_dT9);

        amrex::Real net_log_pf = p_log_pf + Mn49_log_pf - Fe50_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dMn49_log_pf_dT9 - dFe50_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  61.5702885380999 + -48.10155227730567 * tfactors.T9i + -36.173 * tfactors.T913i + 1.28192 * tfactors.T913
                             + -3.04331 * tfactors.T9 + 0.393688 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  48.10155227730567 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -36.173 * tfactors.T943i + (1.0/3.0) * 1.28192 * tfactors.T923i
                                      + -3.04331 + (5.0/3.0) * 0.393688 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Fe50_to_p_Mn49_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe50_to_p_Mn49_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe50_to_He4_Cr46_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Fe50 --> He4 + Cr46

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real Fe50_log_pf, dFe50_log_pf_dT9;

        // interpolating Fe50 partition function
        get_partition_function_cached(Fe50, tfactors.T9, pf_cache, Fe50_log_pf, dFe50_log_pf_dT9);

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Cr46_log_pf, dCr46_log_pf_dT9;

        // interpolating Cr46 partition function
        get_partition_function_cached(Cr46, tfactors.T9, pf_cache, Cr46_log_pf, dCr46_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Cr46_log_pf - Fe50_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCr46_log_pf_dT9 - dFe50_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  81.32553640671367 + -86.2089888369161 * tfactors.T9i + -86.6499 * tfactors.T913i + -6.44312 * tfactors.T913
                             + 1.15026 * tfactors.T9 + -0.173961 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  86.2089888369161 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -86.6499 * tfactors.T943i + (1.0/3.0) * -6.44312 * tfactors.T923i
                                      + 1.15026 + (5.0/3.0) * -0.173961 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Fe50_to_He4_Cr46_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe50_to_He4_Cr46_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe52_to_p_Mn51_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Fe52 --> p + Mn51

        amrex::Real rate{}, drate_dT{};
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Fe52_to_p_Mn51_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe52_to_p_Mn51_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe52_to_He4_Cr48_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Fe52 --> He4 + Cr48

        amrex::Real rate{}, drate_dT{};
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Fe52_to_He4_Cr48_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe52_to_He4_Cr48_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe54_to_p_Mn53_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Fe54 --> p + Mn53

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real Mn53_log_pf, dMn53_log_pf_dT9;

        // interpolating Mn53 partition function
        get_partition_function_cached(Mn53, tfactors.T9, pf_cache, Mn53_log_pf, dMn53_log_pf_dT9);

        amrex::Real Fe54_log_pf, dFe54_log_pf_dT9;

        // interpolating Fe54 partition function
        get_partition_function_cached(Fe54, tfactors.T9, pf_cache, Fe54_log_pf, dFe54_log_pf_dT9);

        amrex::Real net_log_pf = p_log_pf + Mn53_log_pf - Fe54_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dMn53_log_pf_dT9 - dFe54_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  61.36267309371828 + -102.737944497193 * tfactors.T9i + -36.1913 * tfactors.T913i + 1.44056 * tfactors.T913
                             + -1.50734 * tfactors.T9 + 0.142124 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  102.737944497193 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -36.1913 * tfactors.T943i + (1.0/3.0) * 1.44056 * tfactors.T923i
                                      + -1.50734 + (5.0/3.0) * 0.142124 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Fe54_to_p_Mn53_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe54_to_p_Mn53_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe54_to_He4_Cr50_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Fe54 --> He4 + Cr50

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Cr50_log_pf, dCr50_log_pf_dT9;

        // interpolating Cr50 partition function
        get_partition_function_cached(Cr50, tfactors.T9, pf_cache, Cr50_log_pf, dCr50_log_pf_dT9);

        amrex::Real Fe54_log_pf, dFe54_log_pf_dT9;

        // interpolating Fe54 partition function
        get_partition_function_cached(Fe54, tfactors.T9, pf_cache, Fe54_log_pf, dFe54_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Cr50_log_pf - Fe54_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCr50_log_pf_dT9 - dFe54_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  84.69763426194889 + -97.68864234340795 * tfactors.T9i + -86.8349 * tfactors.T913i + 2.37343 * tfactors.T913
                             + -3.65677 * tfactors.T9 + 0.412128 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  97.68864234340795 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -86.8349 * tfactors.T943i + (1.0/3.0) * 2.37343 * tfactors.T923i
                                      + -3.65677 + (5.0/3.0) * 0.412128 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Fe54_to_He4_Cr50_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe54_to_He4_Cr50_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe56_to_He4_Cr52_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Fe56 --> He4 + Cr52

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Fe56_log_pf, dFe56_log_pf_dT9;

        // interpolating Fe56 partition function
        get_partition_function_cached(Fe56, tfactors.T9, pf_cache, Fe56_log_pf, dFe56_log_pf_dT9);

        amrex::Real Cr52_log_pf, dCr52_log_pf_dT9;

        // interpolating Cr52 partition function
        get_partition_function_cached(Cr52, tfactors.T9, pf_cache, Cr52_log_pf, dCr52_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Cr52_log_pf - Fe56_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCr52_log_pf_dT9 - dFe56_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  82.99716009465217 + -88.34015858985633 * tfactors.T9i + -86.9175 * tfactors.T913i + 5.655 * tfactors.T913
                             + -4.53625 * tfactors.T9 + 0.479806 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  88.34015858985633 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -86.9175 * tfactors.T943i + (1.0/3.0) * 5.655 * tfactors.T923i
                                      + -4.53625 + (5.0/3.0) * 0.479806 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Fe56_to_He4_Cr52_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe56_to_He4_Cr52_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co53_to_p_Fe52_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Co53 --> p + Fe52

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real Co53_log_pf, dCo53_log_pf_dT9;

        // interpolating Co53 partition function
        get_partition_function_cached(Co53, tfactors.T9, pf_cache, Co53_log_pf, dCo53_log_pf_dT9);

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real Fe52_log_pf, dFe52_log_pf_dT9;

        // interpolating Fe52 partition function
        get_partition_function_cached(Fe52, tfactors.T9, pf_cache, Fe52_log_pf, dFe52_log_pf_dT9);

        amrex::Real net_log_pf = p_log_pf + Fe52_log_pf - Co53_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dFe52_log_pf_dT9 - dCo53_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  60.90479478903703 + -18.756046851504873 * tfactors.T9i + -37.1457 * tfactors.T913i + -5.67417 * tfactors.T913
                             + -0.559644 * tfactors.T9 + 0.102768 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  18.756046851504873 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -37.1457 * tfactors.T943i + (1.0/3.0) * -5.67417 * tfactors.T923i
                                      + -0.559644 + (5.0/3.0) * 0.102768 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Co53_to_p_Fe52_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co53_to_p_Fe52_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co53_to_He4_Mn49_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Co53 --> He4 + Mn49

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real Co53_log_pf, dCo53_log_pf_dT9;

        // interpolating Co53 partition function
        get_partition_function_cached(Co53, tfactors.T9, pf_cache, Co53_log_pf, dCo53_log_pf_dT9);

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Mn49_log_pf, dMn49_log_pf_dT9;

        // interpolating Mn49 partition function
        get_partition_function_cached(Mn49, tfactors.T9, pf_cache, Mn49_log_pf, dMn49_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Mn49_log_pf - Co53_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dMn49_log_pf_dT9 - dCo53_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  78.4777879593697 + -86.62094923020983 * tfactors.T9i + -89.1861 * tfactors.T913i + -0.631651 * tfactors.T913
                             + -0.826494 * tfactors.T9 + 0.0355619 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  86.62094923020983 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -89.1861 * tfactors.T943i + (1.0/3.0) * -0.631651 * tfactors.T923i
                                      + -0.826494 + (5.0/3.0) * 0.0355619 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Co53_to_He4_Mn49_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co53_to_He4_Mn49_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co55_to_p_Fe54_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Co55 --> p + Fe54

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  58.76958578644693 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -37.1544 * tfactors.T943i + (1.0/3.0) * 0.950364 * tfactors.T923i
                                      + -1.77529 + (5.0/3.0) * 0.198562 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Co55_to_p_Fe54_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co55_to_p_Fe54_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co55_to_He4_Mn51_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Co55 --> He4 + Mn51

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  95.29300562248025 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -89.274 * tfactors.T943i + (1.0/3.0) * -10.4373 * tfactors.T923i
                                      + 1.00492 + (5.0/3.0) * -0.125548 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Co55_to_He4_Mn51_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co55_to_He4_Mn51_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co56_to_n_Co55_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Co56 --> n + Co55

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  116.994640839646 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.86357 * tfactors.T923i
                                      + 0.616591 + (5.0/3.0) * -0.0839313 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Co56_to_n_Co55_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co56_to_n_Co55_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co57_to_n_Co56_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Co57 --> n + Co56

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  132.01901045161702 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.37855 * tfactors.T923i
                                      + 0.299896 + (5.0/3.0) * -0.04382 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Co57_to_n_Co56_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co57_to_n_Co56_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co57_to_p_Fe56_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Co57 --> p + Fe56

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  69.94636137001572 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -37.1625 * tfactors.T943i + (1.0/3.0) * 1.06776 * tfactors.T923i
                                      + -1.31689 + (5.0/3.0) * 0.122089 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Co57_to_p_Fe56_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co57_to_p_Fe56_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co57_to_He4_Mn53_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Co57 --> He4 + Mn53

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Mn53_log_pf, dMn53_log_pf_dT9;

        // interpolating Mn53 partition function
        get_partition_function_cached(Mn53, tfactors.T9, pf_cache, Mn53_log_pf, dMn53_log_pf_dT9);

        amrex::Real Co57_log_pf, dCo57_log_pf_dT9;

        // interpolating Co57 partition function
        get_partition_function_cached(Co57, tfactors.T9, pf_cache, Co57_log_pf, dCo57_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Mn53_log_pf - Co57_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dMn53_log_pf_dT9 - dCo57_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  87.00976353805295 + -82.16365381973904 * tfactors.T9i + -89.3558 * tfactors.T913i + 0.000882861 * tfactors.T913
                             + -2.79327 * tfactors.T9 + 0.309057 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  82.16365381973904 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -89.3558 * tfactors.T943i + (1.0/3.0) * 0.000882861 * tfactors.T923i
                                      + -2.79327 + (5.0/3.0) * 0.309057 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Co57_to_He4_Mn53_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co57_to_He4_Mn53_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni56_to_p_Co55_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ni56 --> p + Co55

        amrex::Real rate{}, drate_dT{};
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ni56_to_p_Co55_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni56_to_p_Co55_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni56_to_He4_Fe52_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ni56 --> He4 + Fe52

        amrex::Real rate{}, drate_dT{};
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ni56_to_He4_Fe52_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni56_to_He4_Fe52_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni58_to_p_Co57_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ni58 --> p + Co57

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  94.83410720454592 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -38.1133 * tfactors.T943i + (1.0/3.0) * 1.77414 * tfactors.T923i
                                      + -1.48268 + (5.0/3.0) * 0.121073 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ni58_to_p_Co57_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni58_to_p_Co57_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni58_to_He4_Fe54_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ni58 --> He4 + Fe54

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  74.25981652709197 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -91.7628 * tfactors.T943i + (1.0/3.0) * 4.23027 * tfactors.T923i
                                      + -3.31305 + (5.0/3.0) * 0.271293 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ni58_to_He4_Fe54_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni58_to_He4_Fe54_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu59_to_p_Ni58_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Cu59 --> p + Ni58

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  39.67086986195644 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -39.0526 * tfactors.T943i + (1.0/3.0) * 1.0436 * tfactors.T923i
                                      + -2.10834 + (5.0/3.0) * 0.239559 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Cu59_to_p_Ni58_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu59_to_p_Ni58_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu59_to_He4_Co55_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Cu59 --> He4 + Co55

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  55.1611006026015 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -94.1404 * tfactors.T943i + (1.0/3.0) * -2.62786 * tfactors.T923i
                                      + -2.12066 + (5.0/3.0) * 0.237999 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Cu59_to_He4_Co55_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu59_to_He4_Co55_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Zn60_to_p_Cu59_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Zn60 --> p + Cu59

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  59.240729222177784 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -39.9808 * tfactors.T943i + (1.0/3.0) * 1.17311 * tfactors.T923i
                                      + -2.90486 + (5.0/3.0) * 0.339644 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Zn60_to_p_Cu59_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn60_to_p_Cu59_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Zn60_to_He4_Ni56_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Zn60 --> He4 + Ni56

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  31.2372260432927 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -96.4898 * tfactors.T943i + (1.0/3.0) * 6.47209 * tfactors.T923i
                                      + -5.2029 + (5.0/3.0) * 0.533391 * tfactors.T923 + 0.833333 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Zn60_to_He4_Ni56_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn60_to_He4_Ni56_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_C12_to_3He4_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // C12 --> 3 He4

        amrex::Real rate{}, drate_dT{};
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_C12_to_3He4_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_C12_to_3He4_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_C12_C12_to_n_Mg23_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // C12 + C12 --> n + Mg23

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Mg23_log_pf, dMg23_log_pf_dT9;

        // interpolating Mg23 partition function
        get_partition_function_cached(Mg23, tfactors.T9, pf_cache, Mg23_log_pf, dMg23_log_pf_dT9);

        amrex::Real C12_log_pf, dC12_log_pf_dT9;

        // setting C12 log(partition function) to 0.0 by default, independent of T
        C12_log_pf = 0.0_rt;
        dC12_log_pf_dT9 = 0.0_rt;

        amrex::Real net_log_pf = n_log_pf + Mg23_log_pf - C12_log_pf - C12_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dMg23_log_pf_dT9 - dC12_log_pf_dT9 - dC12_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from cf88r
        ln_set_rate =  -12.784659727528789 + -30.144454449865492 * tfactors.T9i + 11.4826 * tfactors.T913
                             + 1.82849 * tfactors.T9 + -0.34844 * tfactors.T953;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  30.144454449865492 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 11.4826 * tfactors.T923i
                                      + 1.82849 + (5.0/3.0) * -0.34844 * tfactors.T923;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_C12_C12_to_n_Mg23_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_C12_C12_to_n_Mg23_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_O16_C12_to_n_Si27_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // O16 + C12 --> n + Si27

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real O16_log_pf, dO16_log_pf_dT9;

        // interpolating O16 partition function
        get_partition_function_cached(O16, tfactors.T9, pf_cache, O16_log_pf, dO16_log_pf_dT9);

        amrex::Real Si27_log_pf, dSi27_log_pf_dT9;

        // interpolating Si27 partition function
        get_partition_function_cached(Si27, tfactors.T9, pf_cache, Si27_log_pf, dSi27_log_pf_dT9);

        amrex::Real C12_log_pf, dC12_log_pf_dT9;

        // setting C12 log(partition function) to 0.0 by default, independent of T
        C12_log_pf = 0.0_rt;
        dC12_log_pf_dT9 = 0.0_rt;

        amrex::Real net_log_pf = n_log_pf + Si27_log_pf - C12_log_pf - O16_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dSi27_log_pf_dT9 - dC12_log_pf_dT9 - dO16_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from cf88r
        ln_set_rate =  -132.19144876580245 + -1.4742031459903027 * tfactors.T9i + -293.089 * tfactors.T913i + 414.404 * tfactors.T913
                             + -28.0562 * tfactors.T9 + 1.61807 * tfactors.T953 + -178.28 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  1.4742031459903027 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -293.089 * tfactors.T943i + (1.0/3.0) * 414.404 * tfactors.T923i
                                      + -28.0562 + (5.0/3.0) * 1.61807 * tfactors.T923 + -178.28 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_O16_C12_to_n_Si27_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_O16_C12_to_n_Si27_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ne20_He4_to_n_Mg23_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ne20 + He4 --> n + Mg23

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Mg23_log_pf, dMg23_log_pf_dT9;

        // interpolating Mg23 partition function
        get_partition_function_cached(Mg23, tfactors.T9, pf_cache, Mg23_log_pf, dMg23_log_pf_dT9);

        amrex::Real Ne20_log_pf, dNe20_log_pf_dT9;

        // interpolating Ne20 partition function
        get_partition_function_cached(Ne20, tfactors.T9, pf_cache, Ne20_log_pf, dNe20_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Mg23_log_pf - He4_log_pf - Ne20_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dMg23_log_pf_dT9 - dHe4_log_pf_dT9 - dNe20_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  17.974691741390394 + -83.72270411882882 * tfactors.T9i + 1.83199 * tfactors.T913
                             + -0.290485 * tfactors.T9 + 0.0242929 * tfactors.T953;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  83.72270411882882 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 1.83199 * tfactors.T923i
                                      + -0.290485 + (5.0/3.0) * 0.0242929 * tfactors.T923;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ne20_He4_to_n_Mg23_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne20_He4_to_n_Mg23_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mg24_He4_to_n_Si27_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Mg24 + He4 --> n + Si27

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Mg24_log_pf, dMg24_log_pf_dT9;

        // interpolating Mg24 partition function
        get_partition_function_cached(Mg24, tfactors.T9, pf_cache, Mg24_log_pf, dMg24_log_pf_dT9);

        amrex::Real Si27_log_pf, dSi27_log_pf_dT9;

        // interpolating Si27 partition function
        get_partition_function_cached(Si27, tfactors.T9, pf_cache, Si27_log_pf, dSi27_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Si27_log_pf - He4_log_pf - Mg24_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dSi27_log_pf_dT9 - dHe4_log_pf_dT9 - dMg24_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  18.680909426508013 + -83.49996467707591 * tfactors.T9i + 0.942981 * tfactors.T913
                             + -0.104624 * tfactors.T9 + 0.00723421 * tfactors.T953;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  83.49996467707591 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 0.942981 * tfactors.T923i
                                      + -0.104624 + (5.0/3.0) * 0.00723421 * tfactors.T923;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Mg24_He4_to_n_Si27_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg24_He4_to_n_Si27_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mg24_He4_to_p_Al27_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Mg24 + He4 --> p + Al27

        amrex::Real rate{}, drate_dT{};
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
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  22.818701538287602 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 18.0416 * tfactors.T923i
                                      + -1.54137 + (5.0/3.0) * 0.0847506 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from il10r
        ln_set_rate =  -26.268451753058486 + -19.5392135382876 * tfactors.T9i + 5.18642 * tfactors.T913i + -34.7936 * tfactors.T913
                             + 168.225 * tfactors.T9 + -115.825 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  19.5392135382876 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 5.18642 * tfactors.T943i + (1.0/3.0) * -34.7936 * tfactors.T923i
                                      + 168.225 + (5.0/3.0) * -115.825 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from il10n
        ln_set_rate =  30.057448246941515 + -18.576201538287602 * tfactors.T9i + -26.4162 * tfactors.T913i
                             + -2.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  18.576201538287602 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -26.4162 * tfactors.T943i
                                      + (5.0/3.0) * -2.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Mg24_He4_to_p_Al27_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg24_He4_to_p_Al27_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mg24_He4_to_C12_O16_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Mg24 + He4 --> C12 + O16

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  78.2098115311048 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -133.413 * tfactors.T943i + (1.0/3.0) * 50.1572 * tfactors.T923i
                                      + -3.15987 + (5.0/3.0) * 0.0178251 * tfactors.T923 + -23.7027 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Mg24_He4_to_C12_O16_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg24_He4_to_C12_O16_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Al27_p_to_n_Si27_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Al27 + p --> n + Si27

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real Al27_log_pf, dAl27_log_pf_dT9;

        // interpolating Al27 partition function
        get_partition_function_cached(Al27, tfactors.T9, pf_cache, Al27_log_pf, dAl27_log_pf_dT9);

        amrex::Real Si27_log_pf, dSi27_log_pf_dT9;

        // interpolating Si27 partition function
        get_partition_function_cached(Si27, tfactors.T9, pf_cache, Si27_log_pf, dSi27_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Si27_log_pf - p_log_pf - Al27_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dSi27_log_pf_dT9 - dp_log_pf_dT9 - dAl27_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  17.8425611795665 + -64.9237631387883 * tfactors.T9i + 2.10552 * tfactors.T913
                             + -0.318695 * tfactors.T9 + 0.0257496 * tfactors.T953;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  64.9237631387883 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 2.10552 * tfactors.T923i
                                      + -0.318695 + (5.0/3.0) * 0.0257496 * tfactors.T923;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Al27_p_to_n_Si27_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Al27_p_to_n_Si27_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Al27_p_to_C12_O16_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Al27 + p --> C12 + O16

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  59.8005159928172 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -119.242 * tfactors.T943i + (1.0/3.0) * 13.3667 * tfactors.T923i
                                      + 0.295425 + (5.0/3.0) * -0.267288 * tfactors.T923 + -9.91729 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Al27_p_to_C12_O16_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Al27_p_to_C12_O16_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si27_He4_to_n_S30_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Si27 + He4 --> n + S30

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real S30_log_pf, dS30_log_pf_dT9;

        // interpolating S30 partition function
        get_partition_function_cached(S30, tfactors.T9, pf_cache, S30_log_pf, dS30_log_pf_dT9);

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Si27_log_pf, dSi27_log_pf_dT9;

        // interpolating Si27 partition function
        get_partition_function_cached(Si27, tfactors.T9, pf_cache, Si27_log_pf, dSi27_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + S30_log_pf - He4_log_pf - Si27_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dS30_log_pf_dT9 - dHe4_log_pf_dT9 - dSi27_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  11.342202681948677 + -46.08922632071148 * tfactors.T9i + 1.55858 * tfactors.T913
                             + -0.0243677 * tfactors.T9 + 0.00835247 * tfactors.T953;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  46.08922632071148 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 1.55858 * tfactors.T923i
                                      + -0.0243677 + (5.0/3.0) * 0.00835247 * tfactors.T923;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Si27_He4_to_n_S30_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si27_He4_to_n_S30_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si28_He4_to_p_P31_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Si28 + He4 --> p + P31

        amrex::Real rate{}, drate_dT{};
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
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  25.663579682396225 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 21.521 * tfactors.T923i
                                      + -1.90355 + (5.0/3.0) * 0.092724 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from il10r
        ln_set_rate =  -13.441669089666712 + -24.114989682396224 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  24.114989682396224 * tfactors.T9i * tfactors.T9i
                                      + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from il10n
        ln_set_rate =  60.360230910333286 + -22.237829682396224 * tfactors.T9i + -31.932 * tfactors.T913i + -77.0334 * tfactors.T913
                             + -43.6847 * tfactors.T9 + -4.28955 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  22.237829682396224 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -31.932 * tfactors.T943i + (1.0/3.0) * -77.0334 * tfactors.T923i
                                      + -43.6847 + (5.0/3.0) * -4.28955 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Si28_He4_to_p_P31_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si28_He4_to_p_P31_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si28_He4_to_C12_Ne20_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Si28 + He4 --> C12 + Ne20

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  186.77284344745925 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 514.197 * tfactors.T943i + (1.0/3.0) * -200.896 * tfactors.T923i
                                      + -6.42713 + (5.0/3.0) * 0.758256 * tfactors.T923 + 236.359 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Si28_He4_to_C12_Ne20_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si28_He4_to_C12_Ne20_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si28_He4_to_O16_O16_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Si28 + He4 --> O16 + O16

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  111.60083385054465 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -119.324 * tfactors.T943i + (1.0/3.0) * -32.2497 * tfactors.T923i
                                      + 1.46214 + (5.0/3.0) * -0.200893 * tfactors.T923 + 13.2148 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Si28_He4_to_O16_O16_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si28_He4_to_O16_O16_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si30_p_to_He4_Al27_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Si30 + p --> He4 + Al27

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real Si30_log_pf, dSi30_log_pf_dT9;

        // interpolating Si30 partition function
        get_partition_function_cached(Si30, tfactors.T9, pf_cache, Si30_log_pf, dSi30_log_pf_dT9);

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Al27_log_pf, dAl27_log_pf_dT9;

        // interpolating Al27 partition function
        get_partition_function_cached(Al27, tfactors.T9, pf_cache, Al27_log_pf, dAl27_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Al27_log_pf - p_log_pf - Si30_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dAl27_log_pf_dT9 - dp_log_pf_dT9 - dSi30_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  50.48822991071769 + -27.52646014539429 * tfactors.T9i + -56.5351 * tfactors.T913i + 1.60477 * tfactors.T913
                             + -1.40594 * tfactors.T9 + 0.127353 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  27.52646014539429 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -56.5351 * tfactors.T943i + (1.0/3.0) * 1.60477 * tfactors.T923i
                                      + -1.40594 + (5.0/3.0) * 0.127353 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Si30_p_to_He4_Al27_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si30_p_to_He4_Al27_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P31_p_to_C12_Ne20_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // P31 + p --> C12 + Ne20

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  156.07991376506303 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 361.154 * tfactors.T943i + (1.0/3.0) * -92.643 * tfactors.T923i
                                      + -9.98738 + (5.0/3.0) * 0.892737 * tfactors.T923 + 161.042 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_P31_p_to_C12_Ne20_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P31_p_to_C12_Ne20_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P31_p_to_O16_O16_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // P31 + p --> O16 + O16

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  88.87103716814843 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -145.844 * tfactors.T943i + (1.0/3.0) * 8.72612 * tfactors.T923i
                                      + -0.554035 + (5.0/3.0) * -0.137562 * tfactors.T923 + -6.88807 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_P31_p_to_O16_O16_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P31_p_to_O16_O16_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S32_p_to_He4_P29_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // S32 + p --> He4 + P29

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real P29_log_pf, dP29_log_pf_dT9;

        // interpolating P29 partition function
        get_partition_function_cached(P29, tfactors.T9, pf_cache, P29_log_pf, dP29_log_pf_dT9);

        amrex::Real S32_log_pf, dS32_log_pf_dT9;

        // interpolating S32 partition function
        get_partition_function_cached(S32, tfactors.T9, pf_cache, S32_log_pf, dS32_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + P29_log_pf - p_log_pf - S32_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dP29_log_pf_dT9 - dp_log_pf_dT9 - dS32_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  52.19693496988301 + -48.723680264833554 * tfactors.T9i + -62.3802 * tfactors.T913i + 0.459085 * tfactors.T913
                             + -0.870169 * tfactors.T9 + 0.0631143 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  48.723680264833554 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -62.3802 * tfactors.T943i + (1.0/3.0) * 0.459085 * tfactors.T923i
                                      + -0.870169 + (5.0/3.0) * 0.0631143 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_S32_p_to_He4_P29_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S32_p_to_He4_P29_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S32_He4_to_p_Cl35_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // S32 + He4 --> p + Cl35

        amrex::Real rate{}, drate_dT{};
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
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  25.589703706656966 * tfactors.T9i * tfactors.T9i
                                      + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from il10r
        ln_set_rate =  -57.37716473446963 + -22.187684706656967 * tfactors.T9i + 25.5338 * tfactors.T913
                             + 6.45824 * tfactors.T9 + -0.950294 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  22.187684706656967 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 25.5338 * tfactors.T923i
                                      + 6.45824 + (5.0/3.0) * -0.950294 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from il10n
        ln_set_rate =  32.27223526553037 + -21.654753706656965 * tfactors.T9i + -30.9147 * tfactors.T913i + -1.2345 * tfactors.T913
                             + 22.5118 * tfactors.T9 + -33.0589 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  21.654753706656965 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -30.9147 * tfactors.T943i + (1.0/3.0) * -1.2345 * tfactors.T923i
                                      + 22.5118 + (5.0/3.0) * -33.0589 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from il10r
        ln_set_rate =  2.4434452655303747 + -27.664513706656965 * tfactors.T9i + 5.33756 * tfactors.T913
                             + 1.64418 * tfactors.T9 + -0.246167 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  27.664513706656965 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 5.33756 * tfactors.T923i
                                      + 1.64418 + (5.0/3.0) * -0.246167 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_S32_He4_to_p_Cl35_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S32_He4_to_p_Cl35_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S34_p_to_He4_P31_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // S34 + p --> He4 + P31

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real S34_log_pf, dS34_log_pf_dT9;

        // interpolating S34 partition function
        get_partition_function_cached(S34, tfactors.T9, pf_cache, S34_log_pf, dS34_log_pf_dT9);

        amrex::Real P31_log_pf, dP31_log_pf_dT9;

        // interpolating P31 partition function
        get_partition_function_cached(P31, tfactors.T9, pf_cache, P31_log_pf, dP31_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + P31_log_pf - p_log_pf - S34_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dP31_log_pf_dT9 - dp_log_pf_dT9 - dS34_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  51.02256860140231 + -7.277081980312247 * tfactors.T9i + -62.5433 * tfactors.T913i + 1.69339 * tfactors.T913
                             + -1.27872 * tfactors.T9 + 0.106616 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  7.277081980312247 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -62.5433 * tfactors.T943i + (1.0/3.0) * 1.69339 * tfactors.T923i
                                      + -1.27872 + (5.0/3.0) * 0.106616 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_S34_p_to_He4_P31_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S34_p_to_He4_P31_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl33_p_to_He4_S30_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Cl33 + p --> He4 + S30

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real S30_log_pf, dS30_log_pf_dT9;

        // interpolating S30 partition function
        get_partition_function_cached(S30, tfactors.T9, pf_cache, S30_log_pf, dS30_log_pf_dT9);

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Cl33_log_pf, dCl33_log_pf_dT9;

        // interpolating Cl33 partition function
        get_partition_function_cached(Cl33, tfactors.T9, pf_cache, Cl33_log_pf, dCl33_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + S30_log_pf - p_log_pf - Cl33_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dS30_log_pf_dT9 - dp_log_pf_dT9 - dCl33_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  50.59092547517039 + -24.13733741901471 * tfactors.T9i + -65.211 * tfactors.T913i + 0.403887 * tfactors.T913
                             + -0.808239 * tfactors.T9 + 0.0581461 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  24.13733741901471 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -65.211 * tfactors.T943i + (1.0/3.0) * 0.403887 * tfactors.T923i
                                      + -0.808239 + (5.0/3.0) * 0.0581461 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Cl33_p_to_He4_S30_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl33_p_to_He4_S30_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar36_p_to_He4_Cl33_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ar36 + p --> He4 + Cl33

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real Cl33_log_pf, dCl33_log_pf_dT9;

        // interpolating Cl33 partition function
        get_partition_function_cached(Cl33, tfactors.T9, pf_cache, Cl33_log_pf, dCl33_log_pf_dT9);

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Ar36_log_pf, dAr36_log_pf_dT9;

        // interpolating Ar36 partition function
        get_partition_function_cached(Ar36, tfactors.T9, pf_cache, Ar36_log_pf, dAr36_log_pf_dT9);

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real net_log_pf = He4_log_pf + Cl33_log_pf - p_log_pf - Ar36_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCl33_log_pf_dT9 - dp_log_pf_dT9 - dAr36_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  54.783557788881225 + -50.64428487608636 * tfactors.T9i + -68.1442 * tfactors.T913i + 0.291238 * tfactors.T913
                             + -0.791384 * tfactors.T9 + 0.0524823 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  50.64428487608636 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -68.1442 * tfactors.T943i + (1.0/3.0) * 0.291238 * tfactors.T923i
                                      + -0.791384 + (5.0/3.0) * 0.0524823 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ar36_p_to_He4_Cl33_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar36_p_to_He4_Cl33_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar36_He4_to_p_K39_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ar36 + He4 --> p + K39

        amrex::Real rate{}, drate_dT{};
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
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  14.951286608188834 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -30.0732 * tfactors.T943i + (1.0/3.0) * 7.03263 * tfactors.T923i
                                      + -1.10085 + (5.0/3.0) * 0.133768 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ar36_He4_to_p_K39_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar36_He4_to_p_K39_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar38_p_to_He4_Cl35_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ar38 + p --> He4 + Cl35

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real Ar38_log_pf, dAr38_log_pf_dT9;

        // interpolating Ar38 partition function
        get_partition_function_cached(Ar38, tfactors.T9, pf_cache, Ar38_log_pf, dAr38_log_pf_dT9);

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Cl35_log_pf, dCl35_log_pf_dT9;

        // interpolating Cl35 partition function
        get_partition_function_cached(Cl35, tfactors.T9, pf_cache, Cl35_log_pf, dCl35_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Cl35_log_pf - p_log_pf - Ar38_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCl35_log_pf_dT9 - dp_log_pf_dT9 - dAr38_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  53.758798882415604 + -9.715822523380488 * tfactors.T9i + -68.2848 * tfactors.T913i + 2.5993 * tfactors.T913
                             + -1.59144 * tfactors.T9 + 0.137745 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  9.715822523380488 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -68.2848 * tfactors.T943i + (1.0/3.0) * 2.5993 * tfactors.T923i
                                      + -1.59144 + (5.0/3.0) * 0.137745 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ar38_p_to_He4_Cl35_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar38_p_to_He4_Cl35_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K37_p_to_He4_Ar34_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // K37 + p --> He4 + Ar34

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Ar34_log_pf, dAr34_log_pf_dT9;

        // interpolating Ar34 partition function
        get_partition_function_cached(Ar34, tfactors.T9, pf_cache, Ar34_log_pf, dAr34_log_pf_dT9);

        amrex::Real K37_log_pf, dK37_log_pf_dT9;

        // interpolating K37 partition function
        get_partition_function_cached(K37, tfactors.T9, pf_cache, K37_log_pf, dK37_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Ar34_log_pf - p_log_pf - K37_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dAr34_log_pf_dT9 - dp_log_pf_dT9 - dK37_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  53.03935240918314 + -18.07815432695771 * tfactors.T9i + -70.8666 * tfactors.T913i + 0.263069 * tfactors.T913
                             + -0.802985 * tfactors.T9 + 0.0550287 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  18.07815432695771 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -70.8666 * tfactors.T943i + (1.0/3.0) * 0.263069 * tfactors.T923i
                                      + -0.802985 + (5.0/3.0) * 0.0550287 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_K37_p_to_He4_Ar34_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K37_p_to_He4_Ar34_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K39_He4_to_p_Ca42_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // K39 + He4 --> p + Ca42

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real K39_log_pf, dK39_log_pf_dT9;

        // interpolating K39 partition function
        get_partition_function_cached(K39, tfactors.T9, pf_cache, K39_log_pf, dK39_log_pf_dT9);

        amrex::Real Ca42_log_pf, dCa42_log_pf_dT9;

        // interpolating Ca42 partition function
        get_partition_function_cached(Ca42, tfactors.T9, pf_cache, Ca42_log_pf, dCa42_log_pf_dT9);

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real net_log_pf = p_log_pf + Ca42_log_pf - He4_log_pf - K39_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dCa42_log_pf_dT9 - dHe4_log_pf_dT9 - dK39_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  -80.85570368170251 + -1.4384983175969779 * tfactors.T9i + -31.1381 * tfactors.T913i + 117.868 * tfactors.T913
                             + -27.0455 * tfactors.T9 + 2.75661 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  1.4384983175969779 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -31.1381 * tfactors.T943i + (1.0/3.0) * 117.868 * tfactors.T923i
                                      + -27.0455 + (5.0/3.0) * 2.75661 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_K39_He4_to_p_Ca42_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K39_He4_to_p_Ca42_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca40_p_to_He4_K37_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ca40 + p --> He4 + K37

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real K37_log_pf, dK37_log_pf_dT9;

        // interpolating K37 partition function
        get_partition_function_cached(K37, tfactors.T9, pf_cache, K37_log_pf, dK37_log_pf_dT9);

        amrex::Real Ca40_log_pf, dCa40_log_pf_dT9;

        // interpolating Ca40 partition function
        get_partition_function_cached(Ca40, tfactors.T9, pf_cache, Ca40_log_pf, dCa40_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + K37_log_pf - p_log_pf - Ca40_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dK37_log_pf_dT9 - dp_log_pf_dT9 - dCa40_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  56.9349260689074 + -60.13631651881774 * tfactors.T9i + -73.6776 * tfactors.T913i + 0.493294 * tfactors.T913
                             + -0.902307 * tfactors.T9 + 0.0641391 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  60.13631651881774 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -73.6776 * tfactors.T943i + (1.0/3.0) * 0.493294 * tfactors.T923i
                                      + -0.902307 + (5.0/3.0) * 0.0641391 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ca40_p_to_He4_K37_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca40_p_to_He4_K37_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca40_He4_to_p_Sc43_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ca40 + He4 --> p + Sc43

        amrex::Real rate{}, drate_dT{};
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
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  40.87525788838128 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -32.1734 * tfactors.T943i + (1.0/3.0) * 0.0296879 * tfactors.T923i
                                      + -0.95232 + (5.0/3.0) * 0.129022 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ca40_He4_to_p_Sc43_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca40_He4_to_p_Sc43_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc41_p_to_He4_Ca38_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Sc41 + p --> He4 + Ca38

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Sc41_log_pf, dSc41_log_pf_dT9;

        // interpolating Sc41 partition function
        get_partition_function_cached(Sc41, tfactors.T9, pf_cache, Sc41_log_pf, dSc41_log_pf_dT9);

        amrex::Real Ca38_log_pf, dCa38_log_pf_dT9;

        // interpolating Ca38 partition function
        get_partition_function_cached(Ca38, tfactors.T9, pf_cache, Ca38_log_pf, dCa38_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Ca38_log_pf - p_log_pf - Sc41_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCa38_log_pf_dT9 - dp_log_pf_dT9 - dSc41_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  55.26173122707938 + -19.95750603673758 * tfactors.T9i + -76.3063 * tfactors.T913i + -0.741136 * tfactors.T913
                             + -0.511467 * tfactors.T9 + 0.0171615 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  19.95750603673758 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -76.3063 * tfactors.T943i + (1.0/3.0) * -0.741136 * tfactors.T923i
                                      + -0.511467 + (5.0/3.0) * 0.0171615 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Sc41_p_to_He4_Ca38_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc41_p_to_He4_Ca38_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti44_p_to_He4_Sc41_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ti44 + p --> He4 + Sc41

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Ti44_log_pf, dTi44_log_pf_dT9;

        // interpolating Ti44 partition function
        get_partition_function_cached(Ti44, tfactors.T9, pf_cache, Ti44_log_pf, dTi44_log_pf_dT9);

        amrex::Real Sc41_log_pf, dSc41_log_pf_dT9;

        // interpolating Sc41 partition function
        get_partition_function_cached(Sc41, tfactors.T9, pf_cache, Sc41_log_pf, dSc41_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Sc41_log_pf - p_log_pf - Ti44_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dSc41_log_pf_dT9 - dp_log_pf_dT9 - dTi44_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  60.827665381437605 + -46.90760683186659 * tfactors.T9i + -79.013 * tfactors.T913i + -0.989013 * tfactors.T913
                             + -0.463943 * tfactors.T9 + 0.00902542 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  46.90760683186659 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -79.013 * tfactors.T943i + (1.0/3.0) * -0.989013 * tfactors.T923i
                                      + -0.463943 + (5.0/3.0) * 0.00902542 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ti44_p_to_He4_Sc41_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti44_p_to_He4_Sc41_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti46_p_to_He4_Sc43_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ti46 + p --> He4 + Sc43

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Ti46_log_pf, dTi46_log_pf_dT9;

        // interpolating Ti46 partition function
        get_partition_function_cached(Ti46, tfactors.T9, pf_cache, Ti46_log_pf, dTi46_log_pf_dT9);

        amrex::Real Sc43_log_pf, dSc43_log_pf_dT9;

        // interpolating Sc43 partition function
        get_partition_function_cached(Sc43, tfactors.T9, pf_cache, Sc43_log_pf, dSc43_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Sc43_log_pf - p_log_pf - Ti46_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dSc43_log_pf_dT9 - dp_log_pf_dT9 - dTi46_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  60.7193303156212 + -35.69683001018291 * tfactors.T9i + -79.122 * tfactors.T913i + -0.495539 * tfactors.T913
                             + -0.618414 * tfactors.T9 + 0.0175873 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  35.69683001018291 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -79.122 * tfactors.T943i + (1.0/3.0) * -0.495539 * tfactors.T923i
                                      + -0.618414 + (5.0/3.0) * 0.0175873 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ti46_p_to_He4_Sc43_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti46_p_to_He4_Sc43_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti46_He4_to_p_V49_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ti46 + He4 --> p + V49

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real V49_log_pf, dV49_log_pf_dT9;

        // interpolating V49 partition function
        get_partition_function_cached(V49, tfactors.T9, pf_cache, V49_log_pf, dV49_log_pf_dT9);

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real Ti46_log_pf, dTi46_log_pf_dT9;

        // interpolating Ti46 partition function
        get_partition_function_cached(Ti46, tfactors.T9, pf_cache, Ti46_log_pf, dTi46_log_pf_dT9);

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real net_log_pf = p_log_pf + V49_log_pf - He4_log_pf - Ti46_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dV49_log_pf_dT9 - dHe4_log_pf_dT9 - dTi46_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  -19.293876950988306 + -11.954106481261894 * tfactors.T9i + -34.217 * tfactors.T913i + 47.2584 * tfactors.T913
                             + -9.18871 * tfactors.T9 + 0.87823 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  11.954106481261894 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -34.217 * tfactors.T943i + (1.0/3.0) * 47.2584 * tfactors.T923i
                                      + -9.18871 + (5.0/3.0) * 0.87823 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ti46_He4_to_p_V49_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti46_He4_to_p_V49_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V45_p_to_He4_Ti42_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // V45 + p --> He4 + Ti42

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Ti42_log_pf, dTi42_log_pf_dT9;

        // interpolating Ti42 partition function
        get_partition_function_cached(Ti42, tfactors.T9, pf_cache, Ti42_log_pf, dTi42_log_pf_dT9);

        amrex::Real V45_log_pf, dV45_log_pf_dT9;

        // interpolating V45 partition function
        get_partition_function_cached(V45, tfactors.T9, pf_cache, V45_log_pf, dV45_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Ti42_log_pf - p_log_pf - V45_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dTi42_log_pf_dT9 - dp_log_pf_dT9 - dV45_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  57.630501331425855 + -22.257405483263856 * tfactors.T9i + -81.5598 * tfactors.T913i + -1.03231 * tfactors.T913
                             + -0.471149 * tfactors.T9 + 0.0112292 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  22.257405483263856 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -81.5598 * tfactors.T943i + (1.0/3.0) * -1.03231 * tfactors.T923i
                                      + -0.471149 + (5.0/3.0) * 0.0112292 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_V45_p_to_He4_Ti42_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V45_p_to_He4_Ti42_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V47_p_to_He4_Ti44_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // V47 + p --> He4 + Ti44

        amrex::Real rate{}, drate_dT{};
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
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  6.086479266444236 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 70.2835 * tfactors.T923i
                                      + -7.99061 + (5.0/3.0) * 0.486213 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_V47_p_to_He4_Ti44_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V47_p_to_He4_Ti44_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr48_p_to_He4_V45_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Cr48 + p --> He4 + V45

        amrex::Real rate{}, drate_dT{};
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

        amrex::Real V45_log_pf, dV45_log_pf_dT9;

        // interpolating V45 partition function
        get_partition_function_cached(V45, tfactors.T9, pf_cache, V45_log_pf, dV45_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + V45_log_pf - p_log_pf - Cr48_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dV45_log_pf_dT9 - dp_log_pf_dT9 - dCr48_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  63.03246114679632 + -70.44574720883793 * tfactors.T9i + -84.1758 * tfactors.T913i + -1.12403 * tfactors.T913
                             + -0.469581 * tfactors.T9 + 0.00932056 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  70.44574720883793 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -84.1758 * tfactors.T943i + (1.0/3.0) * -1.12403 * tfactors.T923i
                                      + -0.469581 + (5.0/3.0) * 0.00932056 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Cr48_p_to_He4_V45_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr48_p_to_He4_V45_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr50_p_to_He4_V47_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Cr50 + p --> He4 + V47

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Cr50_log_pf, dCr50_log_pf_dT9;

        // interpolating Cr50 partition function
        get_partition_function_cached(Cr50, tfactors.T9, pf_cache, Cr50_log_pf, dCr50_log_pf_dT9);

        amrex::Real V47_log_pf, dV47_log_pf_dT9;

        // interpolating V47 partition function
        get_partition_function_cached(V47, tfactors.T9, pf_cache, V47_log_pf, dV47_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + V47_log_pf - p_log_pf - Cr50_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dV47_log_pf_dT9 - dp_log_pf_dT9 - dCr50_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  62.57281513112344 + -39.342041242503285 * tfactors.T9i + -84.2732 * tfactors.T913i + -1.10027 * tfactors.T913
                             + -0.487084 * tfactors.T9 + 0.00453269 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  39.342041242503285 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -84.2732 * tfactors.T943i + (1.0/3.0) * -1.10027 * tfactors.T923i
                                      + -0.487084 + (5.0/3.0) * 0.00453269 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Cr50_p_to_He4_V47_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr50_p_to_He4_V47_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr50_He4_to_p_Mn53_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Cr50 + He4 --> p + Mn53

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real Mn53_log_pf, dMn53_log_pf_dT9;

        // interpolating Mn53 partition function
        get_partition_function_cached(Mn53, tfactors.T9, pf_cache, Mn53_log_pf, dMn53_log_pf_dT9);

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Cr50_log_pf, dCr50_log_pf_dT9;

        // interpolating Cr50 partition function
        get_partition_function_cached(Cr50, tfactors.T9, pf_cache, Cr50_log_pf, dCr50_log_pf_dT9);

        amrex::Real net_log_pf = p_log_pf + Mn53_log_pf - He4_log_pf - Cr50_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dMn53_log_pf_dT9 - dHe4_log_pf_dT9 - dCr50_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  -67.54406116823061 + -5.049302153785043 * tfactors.T9i + -36.1913 * tfactors.T913i + 99.0857 * tfactors.T913
                             + -20.8943 * tfactors.T9 + 2.06209 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  5.049302153785043 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -36.1913 * tfactors.T943i + (1.0/3.0) * 99.0857 * tfactors.T923i
                                      + -20.8943 + (5.0/3.0) * 2.06209 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Cr50_He4_to_p_Mn53_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr50_He4_to_p_Mn53_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr52_p_to_He4_V49_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Cr52 + p --> He4 + V49

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real V49_log_pf, dV49_log_pf_dT9;

        // interpolating V49 partition function
        get_partition_function_cached(V49, tfactors.T9, pf_cache, V49_log_pf, dV49_log_pf_dT9);

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Cr52_log_pf, dCr52_log_pf_dT9;

        // interpolating Cr52 partition function
        get_partition_function_cached(Cr52, tfactors.T9, pf_cache, Cr52_log_pf, dCr52_log_pf_dT9);

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real net_log_pf = He4_log_pf + V49_log_pf - p_log_pf - Cr52_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dV49_log_pf_dT9 - dp_log_pf_dT9 - dCr52_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  60.786974428474686 + -30.0934723900715 * tfactors.T9i + -84.3633 * tfactors.T913i + 2.73419 * tfactors.T913
                             + -2.05655 * tfactors.T9 + 0.203834 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  30.0934723900715 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -84.3633 * tfactors.T943i + (1.0/3.0) * 2.73419 * tfactors.T923i
                                      + -2.05655 + (5.0/3.0) * 0.203834 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Cr52_p_to_He4_V49_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr52_p_to_He4_V49_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn49_p_to_He4_Cr46_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Mn49 + p --> He4 + Cr46

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Mn49_log_pf, dMn49_log_pf_dT9;

        // interpolating Mn49 partition function
        get_partition_function_cached(Mn49, tfactors.T9, pf_cache, Mn49_log_pf, dMn49_log_pf_dT9);

        amrex::Real Cr46_log_pf, dCr46_log_pf_dT9;

        // interpolating Cr46 partition function
        get_partition_function_cached(Cr46, tfactors.T9, pf_cache, Cr46_log_pf, dCr46_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Cr46_log_pf - p_log_pf - Mn49_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCr46_log_pf_dT9 - dp_log_pf_dT9 - dMn49_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  60.051147868613754 + -38.10743655961044 * tfactors.T9i + -86.6499 * tfactors.T913i + -1.08241 * tfactors.T913
                             + -0.516569 * tfactors.T9 + 0.0131721 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  38.10743655961044 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -86.6499 * tfactors.T943i + (1.0/3.0) * -1.08241 * tfactors.T923i
                                      + -0.516569 + (5.0/3.0) * 0.0131721 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Mn49_p_to_He4_Cr46_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn49_p_to_He4_Cr46_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn51_p_to_He4_Cr48_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Mn51 + p --> He4 + Cr48

        amrex::Real rate{}, drate_dT{};
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
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  6.47700151570639 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -86.7459 * tfactors.T943i + (1.0/3.0) * 1.05653 * tfactors.T923i
                                      + -1.15757 + (5.0/3.0) * 0.0877546 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Mn51_p_to_He4_Cr48_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn51_p_to_He4_Cr48_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn53_He4_to_n_Co56_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Mn53 + He4 --> n + Co56

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Mn53_log_pf, dMn53_log_pf_dT9;

        // interpolating Mn53 partition function
        get_partition_function_cached(Mn53, tfactors.T9, pf_cache, Mn53_log_pf, dMn53_log_pf_dT9);

        amrex::Real Co56_log_pf, dCo56_log_pf_dT9;

        // interpolating Co56 partition function
        get_partition_function_cached(Co56, tfactors.T9, pf_cache, Co56_log_pf, dCo56_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Co56_log_pf - He4_log_pf - Mn53_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dCo56_log_pf_dT9 - dHe4_log_pf_dT9 - dMn53_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  11.401304558489999 + -49.85535663187797 * tfactors.T9i + -2.81854 * tfactors.T913
                             + 1.84188 * tfactors.T9 + -0.164509 * tfactors.T953;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  49.85535663187797 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -2.81854 * tfactors.T923i
                                      + 1.84188 + (5.0/3.0) * -0.164509 * tfactors.T923;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Mn53_He4_to_n_Co56_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn53_He4_to_n_Co56_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe52_p_to_He4_Mn49_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Fe52 + p --> He4 + Mn49

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Fe52_log_pf, dFe52_log_pf_dT9;

        // interpolating Fe52 partition function
        get_partition_function_cached(Fe52, tfactors.T9, pf_cache, Fe52_log_pf, dFe52_log_pf_dT9);

        amrex::Real Mn49_log_pf, dMn49_log_pf_dT9;

        // interpolating Mn49 partition function
        get_partition_function_cached(Mn49, tfactors.T9, pf_cache, Mn49_log_pf, dMn49_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Mn49_log_pf - p_log_pf - Fe52_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dMn49_log_pf_dT9 - dp_log_pf_dT9 - dFe52_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  65.50849317033268 + -67.86490237870495 * tfactors.T9i + -89.1861 * tfactors.T913i + -1.31711 * tfactors.T913
                             + -0.507338 * tfactors.T9 + 0.0108137 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  67.86490237870495 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -89.1861 * tfactors.T943i + (1.0/3.0) * -1.31711 * tfactors.T923i
                                      + -0.507338 + (5.0/3.0) * 0.0108137 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Fe52_p_to_He4_Mn49_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe52_p_to_He4_Mn49_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe54_p_to_He4_Mn51_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Fe54 + p --> He4 + Mn51

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  36.52341983603332 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -89.274 * tfactors.T943i + (1.0/3.0) * -0.862452 * tfactors.T923i
                                      + -0.635672 + (5.0/3.0) * 0.0196464 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Fe54_p_to_He4_Mn51_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe54_p_to_He4_Mn51_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe54_He4_to_p_Co57_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Fe54 + He4 --> p + Co57

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  20.574290677453938 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -38.1133 * tfactors.T943i + (1.0/3.0) * 29.3541 * tfactors.T923i
                                      + -4.75966 + (5.0/3.0) * 0.40418 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Fe54_He4_to_p_Co57_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe54_He4_to_p_Co57_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe56_p_to_n_Co56_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Fe56 + p --> n + Co56

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  62.07264908160129 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.13331 * tfactors.T923i
                                      + 0.347185 + (5.0/3.0) * -0.0328879 * tfactors.T923;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Fe56_p_to_n_Co56_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe56_p_to_n_Co56_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe56_p_to_He4_Mn53_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Fe56 + p --> He4 + Mn53

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Mn53_log_pf, dMn53_log_pf_dT9;

        // interpolating Mn53 partition function
        get_partition_function_cached(Mn53, tfactors.T9, pf_cache, Mn53_log_pf, dMn53_log_pf_dT9);

        amrex::Real Fe56_log_pf, dFe56_log_pf_dT9;

        // interpolating Fe56 partition function
        get_partition_function_cached(Fe56, tfactors.T9, pf_cache, Fe56_log_pf, dFe56_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Mn53_log_pf - p_log_pf - Fe56_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dMn53_log_pf_dT9 - dp_log_pf_dT9 - dFe56_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  57.383021665674285 + -12.217292449723331 * tfactors.T9i + -89.3558 * tfactors.T913i + 9.25287 * tfactors.T913
                             + -3.70352 * tfactors.T9 + 0.372202 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  12.217292449723331 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -89.3558 * tfactors.T943i + (1.0/3.0) * 9.25287 * tfactors.T923i
                                      + -3.70352 + (5.0/3.0) * 0.372202 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Fe56_p_to_He4_Mn53_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe56_p_to_He4_Mn53_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co53_p_to_He4_Fe50_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Co53 + p --> He4 + Fe50

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real Fe50_log_pf, dFe50_log_pf_dT9;

        // interpolating Fe50 partition function
        get_partition_function_cached(Fe50, tfactors.T9, pf_cache, Fe50_log_pf, dFe50_log_pf_dT9);

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Co53_log_pf, dCo53_log_pf_dT9;

        // interpolating Co53 partition function
        get_partition_function_cached(Co53, tfactors.T9, pf_cache, Co53_log_pf, dCo53_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Fe50_log_pf - p_log_pf - Co53_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dFe50_log_pf_dT9 - dp_log_pf_dT9 - dCo53_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  61.93239942126981 + -38.51939695293296 * tfactors.T9i + -91.5952 * tfactors.T913i + -1.33713 * tfactors.T913
                             + -0.465416 * tfactors.T9 + 0.00761392 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  38.51939695293296 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -91.5952 * tfactors.T943i + (1.0/3.0) * -1.33713 * tfactors.T923i
                                      + -0.465416 + (5.0/3.0) * 0.00761392 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Co53_p_to_He4_Fe50_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co53_p_to_He4_Fe50_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co55_p_to_He4_Fe52_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Co55 + p --> He4 + Fe52

        amrex::Real rate{}, drate_dT{};
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
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  9.676367161770125 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -91.6819 * tfactors.T943i + (1.0/3.0) * -0.329235 * tfactors.T923i
                                      + -0.780924 + (5.0/3.0) * 0.0425179 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Co55_p_to_He4_Fe52_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co55_p_to_He4_Fe52_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co56_p_to_n_Ni56_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Co56 + p --> n + Ni56

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  33.830037058152215 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 1.76846 * tfactors.T923i
                                      + 0.197992 + (5.0/3.0) * -0.017494 * tfactors.T923;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Co56_p_to_n_Ni56_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co56_p_to_n_Ni56_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co56_He4_to_n_Cu59_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Co56 + He4 --> n + Cu59

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  61.8335402370445 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.00882 * tfactors.T923i
                                      + 1.88993 + (5.0/3.0) * -0.202319 * tfactors.T923;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Co56_He4_to_n_Cu59_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co56_He4_to_n_Cu59_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni56_p_to_He4_Co53_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ni56 + p --> He4 + Co53

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real Co53_log_pf, dCo53_log_pf_dT9;

        // interpolating Co53 partition function
        get_partition_function_cached(Co53, tfactors.T9, pf_cache, Co53_log_pf, dCo53_log_pf_dT9);

        amrex::Real He4_log_pf, dHe4_log_pf_dT9;

        // setting He4 log(partition function) to 0.0 by default, independent of T
        He4_log_pf = 0.0_rt;
        dHe4_log_pf_dT9 = 0.0_rt;

        amrex::Real Ni56_log_pf, dNi56_log_pf_dT9;

        // interpolating Ni56 partition function
        get_partition_function_cached(Ni56, tfactors.T9, pf_cache, Ni56_log_pf, dNi56_log_pf_dT9);

        amrex::Real net_log_pf = He4_log_pf + Co53_log_pf - p_log_pf - Ni56_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dHe4_log_pf_dT9 + dCo53_log_pf_dT9 - dp_log_pf_dT9 - dNi56_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  67.61201006170074 + -74.08492409178784 * tfactors.T9i + -94.0605 * tfactors.T913i + -1.48203 * tfactors.T913
                             + -0.476442 * tfactors.T9 + 0.00631587 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  74.08492409178784 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -94.0605 * tfactors.T943i + (1.0/3.0) * -1.48203 * tfactors.T923i
                                      + -0.476442 + (5.0/3.0) * 0.00631587 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ni56_p_to_He4_Co53_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni56_p_to_He4_Co53_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni56_He4_to_p_Cu59_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ni56 + He4 --> p + Cu59

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  28.003503178913892 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -39.9808 * tfactors.T943i + (1.0/3.0) * 13.6773 * tfactors.T923i
                                      + -3.76429 + (5.0/3.0) * 0.438096 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ni56_He4_to_p_Cu59_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni56_He4_to_p_Cu59_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni58_p_to_He4_Co55_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ni58 + p --> He4 + Co55

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  15.490230740645051 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -94.1404 * tfactors.T943i + (1.0/3.0) * 3.39179 * tfactors.T923i
                                      + -1.71062 + (5.0/3.0) * 0.133003 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ni58_p_to_He4_Co55_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni58_p_to_He4_Co55_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni57_to_n_Ni56_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ni57 --> n + Ni56

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  118.91983039605458 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.19665 * tfactors.T923i
                                      + 0.507179 + (5.0/3.0) * -0.074604 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ni57_to_n_Ni56_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni57_to_n_Ni56_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni58_to_n_Ni57_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ni58 --> n + Ni57

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  141.76332431826057 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.90814 * tfactors.T923i
                                      + 0.493188 + (5.0/3.0) * -0.0684633 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ni58_to_n_Ni57_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni58_to_n_Ni57_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co54_to_n_Co53_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Co54 --> n + Co53

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Co53_log_pf, dCo53_log_pf_dT9;

        // interpolating Co53 partition function
        get_partition_function_cached(Co53, tfactors.T9, pf_cache, Co53_log_pf, dCo53_log_pf_dT9);

        amrex::Real Co54_log_pf, dCo54_log_pf_dT9;

        // interpolating Co54 partition function
        get_partition_function_cached(Co54, tfactors.T9, pf_cache, Co54_log_pf, dCo54_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Co53_log_pf - Co54_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dCo53_log_pf_dT9 - dCo54_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  36.98265527996963 + -155.75605226924756 * tfactors.T9i + -1.14795 * tfactors.T913
                             + 0.413083 * tfactors.T9 + -0.0515627 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  155.75605226924756 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.14795 * tfactors.T923i
                                      + 0.413083 + (5.0/3.0) * -0.0515627 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Co54_to_n_Co53_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co54_to_n_Co53_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co55_to_n_Co54_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Co55 --> n + Co54

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Co54_log_pf, dCo54_log_pf_dT9;

        // interpolating Co54 partition function
        get_partition_function_cached(Co54, tfactors.T9, pf_cache, Co54_log_pf, dCo54_log_pf_dT9);

        amrex::Real Co55_log_pf, dCo55_log_pf_dT9;

        // interpolating Co55 partition function
        get_partition_function_cached(Co55, tfactors.T9, pf_cache, Co55_log_pf, dCo55_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Co54_log_pf - Co55_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dCo54_log_pf_dT9 - dCo55_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  35.99671117382933 + -163.52179579618158 * tfactors.T9i + 0.433187 * tfactors.T913
                             + 0.0888642 * tfactors.T9 + -0.0315521 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  163.52179579618158 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 0.433187 * tfactors.T923i
                                      + 0.0888642 + (5.0/3.0) * -0.0315521 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Co55_to_n_Co54_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co55_to_n_Co54_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe51_to_n_Fe50_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Fe51 --> n + Fe50

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Fe50_log_pf, dFe50_log_pf_dT9;

        // interpolating Fe50 partition function
        get_partition_function_cached(Fe50, tfactors.T9, pf_cache, Fe50_log_pf, dFe50_log_pf_dT9);

        amrex::Real Fe51_log_pf, dFe51_log_pf_dT9;

        // interpolating Fe51 partition function
        get_partition_function_cached(Fe51, tfactors.T9, pf_cache, Fe51_log_pf, dFe51_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Fe50_log_pf - Fe51_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dFe50_log_pf_dT9 - dFe51_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  35.13682696994766 + -159.96269008828872 * tfactors.T9i + -1.43737 * tfactors.T913
                             + 0.511909 * tfactors.T9 + -0.0531273 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  159.96269008828872 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.43737 * tfactors.T923i
                                      + 0.511909 + (5.0/3.0) * -0.0531273 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Fe51_to_n_Fe50_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe51_to_n_Fe50_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe52_to_n_Fe51_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Fe52 --> n + Fe51

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Fe51_log_pf, dFe51_log_pf_dT9;

        // interpolating Fe51 partition function
        get_partition_function_cached(Fe51, tfactors.T9, pf_cache, Fe51_log_pf, dFe51_log_pf_dT9);

        amrex::Real Fe52_log_pf, dFe52_log_pf_dT9;

        // interpolating Fe52 partition function
        get_partition_function_cached(Fe52, tfactors.T9, pf_cache, Fe52_log_pf, dFe52_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Fe51_log_pf - Fe52_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dFe51_log_pf_dT9 - dFe52_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  36.49553809630894 + -188.15818776827444 * tfactors.T9i + -0.566649 * tfactors.T913
                             + 0.19001 * tfactors.T9 + -0.0219395 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  188.15818776827444 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.566649 * tfactors.T923i
                                      + 0.19001 + (5.0/3.0) * -0.0219395 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Fe52_to_n_Fe51_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe52_to_n_Fe51_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe53_to_n_Fe52_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Fe53 --> n + Fe52

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  124.01421385132892 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.344319 * tfactors.T923i
                                      + 0.178277 + (5.0/3.0) * -0.0334326 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Fe53_to_n_Fe52_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe53_to_n_Fe52_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe54_to_n_Fe53_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Fe54 --> n + Fe53

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  155.25009527915816 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.10421 * tfactors.T923i
                                      + 0.379905 + (5.0/3.0) * -0.0581878 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Fe54_to_n_Fe53_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe54_to_n_Fe53_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe55_to_n_Fe54_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Fe55 --> n + Fe54

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  107.84100838776591 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -8.66617 * tfactors.T943i + (1.0/3.0) * 26.4472 * tfactors.T923i
                                      + -1.9222 + (5.0/3.0) * 0.0986404 * tfactors.T923 + -8.28317 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Fe55_to_n_Fe54_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe55_to_n_Fe54_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe56_to_n_Fe55_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Fe56 --> n + Fe55

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  130.0322634199283 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 8.06062 * tfactors.T943i + (1.0/3.0) * -14.4809 * tfactors.T923i
                                      + 0.94252 + (5.0/3.0) * -0.0776007 * tfactors.T923 + 7.97093 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Fe56_to_n_Fe55_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe56_to_n_Fe55_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn50_to_n_Mn49_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Mn50 --> n + Mn49

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Mn50_log_pf, dMn50_log_pf_dT9;

        // interpolating Mn50 partition function
        get_partition_function_cached(Mn50, tfactors.T9, pf_cache, Mn50_log_pf, dMn50_log_pf_dT9);

        amrex::Real Mn49_log_pf, dMn49_log_pf_dT9;

        // interpolating Mn49 partition function
        get_partition_function_cached(Mn49, tfactors.T9, pf_cache, Mn49_log_pf, dMn49_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Mn49_log_pf - Mn50_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dMn49_log_pf_dT9 - dMn50_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  36.75800895729457 + -151.7674633456902 * tfactors.T9i + -0.759512 * tfactors.T913
                             + 0.246336 * tfactors.T9 + -0.0261234 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  151.7674633456902 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.759512 * tfactors.T923i
                                      + 0.246336 + (5.0/3.0) * -0.0261234 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Mn50_to_n_Mn49_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn50_to_n_Mn49_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn51_to_n_Mn50_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Mn51 --> n + Mn50

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  158.83832832749732 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.538879 * tfactors.T923i
                                      + 0.284528 + (5.0/3.0) * -0.0459849 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Mn51_to_n_Mn50_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn51_to_n_Mn50_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn52_to_n_Mn51_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Mn52 --> n + Mn51

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  122.30591273875024 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.83611 * tfactors.T923i
                                      + 0.460384 + (5.0/3.0) * -0.0584947 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Mn52_to_n_Mn51_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn52_to_n_Mn51_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn53_to_n_Mn52_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Mn53 --> n + Mn52

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Mn52_log_pf, dMn52_log_pf_dT9;

        // interpolating Mn52 partition function
        get_partition_function_cached(Mn52, tfactors.T9, pf_cache, Mn52_log_pf, dMn52_log_pf_dT9);

        amrex::Real Mn53_log_pf, dMn53_log_pf_dT9;

        // interpolating Mn53 partition function
        get_partition_function_cached(Mn53, tfactors.T9, pf_cache, Mn53_log_pf, dMn53_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Mn52_log_pf - Mn53_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dMn52_log_pf_dT9 - dMn53_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  41.278513591470265 + -139.8370903552252 * tfactors.T9i + -0.91084 * tfactors.T913
                             + 0.162511 * tfactors.T9 + -0.0279619 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  139.8370903552252 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.91084 * tfactors.T923i
                                      + 0.162511 + (5.0/3.0) * -0.0279619 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Mn53_to_n_Mn52_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn53_to_n_Mn52_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr47_to_n_Cr46_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Cr47 --> n + Cr46

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Cr47_log_pf, dCr47_log_pf_dT9;

        // interpolating Cr47 partition function
        get_partition_function_cached(Cr47, tfactors.T9, pf_cache, Cr47_log_pf, dCr47_log_pf_dT9);

        amrex::Real Cr46_log_pf, dCr46_log_pf_dT9;

        // interpolating Cr46 partition function
        get_partition_function_cached(Cr46, tfactors.T9, pf_cache, Cr46_log_pf, dCr46_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Cr46_log_pf - Cr47_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dCr46_log_pf_dT9 - dCr47_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  37.791226447021984 + -152.74235891313134 * tfactors.T9i + -1.5945 * tfactors.T913
                             + 0.39528 * tfactors.T9 + -0.0372417 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  152.74235891313134 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.5945 * tfactors.T923i
                                      + 0.39528 + (5.0/3.0) * -0.0372417 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Cr47_to_n_Cr46_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr47_to_n_Cr46_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr48_to_n_Cr47_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Cr48 --> n + Cr47

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Cr48_log_pf, dCr48_log_pf_dT9;

        // interpolating Cr48 partition function
        get_partition_function_cached(Cr48, tfactors.T9, pf_cache, Cr48_log_pf, dCr48_log_pf_dT9);

        amrex::Real Cr47_log_pf, dCr47_log_pf_dT9;

        // interpolating Cr47 partition function
        get_partition_function_cached(Cr47, tfactors.T9, pf_cache, Cr47_log_pf, dCr47_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Cr47_log_pf - Cr48_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dCr47_log_pf_dT9 - dCr48_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  37.66466103343327 + -189.49386780398908 * tfactors.T9i + -1.36332 * tfactors.T913
                             + 0.288027 * tfactors.T9 + -0.0333017 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  189.49386780398908 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.36332 * tfactors.T923i
                                      + 0.288027 + (5.0/3.0) * -0.0333017 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Cr48_to_n_Cr47_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr48_to_n_Cr47_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr49_to_n_Cr48_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Cr49 --> n + Cr48

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Cr48_log_pf, dCr48_log_pf_dT9;

        // interpolating Cr48 partition function
        get_partition_function_cached(Cr48, tfactors.T9, pf_cache, Cr48_log_pf, dCr48_log_pf_dT9);

        amrex::Real Cr49_log_pf, dCr49_log_pf_dT9;

        // interpolating Cr49 partition function
        get_partition_function_cached(Cr49, tfactors.T9, pf_cache, Cr49_log_pf, dCr49_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Cr48_log_pf - Cr49_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dCr48_log_pf_dT9 - dCr49_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  37.75912618027852 + -122.80734396679019 * tfactors.T9i + -0.000452357 * tfactors.T913
                             + 0.0412838 * tfactors.T9 + -0.0160865 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  122.80734396679019 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.000452357 * tfactors.T923i
                                      + 0.0412838 + (5.0/3.0) * -0.0160865 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Cr49_to_n_Cr48_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr49_to_n_Cr48_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr50_to_n_Cr49_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Cr50 --> n + Cr49

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Cr49_log_pf, dCr49_log_pf_dT9;

        // interpolating Cr49 partition function
        get_partition_function_cached(Cr49, tfactors.T9, pf_cache, Cr49_log_pf, dCr49_log_pf_dT9);

        amrex::Real Cr50_log_pf, dCr50_log_pf_dT9;

        // interpolating Cr50 partition function
        get_partition_function_cached(Cr50, tfactors.T9, pf_cache, Cr50_log_pf, dCr50_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Cr49_log_pf - Cr50_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dCr49_log_pf_dT9 - dCr50_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  43.324496345019064 + -150.86196279664787 * tfactors.T9i + -1.98927 * tfactors.T913
                             + 0.465124 * tfactors.T9 + -0.0582249 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  150.86196279664787 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.98927 * tfactors.T923i
                                      + 0.465124 + (5.0/3.0) * -0.0582249 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Cr50_to_n_Cr49_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr50_to_n_Cr49_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr51_to_n_Cr50_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Cr51 --> n + Cr50

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Cr50_log_pf, dCr50_log_pf_dT9;

        // interpolating Cr50 partition function
        get_partition_function_cached(Cr50, tfactors.T9, pf_cache, Cr50_log_pf, dCr50_log_pf_dT9);

        amrex::Real Cr51_log_pf, dCr51_log_pf_dT9;

        // interpolating Cr51 partition function
        get_partition_function_cached(Cr51, tfactors.T9, pf_cache, Cr51_log_pf, dCr51_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Cr50_log_pf - Cr51_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dCr50_log_pf_dT9 - dCr51_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ks03
        ln_set_rate =  36.724790017844484 + -107.48815208414811 * tfactors.T9i + 0.815037 * tfactors.T913i + -0.708732 * tfactors.T913
                             + 0.386295 * tfactors.T9 + -0.0615661 * tfactors.T953 + 1.5454627 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  107.48815208414811 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 0.815037 * tfactors.T943i + (1.0/3.0) * -0.708732 * tfactors.T923i
                                      + 0.386295 + (5.0/3.0) * -0.0615661 * tfactors.T923 + 1.5454627 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Cr51_to_n_Cr50_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr51_to_n_Cr50_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr52_to_n_Cr51_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Cr52 --> n + Cr51

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Cr51_log_pf, dCr51_log_pf_dT9;

        // interpolating Cr51 partition function
        get_partition_function_cached(Cr51, tfactors.T9, pf_cache, Cr51_log_pf, dCr51_log_pf_dT9);

        amrex::Real Cr52_log_pf, dCr52_log_pf_dT9;

        // interpolating Cr52 partition function
        get_partition_function_cached(Cr52, tfactors.T9, pf_cache, Cr52_log_pf, dCr52_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Cr51_log_pf - Cr52_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dCr51_log_pf_dT9 - dCr52_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ks03
        ln_set_rate =  45.90309315543206 + -139.8057853770977 * tfactors.T9i + 6.95563 * tfactors.T913i + -11.4485 * tfactors.T913
                             + 0.736956 * tfactors.T9 + -0.0695323 * tfactors.T953 + 6.81244 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  139.8057853770977 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 6.95563 * tfactors.T943i + (1.0/3.0) * -11.4485 * tfactors.T923i
                                      + 0.736956 + (5.0/3.0) * -0.0695323 * tfactors.T923 + 6.81244 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Cr52_to_n_Cr51_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr52_to_n_Cr51_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V46_to_n_V45_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // V46 --> n + V45

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real V45_log_pf, dV45_log_pf_dT9;

        // interpolating V45 partition function
        get_partition_function_cached(V45, tfactors.T9, pf_cache, V45_log_pf, dV45_log_pf_dT9);

        amrex::Real V46_log_pf, dV46_log_pf_dT9;

        // interpolating V46 partition function
        get_partition_function_cached(V46, tfactors.T9, pf_cache, V46_log_pf, dV46_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + V45_log_pf - V46_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dV45_log_pf_dT9 - dV46_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  38.338260632832615 + -153.88540394802433 * tfactors.T9i + -1.62632 * tfactors.T913
                             + 0.391829 * tfactors.T9 + -0.0403093 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  153.88540394802433 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.62632 * tfactors.T923i
                                      + 0.391829 + (5.0/3.0) * -0.0403093 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_V46_to_n_V45_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V46_to_n_V45_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V47_to_n_V46_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // V47 --> n + V46

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real V46_log_pf, dV46_log_pf_dT9;

        // interpolating V46 partition function
        get_partition_function_cached(V46, tfactors.T9, pf_cache, V46_log_pf, dV46_log_pf_dT9);

        amrex::Real V47_log_pf, dV47_log_pf_dT9;

        // interpolating V47 partition function
        get_partition_function_cached(V47, tfactors.T9, pf_cache, V47_log_pf, dV47_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + V46_log_pf - V47_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dV46_log_pf_dT9 - dV47_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  38.04950790813785 + -150.88760878174838 * tfactors.T9i + -0.692996 * tfactors.T913
                             + 0.316873 * tfactors.T9 + -0.0417235 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  150.88760878174838 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.692996 * tfactors.T923i
                                      + 0.316873 + (5.0/3.0) * -0.0417235 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_V47_to_n_V46_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V47_to_n_V46_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V48_to_n_V47_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // V48 --> n + V47

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real V47_log_pf, dV47_log_pf_dT9;

        // interpolating V47 partition function
        get_partition_function_cached(V47, tfactors.T9, pf_cache, V47_log_pf, dV47_log_pf_dT9);

        amrex::Real V48_log_pf, dV48_log_pf_dT9;

        // interpolating V48 partition function
        get_partition_function_cached(V48, tfactors.T9, pf_cache, V48_log_pf, dV48_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + V47_log_pf - V48_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dV47_log_pf_dT9 - dV48_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  40.436203741894104 + -122.33770911835394 * tfactors.T9i + -2.40159 * tfactors.T913
                             + 0.594573 * tfactors.T9 + -0.0682896 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  122.33770911835394 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -2.40159 * tfactors.T923i
                                      + 0.594573 + (5.0/3.0) * -0.0682896 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_V48_to_n_V47_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V48_to_n_V47_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V49_to_n_V48_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // V49 --> n + V48

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real V49_log_pf, dV49_log_pf_dT9;

        // interpolating V49 partition function
        get_partition_function_cached(V49, tfactors.T9, pf_cache, V49_log_pf, dV49_log_pf_dT9);

        amrex::Real V48_log_pf, dV48_log_pf_dT9;

        // interpolating V48 partition function
        get_partition_function_cached(V48, tfactors.T9, pf_cache, V48_log_pf, dV48_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + V48_log_pf - V49_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dV48_log_pf_dT9 - dV49_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  41.69782013403119 + -134.0962191953525 * tfactors.T9i + -0.917026 * tfactors.T913
                             + -0.109162 * tfactors.T9 + 0.0127488 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  134.0962191953525 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.917026 * tfactors.T923i
                                      + -0.109162 + (5.0/3.0) * 0.0127488 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_V49_to_n_V48_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V49_to_n_V48_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti43_to_n_Ti42_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ti43 --> n + Ti42

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Ti42_log_pf, dTi42_log_pf_dT9;

        // interpolating Ti42 partition function
        get_partition_function_cached(Ti42, tfactors.T9, pf_cache, Ti42_log_pf, dTi42_log_pf_dT9);

        amrex::Real Ti43_log_pf, dTi43_log_pf_dT9;

        // interpolating Ti43 partition function
        get_partition_function_cached(Ti43, tfactors.T9, pf_cache, Ti43_log_pf, dTi43_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Ti42_log_pf - Ti43_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dTi42_log_pf_dT9 - dTi43_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  35.246514232360624 + -142.53792590287344 * tfactors.T9i + -1.00856 * tfactors.T913
                             + 0.342738 * tfactors.T9 + -0.0367066 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  142.53792590287344 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.00856 * tfactors.T923i
                                      + 0.342738 + (5.0/3.0) * -0.0367066 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ti43_to_n_Ti42_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti43_to_n_Ti42_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti44_to_n_Ti43_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ti44 --> n + Ti43

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Ti43_log_pf, dTi43_log_pf_dT9;

        // interpolating Ti43 partition function
        get_partition_function_cached(Ti43, tfactors.T9, pf_cache, Ti43_log_pf, dTi43_log_pf_dT9);

        amrex::Real Ti44_log_pf, dTi44_log_pf_dT9;

        // interpolating Ti44 partition function
        get_partition_function_cached(Ti44, tfactors.T9, pf_cache, Ti44_log_pf, dTi44_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Ti43_log_pf - Ti44_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dTi43_log_pf_dT9 - dTi44_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  37.32554402156041 + -189.19911304371885 * tfactors.T9i + -1.09388 * tfactors.T913
                             + 0.248209 * tfactors.T9 + -0.02653 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  189.19911304371885 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.09388 * tfactors.T923i
                                      + 0.248209 + (5.0/3.0) * -0.02653 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ti44_to_n_Ti43_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti44_to_n_Ti43_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti45_to_n_Ti44_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ti45 --> n + Ti44

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Ti44_log_pf, dTi44_log_pf_dT9;

        // interpolating Ti44 partition function
        get_partition_function_cached(Ti44, tfactors.T9, pf_cache, Ti44_log_pf, dTi44_log_pf_dT9);

        amrex::Real Ti45_log_pf, dTi45_log_pf_dT9;

        // interpolating Ti45 partition function
        get_partition_function_cached(Ti45, tfactors.T9, pf_cache, Ti45_log_pf, dTi45_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Ti44_log_pf - Ti45_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dTi44_log_pf_dT9 - dTi45_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  36.864811289829866 + -110.62608129453406 * tfactors.T9i + 1.01203 * tfactors.T913
                             + -0.201174 * tfactors.T9 + 0.00360954 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  110.62608129453406 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 1.01203 * tfactors.T923i
                                      + -0.201174 + (5.0/3.0) * 0.00360954 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ti45_to_n_Ti44_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti45_to_n_Ti44_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti46_to_n_Ti45_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ti46 --> n + Ti45

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Ti46_log_pf, dTi46_log_pf_dT9;

        // interpolating Ti46 partition function
        get_partition_function_cached(Ti46, tfactors.T9, pf_cache, Ti46_log_pf, dTi46_log_pf_dT9);

        amrex::Real Ti45_log_pf, dTi45_log_pf_dT9;

        // interpolating Ti45 partition function
        get_partition_function_cached(Ti45, tfactors.T9, pf_cache, Ti45_log_pf, dTi45_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Ti45_log_pf - Ti46_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dTi45_log_pf_dT9 - dTi46_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  42.907047253582405 + -153.05533276685685 * tfactors.T9i + -1.55633 * tfactors.T913
                             + 0.300783 * tfactors.T9 + -0.0373696 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  153.05533276685685 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.55633 * tfactors.T923i
                                      + 0.300783 + (5.0/3.0) * -0.0373696 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ti46_to_n_Ti45_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti46_to_n_Ti45_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc42_to_n_Sc41_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Sc42 --> n + Sc41

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Sc42_log_pf, dSc42_log_pf_dT9;

        // interpolating Sc42 partition function
        get_partition_function_cached(Sc42, tfactors.T9, pf_cache, Sc42_log_pf, dSc42_log_pf_dT9);

        amrex::Real Sc41_log_pf, dSc41_log_pf_dT9;

        // interpolating Sc41 partition function
        get_partition_function_cached(Sc41, tfactors.T9, pf_cache, Sc41_log_pf, dSc41_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Sc41_log_pf - Sc42_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dSc41_log_pf_dT9 - dSc42_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  36.372198803455944 + -134.0316980746222 * tfactors.T9i + -1.01509 * tfactors.T913
                             + 0.362316 * tfactors.T9 + -0.0412181 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  134.0316980746222 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.01509 * tfactors.T923i
                                      + 0.362316 + (5.0/3.0) * -0.0412181 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Sc42_to_n_Sc41_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc42_to_n_Sc41_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc43_to_n_Sc42_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Sc43 --> n + Sc42

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Sc42_log_pf, dSc42_log_pf_dT9;

        // interpolating Sc42 partition function
        get_partition_function_cached(Sc42, tfactors.T9, pf_cache, Sc42_log_pf, dSc42_log_pf_dT9);

        amrex::Real Sc43_log_pf, dSc43_log_pf_dT9;

        // interpolating Sc43 partition function
        get_partition_function_cached(Sc43, tfactors.T9, pf_cache, Sc43_log_pf, dSc43_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Sc42_log_pf - Sc43_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dSc42_log_pf_dT9 - dSc43_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  36.55709480577273 + -140.8604928084524 * tfactors.T9i + -0.921604 * tfactors.T913
                             + 0.333782 * tfactors.T9 + -0.0404241 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  140.8604928084524 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.921604 * tfactors.T923i
                                      + 0.333782 + (5.0/3.0) * -0.0404241 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Sc43_to_n_Sc42_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc43_to_n_Sc42_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca39_to_n_Ca38_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ca39 --> n + Ca38

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Ca39_log_pf, dCa39_log_pf_dT9;

        // interpolating Ca39 partition function
        get_partition_function_cached(Ca39, tfactors.T9, pf_cache, Ca39_log_pf, dCa39_log_pf_dT9);

        amrex::Real Ca38_log_pf, dCa38_log_pf_dT9;

        // interpolating Ca38 partition function
        get_partition_function_cached(Ca38, tfactors.T9, pf_cache, Ca38_log_pf, dCa38_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Ca38_log_pf - Ca39_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dCa38_log_pf_dT9 - dCa39_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  34.369887995652356 + -154.28808072685763 * tfactors.T9i + -0.418748 * tfactors.T913
                             + 0.243656 * tfactors.T9 + -0.0268305 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  154.28808072685763 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.418748 * tfactors.T923i
                                      + 0.243656 + (5.0/3.0) * -0.0268305 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ca39_to_n_Ca38_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca39_to_n_Ca38_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca40_to_n_Ca39_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ca40 --> n + Ca39

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Ca40_log_pf, dCa40_log_pf_dT9;

        // interpolating Ca40 partition function
        get_partition_function_cached(Ca40, tfactors.T9, pf_cache, Ca40_log_pf, dCa40_log_pf_dT9);

        amrex::Real Ca39_log_pf, dCa39_log_pf_dT9;

        // interpolating Ca39 partition function
        get_partition_function_cached(Ca39, tfactors.T9, pf_cache, Ca39_log_pf, dCa39_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Ca39_log_pf - Ca40_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dCa39_log_pf_dT9 - dCa40_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  34.41095439853984 + -181.43687408133584 * tfactors.T9i + 0.436603 * tfactors.T913
                             + 0.098368 * tfactors.T9 + -0.0226294 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  181.43687408133584 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 0.436603 * tfactors.T923i
                                      + 0.098368 + (5.0/3.0) * -0.0226294 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ca40_to_n_Ca39_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca40_to_n_Ca39_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca41_to_n_Ca40_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ca41 --> n + Ca40

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Ca40_log_pf, dCa40_log_pf_dT9;

        // interpolating Ca40 partition function
        get_partition_function_cached(Ca40, tfactors.T9, pf_cache, Ca40_log_pf, dCa40_log_pf_dT9);

        amrex::Real Ca41_log_pf, dCa41_log_pf_dT9;

        // interpolating Ca41 partition function
        get_partition_function_cached(Ca41, tfactors.T9, pf_cache, Ca41_log_pf, dCa41_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Ca40_log_pf - Ca41_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dCa40_log_pf_dT9 - dCa41_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ks03
        ln_set_rate =  24.6568173607497 + -96.96102952477759 * tfactors.T9i + -9.18636 * tfactors.T913i + 21.0095 * tfactors.T913
                             + -1.33326 * tfactors.T9 + 0.0769347 * tfactors.T953 + -7.16245 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  96.96102952477759 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -9.18636 * tfactors.T943i + (1.0/3.0) * 21.0095 * tfactors.T923i
                                      + -1.33326 + (5.0/3.0) * 0.0769347 * tfactors.T923 + -7.16245 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ca41_to_n_Ca40_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca41_to_n_Ca40_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca42_to_n_Ca41_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ca42 --> n + Ca41

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Ca41_log_pf, dCa41_log_pf_dT9;

        // interpolating Ca41 partition function
        get_partition_function_cached(Ca41, tfactors.T9, pf_cache, Ca41_log_pf, dCa41_log_pf_dT9);

        amrex::Real Ca42_log_pf, dCa42_log_pf_dT9;

        // interpolating Ca42 partition function
        get_partition_function_cached(Ca42, tfactors.T9, pf_cache, Ca42_log_pf, dCa42_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Ca41_log_pf - Ca42_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dCa41_log_pf_dT9 - dCa42_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ks03
        ln_set_rate =  46.216884269737946 + -133.26891724956297 * tfactors.T9i + 3.7037 * tfactors.T913i + -9.99246 * tfactors.T913
                             + 1.05894 * tfactors.T9 + -0.0978564 * tfactors.T953 + 4.93604 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  133.26891724956297 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 3.7037 * tfactors.T943i + (1.0/3.0) * -9.99246 * tfactors.T923i
                                      + 1.05894 + (5.0/3.0) * -0.0978564 * tfactors.T923 + 4.93604 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ca42_to_n_Ca41_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca42_to_n_Ca41_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar35_to_n_Ar34_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ar35 --> n + Ar34

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Ar34_log_pf, dAr34_log_pf_dT9;

        // interpolating Ar34 partition function
        get_partition_function_cached(Ar34, tfactors.T9, pf_cache, Ar34_log_pf, dAr34_log_pf_dT9);

        amrex::Real Ar35_log_pf, dAr35_log_pf_dT9;

        // interpolating Ar35 partition function
        get_partition_function_cached(Ar35, tfactors.T9, pf_cache, Ar35_log_pf, dAr35_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Ar34_log_pf - Ar35_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dAr34_log_pf_dT9 - dAr35_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  32.997766916371916 + -147.84536831101275 * tfactors.T9i + 1.14834 * tfactors.T913
                             + -0.142306 * tfactors.T9 + 0.00918793 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  147.84536831101275 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 1.14834 * tfactors.T923i
                                      + -0.142306 + (5.0/3.0) * 0.00918793 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ar35_to_n_Ar34_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar35_to_n_Ar34_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar36_to_n_Ar35_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ar36 --> n + Ar35

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Ar35_log_pf, dAr35_log_pf_dT9;

        // interpolating Ar35 partition function
        get_partition_function_cached(Ar35, tfactors.T9, pf_cache, Ar35_log_pf, dAr35_log_pf_dT9);

        amrex::Real Ar36_log_pf, dAr36_log_pf_dT9;

        // interpolating Ar36 partition function
        get_partition_function_cached(Ar36, tfactors.T9, pf_cache, Ar36_log_pf, dAr36_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Ar35_log_pf - Ar36_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dAr35_log_pf_dT9 - dAr36_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  33.8372707432332 + -177.03342363488906 * tfactors.T9i + 0.542023 * tfactors.T913
                             + 0.0223868 * tfactors.T9 + -0.0120062 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  177.03342363488906 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 0.542023 * tfactors.T923i
                                      + 0.0223868 + (5.0/3.0) * -0.0120062 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ar36_to_n_Ar35_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar36_to_n_Ar35_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar37_to_n_Ar36_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ar37 --> n + Ar36

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Ar37_log_pf, dAr37_log_pf_dT9;

        // interpolating Ar37 partition function
        get_partition_function_cached(Ar37, tfactors.T9, pf_cache, Ar37_log_pf, dAr37_log_pf_dT9);

        amrex::Real Ar36_log_pf, dAr36_log_pf_dT9;

        // interpolating Ar36 partition function
        get_partition_function_cached(Ar36, tfactors.T9, pf_cache, Ar36_log_pf, dAr36_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Ar36_log_pf - Ar37_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dAr36_log_pf_dT9 - dAr37_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ks03
        ln_set_rate =  34.31482400258879 + -101.94248915476915 * tfactors.T9i + -3.1764 * tfactors.T913i + 5.13191 * tfactors.T913
                             + -0.00639688 * tfactors.T9 + -0.0292833 * tfactors.T953 + -1.24683 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  101.94248915476915 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -3.1764 * tfactors.T943i + (1.0/3.0) * 5.13191 * tfactors.T923i
                                      + -0.00639688 + (5.0/3.0) * -0.0292833 * tfactors.T923 + -1.24683 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ar37_to_n_Ar36_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar37_to_n_Ar36_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar38_to_n_Ar37_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ar38 --> n + Ar37

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Ar38_log_pf, dAr38_log_pf_dT9;

        // interpolating Ar38 partition function
        get_partition_function_cached(Ar38, tfactors.T9, pf_cache, Ar38_log_pf, dAr38_log_pf_dT9);

        amrex::Real Ar37_log_pf, dAr37_log_pf_dT9;

        // interpolating Ar37 partition function
        get_partition_function_cached(Ar37, tfactors.T9, pf_cache, Ar37_log_pf, dAr37_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Ar37_log_pf - Ar38_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dAr37_log_pf_dT9 - dAr38_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  39.86712497672913 + -137.37971759787655 * tfactors.T9i + -0.825362 * tfactors.T913
                             + 0.336634 * tfactors.T9 + -0.0509617 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  137.37971759787655 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.825362 * tfactors.T923i
                                      + 0.336634 + (5.0/3.0) * -0.0509617 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ar38_to_n_Ar37_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar38_to_n_Ar37_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si29_to_n_Si28_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Si29 --> n + Si28

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Si29_log_pf, dSi29_log_pf_dT9;

        // interpolating Si29 partition function
        get_partition_function_cached(Si29, tfactors.T9, pf_cache, Si29_log_pf, dSi29_log_pf_dT9);

        amrex::Real Si28_log_pf, dSi28_log_pf_dT9;

        // interpolating Si28 partition function
        get_partition_function_cached(Si28, tfactors.T9, pf_cache, Si28_log_pf, dSi28_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Si28_log_pf - Si29_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dSi28_log_pf_dT9 - dSi29_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ka02r
        ln_set_rate =  29.897307803823217 + -98.71207376605948 * tfactors.T9i + 7.68863 * tfactors.T913
                             + -1.7991 * tfactors.T9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  98.71207376605948 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 7.68863 * tfactors.T923i
                                      + -1.7991;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from ka02
        ln_set_rate =  31.757037803823216 + -98.33207376605948 * tfactors.T9i
                             + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  98.33207376605948 * tfactors.T9i * tfactors.T9i
                                      + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Si29_to_n_Si28_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si29_to_n_Si28_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si30_to_n_Si29_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Si30 --> n + Si29

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real Si30_log_pf, dSi30_log_pf_dT9;

        // interpolating Si30 partition function
        get_partition_function_cached(Si30, tfactors.T9, pf_cache, Si30_log_pf, dSi30_log_pf_dT9);

        amrex::Real Si29_log_pf, dSi29_log_pf_dT9;

        // interpolating Si29 partition function
        get_partition_function_cached(Si29, tfactors.T9, pf_cache, Si29_log_pf, dSi29_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + Si29_log_pf - Si30_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dSi29_log_pf_dT9 - dSi30_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ka02r
        ln_set_rate =  33.97093028791906 + -123.2940034088064 * tfactors.T9i + 5.50678 * tfactors.T913
                             + -2.85656 * tfactors.T9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  123.2940034088064 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 5.50678 * tfactors.T923i
                                      + -2.85656;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from ka02n
        ln_set_rate =  36.17208028791906 + -123.1146374088064 * tfactors.T9i + 0.650904 * tfactors.T913
                             + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  123.1146374088064 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 0.650904 * tfactors.T923i
                                      + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Si30_to_n_Si29_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si30_to_n_Si29_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S31_to_n_S30_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // S31 --> n + S30

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real S30_log_pf, dS30_log_pf_dT9;

        // interpolating S30 partition function
        get_partition_function_cached(S30, tfactors.T9, pf_cache, S30_log_pf, dS30_log_pf_dT9);

        amrex::Real S31_log_pf, dS31_log_pf_dT9;

        // interpolating S31 partition function
        get_partition_function_cached(S31, tfactors.T9, pf_cache, S31_log_pf, dS31_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + S30_log_pf - S31_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dS30_log_pf_dT9 - dS31_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  33.31388497305849 + -151.49232022100603 * tfactors.T9i + 1.62298 * tfactors.T913
                             + -0.278802 * tfactors.T9 + 0.0210647 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  151.49232022100603 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 1.62298 * tfactors.T923i
                                      + -0.278802 + (5.0/3.0) * 0.0210647 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_S31_to_n_S30_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S31_to_n_S30_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S32_to_n_S31_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // S32 --> n + S31

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real S31_log_pf, dS31_log_pf_dT9;

        // interpolating S31 partition function
        get_partition_function_cached(S31, tfactors.T9, pf_cache, S31_log_pf, dS31_log_pf_dT9);

        amrex::Real S32_log_pf, dS32_log_pf_dT9;

        // interpolating S32 partition function
        get_partition_function_cached(S32, tfactors.T9, pf_cache, S32_log_pf, dS32_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + S31_log_pf - S32_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dS31_log_pf_dT9 - dS32_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ths8r
        ln_set_rate =  31.939203653425814 + -174.58214440988277 * tfactors.T9i + 1.71463 * tfactors.T913
                             + -0.221804 * tfactors.T9 + 0.00880104 * tfactors.T953 + 1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  174.58214440988277 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 1.71463 * tfactors.T923i
                                      + -0.221804 + (5.0/3.0) * 0.00880104 * tfactors.T923 + 1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_S32_to_n_S31_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S32_to_n_S31_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S33_to_n_S32_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // S33 --> n + S32

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real S32_log_pf, dS32_log_pf_dT9;

        // interpolating S32 partition function
        get_partition_function_cached(S32, tfactors.T9, pf_cache, S32_log_pf, dS32_log_pf_dT9);

        amrex::Real S33_log_pf, dS33_log_pf_dT9;

        // interpolating S33 partition function
        get_partition_function_cached(S33, tfactors.T9, pf_cache, S33_log_pf, dS33_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + S32_log_pf - S33_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dS32_log_pf_dT9 - dS33_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ks03
        ln_set_rate =  34.741454198619536 + -100.0832318567133 * tfactors.T9i + -15.0178 * tfactors.T913i + 16.3567 * tfactors.T913
                             + -0.436839 * tfactors.T9 + -0.00574462 * tfactors.T953 + -8.28034 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  100.0832318567133 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -15.0178 * tfactors.T943i + (1.0/3.0) * 16.3567 * tfactors.T923i
                                      + -0.436839 + (5.0/3.0) * -0.00574462 * tfactors.T923 + -8.28034 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_S33_to_n_S32_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S33_to_n_S32_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S34_to_n_S33_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // S34 --> n + S33

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real n_log_pf, dn_log_pf_dT9;

        // setting n log(partition function) to 0.0 by default, independent of T
        n_log_pf = 0.0_rt;
        dn_log_pf_dT9 = 0.0_rt;

        amrex::Real S33_log_pf, dS33_log_pf_dT9;

        // interpolating S33 partition function
        get_partition_function_cached(S33, tfactors.T9, pf_cache, S33_log_pf, dS33_log_pf_dT9);

        amrex::Real S34_log_pf, dS34_log_pf_dT9;

        // interpolating S34 partition function
        get_partition_function_cached(S34, tfactors.T9, pf_cache, S34_log_pf, dS34_log_pf_dT9);

        amrex::Real net_log_pf = n_log_pf + S33_log_pf - S34_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dn_log_pf_dT9 + dS33_log_pf_dT9 - dS34_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from ks03
        ln_set_rate =  20.889721414503413 + -132.24411175059305 * tfactors.T9i + -23.1351 * tfactors.T913i + 43.3623 * tfactors.T913
                             + -2.74733 * tfactors.T9 + 0.159045 * tfactors.T953 + -18.1332 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  132.24411175059305 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -23.1351 * tfactors.T943i + (1.0/3.0) * 43.3623 * tfactors.T923i
                                      + -2.74733 + (5.0/3.0) * 0.159045 * tfactors.T923 + -18.1332 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_S34_to_n_S33_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S34_to_n_S33_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Na23_p_to_C12_C12_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Na23 + p --> C12 + C12

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  26.004360836807496 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -84.165 * tfactors.T943i + (1.0/3.0) * -1.4191 * tfactors.T923i
                                      + -0.114619 + (5.0/3.0) * -0.070307 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Na23_p_to_C12_C12_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Na23_p_to_C12_C12_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mg24_to_p_Na23_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Mg24 --> p + Na23

        amrex::Real rate{}, drate_dT{};
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Mg24_to_p_Na23_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg24_to_p_Na23_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ne20_He4_to_C12_C12_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ne20 + He4 --> C12 + C12

        amrex::Real rate{}, drate_dT{};
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

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  53.57824966896331 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -84.165 * tfactors.T943i + (1.0/3.0) * -1.56627 * tfactors.T923i
                                      + -0.0736084 + (5.0/3.0) * -0.072797 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ne20_He4_to_C12_C12_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne20_He4_to_C12_C12_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ne20_He4_to_p_Na23_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Ne20 + He4 --> p + Na23

        amrex::Real rate{}, drate_dT{};
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
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  29.889658832139027 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 19.7297 * tfactors.T923i
                                      + -2.20987 + (5.0/3.0) * 0.153374 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from il10r
        ln_set_rate =  0.24516522151364528 + -29.434918832139026 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  29.434918832139026 * tfactors.T9i * tfactors.T9i
                                      + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        // ReacLib set derived from il10n
        ln_set_rate =  19.202935721513644 + -27.573888832139026 * tfactors.T9i + -20.0024 * tfactors.T913i + 11.5988 * tfactors.T913
                             + -1.37398 * tfactors.T9 + -1.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += net_log_pf;
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  27.573888832139026 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -20.0024 * tfactors.T943i + (1.0/3.0) * 11.5988 * tfactors.T923i
                                      + -1.37398 + (5.0/3.0) * -1.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += net_dlog_pf_dT9;
            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Ne20_He4_to_p_Na23_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne20_He4_to_p_Na23_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mg24_to_He4_Ne20_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // Mg24 --> He4 + Ne20

        amrex::Real rate{}, drate_dT{};
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_Mg24_to_He4_Ne20_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg24_to_He4_Ne20_derived) = drate_dT;
        }
    }


    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void
    fill_rates([[maybe_unused]] const tf_t& tfactors,
               [[maybe_unused]] T& rate_eval)
    {

        [[maybe_unused]] amrex::Real rate;
        [[maybe_unused]] amrex::Real drate_dT;

        part_fun::pf_cache_t pf_cache{};

        pf_cache.index_temp_array_1 = interp_net::find_index(tfactors.T9, part_fun::temp_array_1);

        {
            // O16_to_He4_C12_derived
            rate_O16_to_He4_C12_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ne20_to_He4_O16_derived
            rate_Ne20_to_He4_O16_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mg24_to_n_Mg23_derived
            rate_Mg24_to_n_Mg23_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Si27_to_He4_Mg23_derived
            rate_Si27_to_He4_Mg23_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Si28_to_n_Si27_derived
            rate_Si28_to_n_Si27_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Si28_to_p_Al27_derived
            rate_Si28_to_p_Al27_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Si28_to_He4_Mg24_derived
            rate_Si28_to_He4_Mg24_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // P29_to_p_Si28_derived
            rate_P29_to_p_Si28_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // P31_to_p_Si30_derived
            rate_P31_to_p_Si30_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // P31_to_He4_Al27_derived
            rate_P31_to_He4_Al27_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // S30_to_p_P29_derived
            rate_S30_to_p_P29_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // S32_to_p_P31_derived
            rate_S32_to_p_P31_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // S32_to_He4_Si28_derived
            rate_S32_to_He4_Si28_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // S34_to_He4_Si30_derived
            rate_S34_to_He4_Si30_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cl33_to_p_S32_derived
            rate_Cl33_to_p_S32_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cl33_to_He4_P29_derived
            rate_Cl33_to_He4_P29_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cl35_to_p_S34_derived
            rate_Cl35_to_p_S34_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cl35_to_He4_P31_derived
            rate_Cl35_to_He4_P31_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ar34_to_p_Cl33_derived
            rate_Ar34_to_p_Cl33_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ar34_to_He4_S30_derived
            rate_Ar34_to_He4_S30_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ar36_to_p_Cl35_derived
            rate_Ar36_to_p_Cl35_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ar36_to_He4_S32_derived
            rate_Ar36_to_He4_S32_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ar38_to_He4_S34_derived
            rate_Ar38_to_He4_S34_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // K37_to_p_Ar36_derived
            rate_K37_to_p_Ar36_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // K37_to_He4_Cl33_derived
            rate_K37_to_He4_Cl33_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // K39_to_p_Ar38_derived
            rate_K39_to_p_Ar38_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // K39_to_He4_Cl35_derived
            rate_K39_to_He4_Cl35_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca38_to_p_K37_derived
            rate_Ca38_to_p_K37_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca38_to_He4_Ar34_derived
            rate_Ca38_to_He4_Ar34_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca40_to_p_K39_derived
            rate_Ca40_to_p_K39_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca40_to_He4_Ar36_derived
            rate_Ca40_to_He4_Ar36_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca42_to_He4_Ar38_derived
            rate_Ca42_to_He4_Ar38_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc41_to_p_Ca40_derived
            rate_Sc41_to_p_Ca40_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc41_to_He4_K37_derived
            rate_Sc41_to_He4_K37_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc43_to_p_Ca42_derived
            rate_Sc43_to_p_Ca42_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc43_to_He4_K39_derived
            rate_Sc43_to_He4_K39_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti42_to_p_Sc41_derived
            rate_Ti42_to_p_Sc41_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti42_to_He4_Ca38_derived
            rate_Ti42_to_He4_Ca38_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti44_to_p_Sc43_derived
            rate_Ti44_to_p_Sc43_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti44_to_He4_Ca40_derived
            rate_Ti44_to_He4_Ca40_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti46_to_He4_Ca42_derived
            rate_Ti46_to_He4_Ca42_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V45_to_p_Ti44_derived
            rate_V45_to_p_Ti44_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V45_to_He4_Sc41_derived
            rate_V45_to_He4_Sc41_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V47_to_p_Ti46_derived
            rate_V47_to_p_Ti46_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V47_to_He4_Sc43_derived
            rate_V47_to_He4_Sc43_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr46_to_p_V45_derived
            rate_Cr46_to_p_V45_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr46_to_He4_Ti42_derived
            rate_Cr46_to_He4_Ti42_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr48_to_p_V47_derived
            rate_Cr48_to_p_V47_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr48_to_He4_Ti44_derived
            rate_Cr48_to_He4_Ti44_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr50_to_p_V49_derived
            rate_Cr50_to_p_V49_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr50_to_He4_Ti46_derived
            rate_Cr50_to_He4_Ti46_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn49_to_p_Cr48_derived
            rate_Mn49_to_p_Cr48_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn49_to_He4_V45_derived
            rate_Mn49_to_He4_V45_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn51_to_p_Cr50_derived
            rate_Mn51_to_p_Cr50_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn51_to_He4_V47_derived
            rate_Mn51_to_He4_V47_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn53_to_p_Cr52_derived
            rate_Mn53_to_p_Cr52_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn53_to_He4_V49_derived
            rate_Mn53_to_He4_V49_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe50_to_p_Mn49_derived
            rate_Fe50_to_p_Mn49_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe50_to_He4_Cr46_derived
            rate_Fe50_to_He4_Cr46_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe52_to_p_Mn51_derived
            rate_Fe52_to_p_Mn51_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe52_to_He4_Cr48_derived
            rate_Fe52_to_He4_Cr48_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe54_to_p_Mn53_derived
            rate_Fe54_to_p_Mn53_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe54_to_He4_Cr50_derived
            rate_Fe54_to_He4_Cr50_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe56_to_He4_Cr52_derived
            rate_Fe56_to_He4_Cr52_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co53_to_p_Fe52_derived
            rate_Co53_to_p_Fe52_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co53_to_He4_Mn49_derived
            rate_Co53_to_He4_Mn49_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co55_to_p_Fe54_derived
            rate_Co55_to_p_Fe54_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co55_to_He4_Mn51_derived
            rate_Co55_to_He4_Mn51_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co56_to_n_Co55_derived
            rate_Co56_to_n_Co55_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co57_to_n_Co56_derived
            rate_Co57_to_n_Co56_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co57_to_p_Fe56_derived
            rate_Co57_to_p_Fe56_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co57_to_He4_Mn53_derived
            rate_Co57_to_He4_Mn53_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ni56_to_p_Co55_derived
            rate_Ni56_to_p_Co55_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ni56_to_He4_Fe52_derived
            rate_Ni56_to_He4_Fe52_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ni58_to_p_Co57_derived
            rate_Ni58_to_p_Co57_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ni58_to_He4_Fe54_derived
            rate_Ni58_to_He4_Fe54_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cu59_to_p_Ni58_derived
            rate_Cu59_to_p_Ni58_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cu59_to_He4_Co55_derived
            rate_Cu59_to_He4_Co55_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Zn60_to_p_Cu59_derived
            rate_Zn60_to_p_Cu59_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Zn60_to_He4_Ni56_derived
            rate_Zn60_to_He4_Ni56_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // C12_to_3He4_derived
            rate_C12_to_3He4_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // C12_C12_to_n_Mg23_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_C12_C12);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_C12);
            }
#endif
            rate_C12_C12_to_n_Mg23_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // O16_C12_to_n_Si27_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_C12_O16);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_O16);
            }
#endif
            rate_O16_C12_to_n_Si27_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ne20_He4_to_n_Mg23_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ne20);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ne20);
            }
#endif
            rate_Ne20_He4_to_n_Mg23_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Mg24_He4_to_n_Si27_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mg24);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mg24);
            }
#endif
            rate_Mg24_He4_to_n_Si27_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Mg24_He4_to_p_Al27_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mg24);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mg24);
            }
#endif
            rate_Mg24_He4_to_p_Al27_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Mg24_He4_to_C12_O16_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mg24);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mg24);
            }
#endif
            rate_Mg24_He4_to_C12_O16_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Al27_p_to_n_Si27_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Al27);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Al27);
            }
#endif
            rate_Al27_p_to_n_Si27_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Al27_p_to_C12_O16_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Al27);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Al27);
            }
#endif
            rate_Al27_p_to_C12_O16_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Si27_He4_to_n_S30_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si27);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si27);
            }
#endif
            rate_Si27_He4_to_n_S30_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Si28_He4_to_p_P31_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si28);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si28);
            }
#endif
            rate_Si28_He4_to_p_P31_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Si28_He4_to_C12_Ne20_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si28);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si28);
            }
#endif
            rate_Si28_He4_to_C12_Ne20_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Si28_He4_to_O16_O16_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si28);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si28);
            }
#endif
            rate_Si28_He4_to_O16_O16_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Si30_p_to_He4_Al27_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Si30);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Si30);
            }
#endif
            rate_Si30_p_to_He4_Al27_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // P31_p_to_C12_Ne20_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_P31);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_P31);
            }
#endif
            rate_P31_p_to_C12_Ne20_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // P31_p_to_O16_O16_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_P31);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_P31);
            }
#endif
            rate_P31_p_to_O16_O16_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // S32_p_to_He4_P29_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_S32);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_S32);
            }
#endif
            rate_S32_p_to_He4_P29_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // S32_He4_to_p_Cl35_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_S32);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_S32);
            }
#endif
            rate_S32_He4_to_p_Cl35_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // S34_p_to_He4_P31_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_S34);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_S34);
            }
#endif
            rate_S34_p_to_He4_P31_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cl33_p_to_He4_S30_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cl33);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cl33);
            }
#endif
            rate_Cl33_p_to_He4_S30_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ar36_p_to_He4_Cl33_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ar36);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ar36);
            }
#endif
            rate_Ar36_p_to_He4_Cl33_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ar36_He4_to_p_K39_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ar36);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ar36);
            }
#endif
            rate_Ar36_He4_to_p_K39_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ar38_p_to_He4_Cl35_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ar38);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ar38);
            }
#endif
            rate_Ar38_p_to_He4_Cl35_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // K37_p_to_He4_Ar34_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_K37);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_K37);
            }
#endif
            rate_K37_p_to_He4_Ar34_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // K39_He4_to_p_Ca42_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_K39);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_K39);
            }
#endif
            rate_K39_He4_to_p_Ca42_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ca40_p_to_He4_K37_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ca40);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca40);
            }
#endif
            rate_Ca40_p_to_He4_K37_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ca40_He4_to_p_Sc43_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca40);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca40);
            }
#endif
            rate_Ca40_He4_to_p_Sc43_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Sc41_p_to_He4_Ca38_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Sc41);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc41);
            }
#endif
            rate_Sc41_p_to_He4_Ca38_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ti44_p_to_He4_Sc41_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti44);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti44);
            }
#endif
            rate_Ti44_p_to_He4_Sc41_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ti46_p_to_He4_Sc43_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti46);
            }
#endif
            rate_Ti46_p_to_He4_Sc43_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ti46_He4_to_p_V49_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti46);
            }
#endif
            rate_Ti46_He4_to_p_V49_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // V45_p_to_He4_Ti42_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_V45);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V45);
            }
#endif
            rate_V45_p_to_He4_Ti42_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // V47_p_to_He4_Ti44_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_V47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V47);
            }
#endif
            rate_V47_p_to_He4_Ti44_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cr48_p_to_He4_V45_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cr48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr48);
            }
#endif
            rate_Cr48_p_to_He4_V45_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cr50_p_to_He4_V47_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cr50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr50);
            }
#endif
            rate_Cr50_p_to_He4_V47_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cr50_He4_to_p_Mn53_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr50);
            }
#endif
            rate_Cr50_He4_to_p_Mn53_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Cr52_p_to_He4_V49_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cr52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr52);
            }
#endif
            rate_Cr52_p_to_He4_V49_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Mn49_p_to_He4_Cr46_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mn49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn49);
            }
#endif
            rate_Mn49_p_to_He4_Cr46_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Mn51_p_to_He4_Cr48_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mn51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn51);
            }
#endif
            rate_Mn51_p_to_He4_Cr48_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Mn53_He4_to_n_Co56_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn53);
            }
#endif
            rate_Mn53_He4_to_n_Co56_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Fe52_p_to_He4_Mn49_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe52);
            }
#endif
            rate_Fe52_p_to_He4_Mn49_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Fe54_p_to_He4_Mn51_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe54);
            }
#endif
            rate_Fe54_p_to_He4_Mn51_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Fe54_He4_to_p_Co57_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Fe54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe54);
            }
#endif
            rate_Fe54_He4_to_p_Co57_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Fe56_p_to_n_Co56_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe56);
            }
#endif
            rate_Fe56_p_to_n_Co56_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Fe56_p_to_He4_Mn53_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe56);
            }
#endif
            rate_Fe56_p_to_He4_Mn53_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Co53_p_to_He4_Fe50_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Co53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co53);
            }
#endif
            rate_Co53_p_to_He4_Fe50_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Co55_p_to_He4_Fe52_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Co55);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co55);
            }
#endif
            rate_Co55_p_to_He4_Fe52_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Co56_p_to_n_Ni56_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Co56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co56);
            }
#endif
            rate_Co56_p_to_n_Ni56_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Co56_He4_to_n_Cu59_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Co56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co56);
            }
#endif
            rate_Co56_He4_to_n_Cu59_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni56_p_to_He4_Co53_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni56);
            }
#endif
            rate_Ni56_p_to_He4_Co53_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni56_He4_to_p_Cu59_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ni56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ni56);
            }
#endif
            rate_Ni56_He4_to_p_Cu59_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni58_p_to_He4_Co55_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni58);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni58);
            }
#endif
            rate_Ni58_p_to_He4_Co55_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ni57_to_n_Ni56_derived
            rate_Ni57_to_n_Ni56_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ni58_to_n_Ni57_derived
            rate_Ni58_to_n_Ni57_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co54_to_n_Co53_derived
            rate_Co54_to_n_Co53_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Co55_to_n_Co54_derived
            rate_Co55_to_n_Co54_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe51_to_n_Fe50_derived
            rate_Fe51_to_n_Fe50_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe52_to_n_Fe51_derived
            rate_Fe52_to_n_Fe51_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe53_to_n_Fe52_derived
            rate_Fe53_to_n_Fe52_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe54_to_n_Fe53_derived
            rate_Fe54_to_n_Fe53_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe55_to_n_Fe54_derived
            rate_Fe55_to_n_Fe54_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Fe56_to_n_Fe55_derived
            rate_Fe56_to_n_Fe55_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn50_to_n_Mn49_derived
            rate_Mn50_to_n_Mn49_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn51_to_n_Mn50_derived
            rate_Mn51_to_n_Mn50_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn52_to_n_Mn51_derived
            rate_Mn52_to_n_Mn51_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Mn53_to_n_Mn52_derived
            rate_Mn53_to_n_Mn52_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr47_to_n_Cr46_derived
            rate_Cr47_to_n_Cr46_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr48_to_n_Cr47_derived
            rate_Cr48_to_n_Cr47_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr49_to_n_Cr48_derived
            rate_Cr49_to_n_Cr48_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr50_to_n_Cr49_derived
            rate_Cr50_to_n_Cr49_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr51_to_n_Cr50_derived
            rate_Cr51_to_n_Cr50_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Cr52_to_n_Cr51_derived
            rate_Cr52_to_n_Cr51_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V46_to_n_V45_derived
            rate_V46_to_n_V45_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V47_to_n_V46_derived
            rate_V47_to_n_V46_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V48_to_n_V47_derived
            rate_V48_to_n_V47_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // V49_to_n_V48_derived
            rate_V49_to_n_V48_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti43_to_n_Ti42_derived
            rate_Ti43_to_n_Ti42_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti44_to_n_Ti43_derived
            rate_Ti44_to_n_Ti43_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti45_to_n_Ti44_derived
            rate_Ti45_to_n_Ti44_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ti46_to_n_Ti45_derived
            rate_Ti46_to_n_Ti45_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc42_to_n_Sc41_derived
            rate_Sc42_to_n_Sc41_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Sc43_to_n_Sc42_derived
            rate_Sc43_to_n_Sc42_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca39_to_n_Ca38_derived
            rate_Ca39_to_n_Ca38_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca40_to_n_Ca39_derived
            rate_Ca40_to_n_Ca39_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca41_to_n_Ca40_derived
            rate_Ca41_to_n_Ca40_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ca42_to_n_Ca41_derived
            rate_Ca42_to_n_Ca41_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ar35_to_n_Ar34_derived
            rate_Ar35_to_n_Ar34_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ar36_to_n_Ar35_derived
            rate_Ar36_to_n_Ar35_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ar37_to_n_Ar36_derived
            rate_Ar37_to_n_Ar36_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ar38_to_n_Ar37_derived
            rate_Ar38_to_n_Ar37_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Si29_to_n_Si28_derived
            rate_Si29_to_n_Si28_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Si30_to_n_Si29_derived
            rate_Si30_to_n_Si29_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // S31_to_n_S30_derived
            rate_S31_to_n_S30_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // S32_to_n_S31_derived
            rate_S32_to_n_S31_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // S33_to_n_S32_derived
            rate_S33_to_n_S32_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // S34_to_n_S33_derived
            rate_S34_to_n_S33_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Na23_p_to_C12_C12_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Na23);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Na23);
            }
#endif
            rate_Na23_p_to_C12_C12_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Mg24_to_p_Na23_derived
            rate_Mg24_to_p_Na23_derived(rate_eval, tfactors, pf_cache);
        }

        {
            // Ne20_He4_to_C12_C12_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ne20);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ne20);
            }
#endif
            rate_Ne20_He4_to_C12_C12_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Ne20_He4_to_p_Na23_derived
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ne20);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ne20);
            }
#endif
            rate_Ne20_He4_to_p_Na23_derived(rate_eval, tfactors, log_scor, dlog_scor_dT, pf_cache);
        }

        {
            // Mg24_to_He4_Ne20_derived
            rate_Mg24_to_He4_Ne20_derived(rate_eval, tfactors, pf_cache);
        }


    }

}
#endif
