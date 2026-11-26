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
    void rate_B8_to_p_Be7_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // B8 --> p + Be7

        amrex::Real rate{}, drate_dT{};
        // Evaluate partition function terms

        amrex::Real p_log_pf, dp_log_pf_dT9;

        // setting p log(partition function) to 0.0 by default, independent of T
        p_log_pf = 0.0_rt;
        dp_log_pf_dT9 = 0.0_rt;

        amrex::Real Be7_log_pf, dBe7_log_pf_dT9;

        // setting Be7 log(partition function) to 0.0 by default, independent of T
        Be7_log_pf = 0.0_rt;
        dBe7_log_pf_dT9 = 0.0_rt;

        amrex::Real B8_log_pf, dB8_log_pf_dT9;

        // setting B8 log(partition function) to 0.0 by default, independent of T
        B8_log_pf = 0.0_rt;
        dB8_log_pf_dT9 = 0.0_rt;

        amrex::Real net_log_pf = p_log_pf + Be7_log_pf - B8_log_pf;
        [[maybe_unused]] amrex::Real net_dlog_pf_dT9 = dp_log_pf_dT9 + dBe7_log_pf_dT9 - dB8_log_pf_dT9;

        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ReacLib set derived from nacrr
        ln_set_rate =  31.03415329791634 + -8.927520483432566 * tfactors.T9i;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  8.927520483432566 * tfactors.T9i * tfactors.T9i;

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

        // ReacLib set derived from nacrn
        ln_set_rate =  35.83166329791634 + -1.5825204834325672 * tfactors.T9i + -10.264 * tfactors.T913i + -0.203472 * tfactors.T913
                             + 0.121083 * tfactors.T9 + -0.00700063 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        ln_set_rate += net_log_pf;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  1.5825204834325672 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -10.264 * tfactors.T943i + (1.0/3.0) * -0.203472 * tfactors.T923i
                                      + 0.121083 + (5.0/3.0) * -0.00700063 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_B8_to_p_Be7_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_B8_to_p_Be7_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_F17_to_p_O16_derived(T& rate_eval, const tf_t& tfactors, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // F17 --> p + O16

        amrex::Real rate{}, drate_dT{};
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
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9_rt;
            }
        }

        rate_eval.screened_rates(k_F17_to_p_O16_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_F17_to_p_O16_derived) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_F17_He4_to_p_Ne20_derived(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, [[maybe_unused]] part_fun::pf_cache_t& pf_cache) {

        // F17 + He4 --> p + Ne20

        amrex::Real rate{}, drate_dT{};
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
        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  0.004848799532844339 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -43.18 * tfactors.T943i + (1.0/3.0) * 4.46827 * tfactors.T923i
                                      + -1.63915 + (5.0/3.0) * 0.123483 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_F17_He4_to_p_Ne20_derived) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_F17_He4_to_p_Ne20_derived) = drate_dT;
        }
    }


}
#endif
