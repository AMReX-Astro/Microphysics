#ifndef REACLIB_RATES_H
#define REACLIB_RATES_H

#include <AMReX.H>
#include <AMReX_Print.H>

#include <tfactors.H>
#include <actual_network.H>
#include <rate_type.H>


namespace reaclib_rates {

    using namespace Rates;
    using namespace Species;

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_C12_He4_to_O16_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // C12 + He4 --> O16

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // nac2 
        ln_set_rate =  254.634 + -1.84097 * tfactors.T9i + 103.411 * tfactors.T913i + -420.567 * tfactors.T913
                             + 64.0874 * tfactors.T9 + -12.4624 * tfactors.T953 + 137.303 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  1.84097 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 103.411 * tfactors.T943i + (1.0/3.0) * -420.567 * tfactors.T923i
                                      + 64.0874 + (5.0/3.0) * -12.4624 * tfactors.T923 + 137.303 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // nac2 
        ln_set_rate =  69.6526 + -1.39254 * tfactors.T9i + 58.9128 * tfactors.T913i + -148.273 * tfactors.T913
                             + 9.08324 * tfactors.T9 + -0.541041 * tfactors.T953 + 70.3554 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  1.39254 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 58.9128 * tfactors.T943i + (1.0/3.0) * -148.273 * tfactors.T923i
                                      + 9.08324 + (5.0/3.0) * -0.541041 * tfactors.T923 + 70.3554 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_C12_He4_to_O16_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_C12_He4_to_O16_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_O16_He4_to_Ne20_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // O16 + He4 --> Ne20

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // co10r
        ln_set_rate =  3.88571 + -10.3585 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  10.3585 * tfactors.T9i * tfactors.T9i
                                      + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // co10n
        ln_set_rate =  23.903 + -39.7262 * tfactors.T913i + -0.210799 * tfactors.T913
                             + 0.442879 * tfactors.T9 + -0.0797753 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -39.7262 * tfactors.T943i + (1.0/3.0) * -0.210799 * tfactors.T923i
                                      + 0.442879 + (5.0/3.0) * -0.0797753 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // co10r
        ln_set_rate =  9.50848 + -12.7643 * tfactors.T9i + -3.65925 * tfactors.T913
                             + 0.714224 * tfactors.T9 + -0.00107508 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  12.7643 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -3.65925 * tfactors.T923i
                                      + 0.714224 + (5.0/3.0) * -0.00107508 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_O16_He4_to_Ne20_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_O16_He4_to_Ne20_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mg23_n_to_Mg24_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Mg23 + n --> Mg24

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  7.00613 + 2.66964 * tfactors.T913
                             + -0.448904 * tfactors.T9 + 0.0326505 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 2.66964 * tfactors.T923i
                                      + -0.448904 + (5.0/3.0) * 0.0326505 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Mg23_n_to_Mg24_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg23_n_to_Mg24_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mg23_He4_to_Si27_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Mg23 + He4 --> Si27

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  43.7008 + -53.203 * tfactors.T913i + -4.6318 * tfactors.T913
                             + -0.130951 * tfactors.T9 + 0.014691 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -53.203 * tfactors.T943i + (1.0/3.0) * -4.6318 * tfactors.T923i
                                      + -0.130951 + (5.0/3.0) * 0.014691 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Mg23_He4_to_Si27_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg23_He4_to_Si27_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mg24_He4_to_Si28_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Mg24 + He4 --> Si28

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // st08r
        ln_set_rate =  8.03977 + -15.629 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  15.629 * tfactors.T9i * tfactors.T9i
                                      + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // st08r
        ln_set_rate =  -50.5494 + -12.8332 * tfactors.T9i + 21.3721 * tfactors.T913i + 37.7649 * tfactors.T913
                             + -4.10635 * tfactors.T9 + 0.249618 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  12.8332 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 21.3721 * tfactors.T943i + (1.0/3.0) * 37.7649 * tfactors.T923i
                                      + -4.10635 + (5.0/3.0) * 0.249618 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Mg24_He4_to_Si28_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg24_He4_to_Si28_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Al27_p_to_Si28_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Al27 + p --> Si28

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10r
        ln_set_rate =  -13.6664 + -1.90396 * tfactors.T9i + 23.8634 * tfactors.T913
                             + -3.70135 * tfactors.T9 + 0.28964 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  1.90396 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 23.8634 * tfactors.T923i
                                      + -3.70135 + (5.0/3.0) * 0.28964 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10r
        ln_set_rate =  86.0234 + -0.387313 * tfactors.T9i + -26.8327 * tfactors.T913i + -116.137 * tfactors.T913
                             + 0.00950567 * tfactors.T9 + 0.00999755 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  0.387313 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -26.8327 * tfactors.T943i + (1.0/3.0) * -116.137 * tfactors.T923i
                                      + 0.00950567 + (5.0/3.0) * 0.00999755 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10n
        ln_set_rate =  21.1065 + -23.2205 * tfactors.T913i
                             + -2.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -23.2205 * tfactors.T943i
                                      + (5.0/3.0) * -2.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Al27_p_to_Si28_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Al27_p_to_Si28_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Al27_He4_to_P31_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Al27 + He4 --> P31

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  47.2333 + -56.5351 * tfactors.T913i + -0.896208 * tfactors.T913
                             + -1.72024 * tfactors.T9 + 0.185409 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -56.5351 * tfactors.T943i + (1.0/3.0) * -0.896208 * tfactors.T923i
                                      + -1.72024 + (5.0/3.0) * 0.185409 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Al27_He4_to_P31_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Al27_He4_to_P31_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si27_n_to_Si28_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Si27 + n --> Si28

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  9.14054 + 0.559804 * tfactors.T913
                             + -0.0233808 * tfactors.T9 + -0.00121313 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 0.559804 * tfactors.T923i
                                      + -0.0233808 + (5.0/3.0) * -0.00121313 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Si27_n_to_Si28_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si27_n_to_Si28_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si28_p_to_P29_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Si28 + p --> P29

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10n
        ln_set_rate =  16.1779 + -23.8173 * tfactors.T913i + 7.08203 * tfactors.T913
                             + -1.44753 * tfactors.T9 + 0.0804296 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -23.8173 * tfactors.T943i + (1.0/3.0) * 7.08203 * tfactors.T923i
                                      + -1.44753 + (5.0/3.0) * 0.0804296 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10r
        ln_set_rate =  5.73975 + -4.14232 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  4.14232 * tfactors.T9i * tfactors.T9i
                                      + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Si28_p_to_P29_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si28_p_to_P29_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si28_He4_to_S32_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Si28 + He4 --> S32

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  47.9212 + -59.4896 * tfactors.T913i + 4.47205 * tfactors.T913
                             + -4.78989 * tfactors.T9 + 0.557201 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -59.4896 * tfactors.T943i + (1.0/3.0) * 4.47205 * tfactors.T923i
                                      + -4.78989 + (5.0/3.0) * 0.557201 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Si28_He4_to_S32_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si28_He4_to_S32_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si30_p_to_P31_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Si30 + p --> P31

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // de20r
        ln_set_rate =  12.3695 + -6.64105 * tfactors.T9i
                             + -1.1191 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  6.64105 * tfactors.T9i * tfactors.T9i
                                      + -1.1191 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // de20r
        ln_set_rate =  8.79766 + -5.18231 * tfactors.T9i
                             + 1.2883 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  5.18231 * tfactors.T9i * tfactors.T9i
                                      + 1.2883 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // de20r
        ln_set_rate =  -334.266 + -1.13327 * tfactors.T9i
                             + -78.547 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  1.13327 * tfactors.T9i * tfactors.T9i
                                      + -78.547 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // de20r
        ln_set_rate =  -4.24208 + -1.25174 * tfactors.T9i
                             + 4.99034 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  1.25174 * tfactors.T9i * tfactors.T9i
                                      + 4.99034 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // de20r
        ln_set_rate =  -18.9198 + -0.945261 * tfactors.T9i
                             + 1.30331 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  0.945261 * tfactors.T9i * tfactors.T9i
                                      + 1.30331 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // de20r
        ln_set_rate =  -1138.34 + -95.8769 * tfactors.T9i
                             + -896.758 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  95.8769 * tfactors.T9i * tfactors.T9i
                                      + -896.758 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // de20r
        ln_set_rate =  9.96544 + -5.58963 * tfactors.T9i
                             + -1.57082 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  5.58963 * tfactors.T9i * tfactors.T9i
                                      + -1.57082 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // de20r
        ln_set_rate =  -34.8594 + -0.592934 * tfactors.T9i
                             + -1.62839 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  0.592934 * tfactors.T9i * tfactors.T9i
                                      + -1.62839 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Si30_p_to_P31_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si30_p_to_P31_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si30_He4_to_S34_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Si30 + He4 --> S34

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  47.5008 + -59.6559 * tfactors.T913i + 3.70141 * tfactors.T913
                             + -3.12537 * tfactors.T9 + 0.307626 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -59.6559 * tfactors.T943i + (1.0/3.0) * 3.70141 * tfactors.T923i
                                      + -3.12537 + (5.0/3.0) * 0.307626 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Si30_He4_to_S34_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si30_He4_to_S34_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P29_p_to_S30_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // P29 + p --> S30

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10r
        ln_set_rate =  -26.0 + -3.71938 * tfactors.T9i + 8.48834 * tfactors.T913i + 25.65 * tfactors.T913
                             + -3.79773 * tfactors.T9 + 0.320391 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  3.71938 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 8.48834 * tfactors.T943i + (1.0/3.0) * 25.65 * tfactors.T923i
                                      + -3.79773 + (5.0/3.0) * 0.320391 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10n
        ln_set_rate =  20.9731 + -25.6007 * tfactors.T913i
                             + -2.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -25.6007 * tfactors.T943i
                                      + (5.0/3.0) * -2.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_P29_p_to_S30_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P29_p_to_S30_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P29_He4_to_Cl33_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // P29 + He4 --> Cl33

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  41.9979 + -62.3802 * tfactors.T913i + 0.593062 * tfactors.T913
                             + -1.14226 * tfactors.T9 + 0.0934776 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -62.3802 * tfactors.T943i + (1.0/3.0) * 0.593062 * tfactors.T923i
                                      + -1.14226 + (5.0/3.0) * 0.0934776 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_P29_He4_to_Cl33_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P29_He4_to_Cl33_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P31_p_to_S32_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // P31 + p --> S32

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10r
        ln_set_rate =  0.821556 + -3.77704 * tfactors.T9i + 8.09341 * tfactors.T913
                             + -0.615971 * tfactors.T9 + 0.031159 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  3.77704 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 8.09341 * tfactors.T923i
                                      + -0.615971 + (5.0/3.0) * 0.031159 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10r
        ln_set_rate =  -2.66839 + -2.25958 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  2.25958 * tfactors.T9i * tfactors.T9i
                                      + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10n
        ln_set_rate =  19.2596 + -25.3278 * tfactors.T913i + 6.4931 * tfactors.T913
                             + -9.27513 * tfactors.T9 + -0.610439 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -25.3278 * tfactors.T943i + (1.0/3.0) * 6.4931 * tfactors.T923i
                                      + -9.27513 + (5.0/3.0) * -0.610439 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_P31_p_to_S32_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P31_p_to_S32_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P31_He4_to_Cl35_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // P31 + He4 --> Cl35

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  50.451 + -62.5433 * tfactors.T913i + -2.95026 * tfactors.T913
                             + -0.89652 * tfactors.T9 + 0.0774126 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -62.5433 * tfactors.T943i + (1.0/3.0) * -2.95026 * tfactors.T923i
                                      + -0.89652 + (5.0/3.0) * 0.0774126 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_P31_He4_to_Cl35_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P31_He4_to_Cl35_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S30_He4_to_Ar34_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // S30 + He4 --> Ar34

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  43.1127 + -65.211 * tfactors.T913i + -1.41447 * tfactors.T913
                             + -0.542976 * tfactors.T9 + 0.0211165 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -65.211 * tfactors.T943i + (1.0/3.0) * -1.41447 * tfactors.T923i
                                      + -0.542976 + (5.0/3.0) * 0.0211165 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_S30_He4_to_Ar34_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S30_He4_to_Ar34_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S32_p_to_Cl33_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // S32 + p --> Cl33

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10n
        ln_set_rate =  52.47 + -29.7741 * tfactors.T913i + -87.4473 * tfactors.T913
                             + 182.189 * tfactors.T9 + -128.625 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -29.7741 * tfactors.T943i + (1.0/3.0) * -87.4473 * tfactors.T923i
                                      + 182.189 + (5.0/3.0) * -128.625 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10r
        ln_set_rate =  69.3346 + -3.00371 * tfactors.T9i + -33.7204 * tfactors.T913i + -32.7355 * tfactors.T913
                             + 3.92526 * tfactors.T9 + -0.250479 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  3.00371 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -33.7204 * tfactors.T943i + (1.0/3.0) * -32.7355 * tfactors.T923i
                                      + 3.92526 + (5.0/3.0) * -0.250479 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10r
        ln_set_rate =  -27.2382 + -0.874107 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  0.874107 * tfactors.T9i * tfactors.T9i
                                      + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_S32_p_to_Cl33_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S32_p_to_Cl33_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S32_He4_to_Ar36_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // S32 + He4 --> Ar36

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  48.901 + -65.3709 * tfactors.T913i + 5.68294 * tfactors.T913
                             + -5.00388 * tfactors.T9 + 0.571407 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -65.3709 * tfactors.T943i + (1.0/3.0) * 5.68294 * tfactors.T923i
                                      + -5.00388 + (5.0/3.0) * 0.571407 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_S32_He4_to_Ar36_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S32_He4_to_Ar36_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S34_p_to_Cl35_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // S34 + p --> Cl35

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // se19r
        ln_set_rate =  -7.70046 + -2.02615 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  2.02615 * tfactors.T9i * tfactors.T9i
                                      + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // se19r
        ln_set_rate =  8.25443 + -5.75583 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  5.75583 * tfactors.T9i * tfactors.T9i
                                      + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // se19r
        ln_set_rate =  -14.8781 + -1.40414 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  1.40414 * tfactors.T9i * tfactors.T9i
                                      + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // se19r
        ln_set_rate =  7.00142 + -5.10598 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  5.10598 * tfactors.T9i * tfactors.T9i
                                      + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // se19r
        ln_set_rate =  -17.9366 + -1.17205 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  1.17205 * tfactors.T9i * tfactors.T9i
                                      + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // se19r
        ln_set_rate =  5.80604 + -4.6418 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  4.6418 * tfactors.T9i * tfactors.T9i
                                      + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // se19r
        ln_set_rate =  0.43449 + -3.16803 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  3.16803 * tfactors.T9i * tfactors.T9i
                                      + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // se19r
        ln_set_rate =  -60.0247 + -27.689 * tfactors.T9i + 454.53 * tfactors.T913i + -375.557 * tfactors.T913
                             + 12.6533 * tfactors.T9 + -0.408677 * tfactors.T953 + 245.118 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  27.689 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 454.53 * tfactors.T943i + (1.0/3.0) * -375.557 * tfactors.T923i
                                      + 12.6533 + (5.0/3.0) * -0.408677 * tfactors.T923 + 245.118 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_S34_p_to_Cl35_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S34_p_to_Cl35_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S34_He4_to_Ar38_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // S34 + He4 --> Ar38

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // cn18 
        ln_set_rate =  228.729 + -17.9675 * tfactors.T9i + 635.332 * tfactors.T913i + -902.05 * tfactors.T913
                             + 42.5623 * tfactors.T9 + -2.08434 * tfactors.T953 + 477.35 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  17.9675 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 635.332 * tfactors.T943i + (1.0/3.0) * -902.05 * tfactors.T923i
                                      + 42.5623 + (5.0/3.0) * -2.08434 * tfactors.T923 + 477.35 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_S34_He4_to_Ar38_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S34_He4_to_Ar38_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl33_p_to_Ar34_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cl33 + p --> Ar34

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  35.1297 + -27.8815 * tfactors.T913i + -3.18731 * tfactors.T913
                             + -1.76254 * tfactors.T9 + 0.264735 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -27.8815 * tfactors.T943i + (1.0/3.0) * -3.18731 * tfactors.T923i
                                      + -1.76254 + (5.0/3.0) * 0.264735 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Cl33_p_to_Ar34_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl33_p_to_Ar34_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl33_He4_to_K37_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cl33 + He4 --> K37

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  46.6859 + -68.1442 * tfactors.T913i + -2.72746 * tfactors.T913
                             + -0.317537 * tfactors.T9 + 0.0162782 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -68.1442 * tfactors.T943i + (1.0/3.0) * -2.72746 * tfactors.T923i
                                      + -0.317537 + (5.0/3.0) * 0.0162782 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Cl33_He4_to_K37_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl33_He4_to_K37_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl35_p_to_Ar36_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cl35 + p --> Ar36

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10r
        ln_set_rate =  -42.5249 + -0.564651 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  0.564651 * tfactors.T9i * tfactors.T9i
                                      + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10n
        ln_set_rate =  35.6868 + -27.8971 * tfactors.T913i + -16.2304 * tfactors.T913
                             + 35.255 * tfactors.T9 + -25.8411 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -27.8971 * tfactors.T943i + (1.0/3.0) * -16.2304 * tfactors.T923i
                                      + 35.255 + (5.0/3.0) * -25.8411 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10r
        ln_set_rate =  -7.84699 + -3.65092 * tfactors.T9i + 18.0179 * tfactors.T913
                             + -2.86304 * tfactors.T9 + 0.250854 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  3.65092 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 18.0179 * tfactors.T923i
                                      + -2.86304 + (5.0/3.0) * 0.250854 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10r
        ln_set_rate =  -9.03294 + -2.00996 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  2.00996 * tfactors.T9i * tfactors.T9i
                                      + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Cl35_p_to_Ar36_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl35_p_to_Ar36_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl35_He4_to_K39_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cl35 + He4 --> K39

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  52.718 + -68.2848 * tfactors.T913i + 0.0178545 * tfactors.T913
                             + -2.06783 * tfactors.T9 + 0.199659 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -68.2848 * tfactors.T943i + (1.0/3.0) * 0.0178545 * tfactors.T923i
                                      + -2.06783 + (5.0/3.0) * 0.199659 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Cl35_He4_to_K39_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl35_He4_to_K39_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar34_He4_to_Ca38_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ar34 + He4 --> Ca38

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  48.8813 + -70.8666 * tfactors.T913i + -5.98563 * tfactors.T913
                             + 0.582779 * tfactors.T9 + -0.093842 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -70.8666 * tfactors.T943i + (1.0/3.0) * -5.98563 * tfactors.T923i
                                      + 0.582779 + (5.0/3.0) * -0.093842 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ar34_He4_to_Ca38_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar34_He4_to_Ca38_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar36_p_to_K37_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ar36 + p --> K37

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10r
        ln_set_rate =  3.18278 + -8.98031 * tfactors.T9i + 7.50962 * tfactors.T913
                             + -1.68675 * tfactors.T9 + 0.189725 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  8.98031 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 7.50962 * tfactors.T923i
                                      + -1.68675 + (5.0/3.0) * 0.189725 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10r
        ln_set_rate =  4.69769 + -3.6251 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  3.6251 * tfactors.T9i * tfactors.T9i
                                      + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10n
        ln_set_rate =  21.4923 + -28.9682 * tfactors.T913i
                             + -1.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -28.9682 * tfactors.T943i
                                      + (5.0/3.0) * -1.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ar36_p_to_K37_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar36_p_to_K37_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar36_He4_to_Ca40_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ar36 + He4 --> Ca40

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  52.3486 + -71.0046 * tfactors.T913i + 4.0656 * tfactors.T913
                             + -5.26509 * tfactors.T9 + 0.683546 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -71.0046 * tfactors.T943i + (1.0/3.0) * 4.0656 * tfactors.T923i
                                      + -5.26509 + (5.0/3.0) * 0.683546 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ar36_He4_to_Ca40_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar36_He4_to_Ca40_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar38_p_to_K39_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ar38 + p --> K39

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  35.2834 + -29.0021 * tfactors.T913i + -0.525968 * tfactors.T913
                             + -1.94216 * tfactors.T9 + 0.267346 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -29.0021 * tfactors.T943i + (1.0/3.0) * -0.525968 * tfactors.T923i
                                      + -1.94216 + (5.0/3.0) * 0.267346 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ar38_p_to_K39_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar38_p_to_K39_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar38_He4_to_Ca42_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ar38 + He4 --> Ca42

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  52.6394 + -71.1296 * tfactors.T913i + 2.75299 * tfactors.T913
                             + -2.43406 * tfactors.T9 + 0.199511 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -71.1296 * tfactors.T943i + (1.0/3.0) * 2.75299 * tfactors.T923i
                                      + -2.43406 + (5.0/3.0) * 0.199511 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ar38_He4_to_Ca42_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar38_He4_to_Ca42_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K37_p_to_Ca38_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // K37 + p --> Ca38

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  36.2122 + -30.0596 * tfactors.T913i + -2.56332 * tfactors.T913
                             + -2.026 * tfactors.T9 + 0.289027 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -30.0596 * tfactors.T943i + (1.0/3.0) * -2.56332 * tfactors.T923i
                                      + -2.026 + (5.0/3.0) * 0.289027 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_K37_p_to_Ca38_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K37_p_to_Ca38_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K37_He4_to_Sc41_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // K37 + He4 --> Sc41

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  50.601 + -73.6776 * tfactors.T913i + -5.3733 * tfactors.T913
                             + 0.685776 * tfactors.T9 + -0.123735 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -73.6776 * tfactors.T943i + (1.0/3.0) * -5.3733 * tfactors.T923i
                                      + 0.685776 + (5.0/3.0) * -0.123735 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_K37_He4_to_Sc41_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K37_He4_to_Sc41_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K39_p_to_Ca40_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // K39 + p --> Ca40

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // lo18r
        ln_set_rate =  2761.38 + -5.22234 * tfactors.T9i + 802.18 * tfactors.T913i + -4010.27 * tfactors.T913
                             + 1136.19 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  5.22234 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 802.18 * tfactors.T943i + (1.0/3.0) * -4010.27 * tfactors.T923i
                                      + 1136.19 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // lo18r
        ln_set_rate =  588.099 + -12.5647 * tfactors.T9i + 641.844 * tfactors.T913i + -1248.49 * tfactors.T913
                             + 564.926 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  12.5647 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 641.844 * tfactors.T943i + (1.0/3.0) * -1248.49 * tfactors.T923i
                                      + 564.926 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // lo18r
        ln_set_rate =  102.252 + -1.66508 * tfactors.T9i + 41.1723 * tfactors.T913i + -149.299 * tfactors.T913
                             + 10.5229 * tfactors.T9 + -0.68208 * tfactors.T953 + 59.2367 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  1.66508 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 41.1723 * tfactors.T943i + (1.0/3.0) * -149.299 * tfactors.T923i
                                      + 10.5229 + (5.0/3.0) * -0.68208 * tfactors.T923 + 59.2367 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_K39_p_to_Ca40_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K39_p_to_Ca40_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K39_He4_to_Sc43_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // K39 + He4 --> Sc43

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  54.1202 + -73.8006 * tfactors.T913i + 1.87885 * tfactors.T913
                             + -2.75862 * tfactors.T9 + 0.279678 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -73.8006 * tfactors.T943i + (1.0/3.0) * 1.87885 * tfactors.T923i
                                      + -2.75862 + (5.0/3.0) * 0.279678 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_K39_He4_to_Sc43_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K39_He4_to_Sc43_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca38_He4_to_Ti42_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ca38 + He4 --> Ti42

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  58.4869 + -76.3063 * tfactors.T913i + -10.7213 * tfactors.T913
                             + 0.926161 * tfactors.T9 + -0.0855699 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -76.3063 * tfactors.T943i + (1.0/3.0) * -10.7213 * tfactors.T923i
                                      + 0.926161 + (5.0/3.0) * -0.0855699 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ca38_He4_to_Ti42_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca38_He4_to_Ti42_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca40_p_to_Sc41_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ca40 + p --> Sc41

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10n
        ln_set_rate =  20.0077 + -31.1711 * tfactors.T913i
                             + -1.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -31.1711 * tfactors.T943i
                                      + (5.0/3.0) * -1.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10r
        ln_set_rate =  50.687 + -33.5123 * tfactors.T9i + -33.7287 * tfactors.T913
                             + 5.22449 * tfactors.T9 + -0.357104 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  33.5123 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -33.7287 * tfactors.T923i
                                      + 5.22449 + (5.0/3.0) * -0.357104 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10r
        ln_set_rate =  5.62969 + -7.32236 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  7.32236 * tfactors.T9i * tfactors.T9i
                                      + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ca40_p_to_Sc41_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca40_p_to_Sc41_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca40_He4_to_Ti44_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ca40 + He4 --> Ti44

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // chw0 
        ln_set_rate =  53.75 + -76.4273 * tfactors.T913i + 3.87451 * tfactors.T913
                             + -3.61477 * tfactors.T9 + 0.367451 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -76.4273 * tfactors.T943i + (1.0/3.0) * 3.87451 * tfactors.T923i
                                      + -3.61477 + (5.0/3.0) * 0.367451 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ca40_He4_to_Ti44_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca40_He4_to_Ti44_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca42_p_to_Sc43_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ca42 + p --> Sc43

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  35.8192 + -31.1381 * tfactors.T913i + -0.151667 * tfactors.T913
                             + -1.76926 * tfactors.T9 + 0.219569 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -31.1381 * tfactors.T943i + (1.0/3.0) * -0.151667 * tfactors.T923i
                                      + -1.76926 + (5.0/3.0) * 0.219569 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ca42_p_to_Sc43_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca42_p_to_Sc43_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca42_He4_to_Ti46_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ca42 + He4 --> Ti46

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  54.833 + -76.5379 * tfactors.T913i + 2.80099 * tfactors.T913
                             + -2.34788 * tfactors.T9 + 0.150014 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -76.5379 * tfactors.T943i + (1.0/3.0) * 2.80099 * tfactors.T923i
                                      + -2.34788 + (5.0/3.0) * 0.150014 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ca42_He4_to_Ti46_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca42_He4_to_Ti46_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc41_p_to_Ti42_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Sc41 + p --> Ti42

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  38.3524 + -32.1615 * tfactors.T913i + -3.68811 * tfactors.T913
                             + -2.01234 * tfactors.T9 + 0.3232 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -32.1615 * tfactors.T943i + (1.0/3.0) * -3.68811 * tfactors.T923i
                                      + -2.01234 + (5.0/3.0) * 0.3232 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Sc41_p_to_Ti42_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc41_p_to_Ti42_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc41_He4_to_V45_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Sc41 + He4 --> V45

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  53.2525 + -79.013 * tfactors.T913i + -6.64623 * tfactors.T913
                             + 1.10793 * tfactors.T9 + -0.169427 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -79.013 * tfactors.T943i + (1.0/3.0) * -6.64623 * tfactors.T923i
                                      + 1.10793 + (5.0/3.0) * -0.169427 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Sc41_He4_to_V45_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc41_He4_to_V45_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc43_p_to_Ti44_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Sc43 + p --> Ti44

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  36.8432 + -32.1734 * tfactors.T913i + -1.77078 * tfactors.T913
                             + -2.21706 * tfactors.T9 + 0.298499 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -32.1734 * tfactors.T943i + (1.0/3.0) * -1.77078 * tfactors.T923i
                                      + -2.21706 + (5.0/3.0) * 0.298499 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Sc43_p_to_Ti44_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc43_p_to_Ti44_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc43_He4_to_V47_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Sc43 + He4 --> V47

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  59.0195 + -79.122 * tfactors.T913i + -7.07006 * tfactors.T913
                             + 0.424183 * tfactors.T9 + -0.0665231 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -79.122 * tfactors.T943i + (1.0/3.0) * -7.07006 * tfactors.T923i
                                      + 0.424183 + (5.0/3.0) * -0.0665231 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Sc43_He4_to_V47_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc43_He4_to_V47_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti42_He4_to_Cr46_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ti42 + He4 --> Cr46

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  54.5788 + -81.5598 * tfactors.T913i + -7.76887 * tfactors.T913
                             + 1.52445 * tfactors.T9 + -0.287125 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -81.5598 * tfactors.T943i + (1.0/3.0) * -7.76887 * tfactors.T923i
                                      + 1.52445 + (5.0/3.0) * -0.287125 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ti42_He4_to_Cr46_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti42_He4_to_Cr46_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti44_p_to_V45_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ti44 + p --> V45

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // nfisn
        ln_set_rate =  14.2807 + -1.40579 * tfactors.T9i + -22.1209 * tfactors.T913i + 2.87067 * tfactors.T913
                             + -0.0200651 * tfactors.T9 + -0.0194061 * tfactors.T953 + -6.51065 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  1.40579 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -22.1209 * tfactors.T943i + (1.0/3.0) * 2.87067 * tfactors.T923i
                                      + -0.0200651 + (5.0/3.0) * -0.0194061 * tfactors.T923 + -6.51065 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // nfisn
        ln_set_rate =  18.7183 + 0.159399 * tfactors.T9i + -41.059 * tfactors.T913i + 12.0361 * tfactors.T913
                             + -0.649888 * tfactors.T9 + 0.0175707 * tfactors.T953 + -6.46732 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  -0.159399 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -41.059 * tfactors.T943i + (1.0/3.0) * 12.0361 * tfactors.T923i
                                      + -0.649888 + (5.0/3.0) * 0.0175707 * tfactors.T923 + -6.46732 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // nfisn
        ln_set_rate =  21.7795 + -1.8856 * tfactors.T9i + -155.515 * tfactors.T913i + 135.903 * tfactors.T913
                             + -3.34812 * tfactors.T9 + -0.0139974 * tfactors.T953 + -87.4481 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  1.8856 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -155.515 * tfactors.T943i + (1.0/3.0) * 135.903 * tfactors.T923i
                                      + -3.34812 + (5.0/3.0) * -0.0139974 * tfactors.T923 + -87.4481 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // nfisn
        ln_set_rate =  -118.107 + -34.8499 * tfactors.T9i + 84.8583 * tfactors.T913i + 64.855 * tfactors.T913
                             + -6.59381 * tfactors.T9 + 0.413146 * tfactors.T953 + -1.42502 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  34.8499 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 84.8583 * tfactors.T943i + (1.0/3.0) * 64.855 * tfactors.T923i
                                      + -6.59381 + (5.0/3.0) * 0.413146 * tfactors.T923 + -1.42502 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // nfisn
        ln_set_rate =  2.92559 + -2.45693 * tfactors.T9i + -27.0951 * tfactors.T913i + 22.7374 * tfactors.T913
                             + -2.18213 * tfactors.T9 + 0.223757 * tfactors.T953 + -15.3241 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  2.45693 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -27.0951 * tfactors.T943i + (1.0/3.0) * 22.7374 * tfactors.T923i
                                      + -2.18213 + (5.0/3.0) * 0.223757 * tfactors.T923 + -15.3241 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ti44_p_to_V45_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti44_p_to_V45_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti44_He4_to_Cr48_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ti44 + He4 --> Cr48

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  64.7958 + -81.667 * tfactors.T913i + -10.6333 * tfactors.T913
                             + -0.672613 * tfactors.T9 + 0.161209 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -81.667 * tfactors.T943i + (1.0/3.0) * -10.6333 * tfactors.T923i
                                      + -0.672613 + (5.0/3.0) * 0.161209 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ti44_He4_to_Cr48_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti44_He4_to_Cr48_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti46_p_to_V47_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ti46 + p --> V47

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  35.7683 + -33.2034 * tfactors.T913i + 0.505619 * tfactors.T913
                             + -1.73654 * tfactors.T9 + 0.207342 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -33.2034 * tfactors.T943i + (1.0/3.0) * 0.505619 * tfactors.T923i
                                      + -1.73654 + (5.0/3.0) * 0.207342 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ti46_p_to_V47_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti46_p_to_V47_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti46_He4_to_Cr50_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ti46 + He4 --> Cr50

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  57.0372 + -81.7658 * tfactors.T913i + 3.48637 * tfactors.T913
                             + -3.9188 * tfactors.T9 + 0.440356 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -81.7658 * tfactors.T943i + (1.0/3.0) * 3.48637 * tfactors.T923i
                                      + -3.9188 + (5.0/3.0) * 0.440356 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ti46_He4_to_Cr50_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti46_He4_to_Cr50_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V45_p_to_Cr46_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // V45 + p --> Cr46

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  35.927 + -34.1967 * tfactors.T913i + -1.70253 * tfactors.T913
                             + -2.38375 * tfactors.T9 + 0.33421 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -34.1967 * tfactors.T943i + (1.0/3.0) * -1.70253 * tfactors.T923i
                                      + -2.38375 + (5.0/3.0) * 0.33421 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_V45_p_to_Cr46_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V45_p_to_Cr46_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V45_He4_to_Mn49_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // V45 + He4 --> Mn49

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  51.8178 + -84.1758 * tfactors.T913i + -2.09835 * tfactors.T913
                             + -0.662318 * tfactors.T9 + 0.0293536 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -84.1758 * tfactors.T943i + (1.0/3.0) * -2.09835 * tfactors.T923i
                                      + -0.662318 + (5.0/3.0) * 0.0293536 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_V45_He4_to_Mn49_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V45_He4_to_Mn49_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V47_p_to_Cr48_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // V47 + p --> Cr48

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // nfisn
        ln_set_rate =  42.6798 + -6.0593 * tfactors.T9i + -34.0548 * tfactors.T913i + -3.41973 * tfactors.T913
                             + 1.16501 * tfactors.T9 + -0.105543 * tfactors.T953 + -7.70886 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  6.0593 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -34.0548 * tfactors.T943i + (1.0/3.0) * -3.41973 * tfactors.T923i
                                      + 1.16501 + (5.0/3.0) * -0.105543 * tfactors.T923 + -7.70886 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // nfisn
        ln_set_rate =  511.463 + -5.29491 * tfactors.T9i + 317.171 * tfactors.T913i + -911.679 * tfactors.T913
                             + 94.4245 * tfactors.T9 + -10.1973 * tfactors.T953 + 330.727 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  5.29491 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 317.171 * tfactors.T943i + (1.0/3.0) * -911.679 * tfactors.T923i
                                      + 94.4245 + (5.0/3.0) * -10.1973 * tfactors.T923 + 330.727 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // nfisn
        ln_set_rate =  23.8315 + 0.246665 * tfactors.T9i + -45.9868 * tfactors.T913i + 13.6822 * tfactors.T913
                             + -0.376902 * tfactors.T9 + -0.0194875 * tfactors.T953 + -8.42325 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  -0.246665 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -45.9868 * tfactors.T943i + (1.0/3.0) * 13.6822 * tfactors.T923i
                                      + -0.376902 + (5.0/3.0) * -0.0194875 * tfactors.T923 + -8.42325 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // nfisn
        ln_set_rate =  40.5626 + -0.514414 * tfactors.T9i + -110.655 * tfactors.T913i + 83.0232 * tfactors.T913
                             + -19.7762 * tfactors.T9 + 3.03961 * tfactors.T953 + -49.4742 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  0.514414 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -110.655 * tfactors.T943i + (1.0/3.0) * 83.0232 * tfactors.T923i
                                      + -19.7762 + (5.0/3.0) * 3.03961 * tfactors.T923 + -49.4742 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_V47_p_to_Cr48_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V47_p_to_Cr48_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V47_He4_to_Mn51_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // V47 + He4 --> Mn51

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  56.8618 + -84.2732 * tfactors.T913i + -2.98061 * tfactors.T913
                             + -0.531361 * tfactors.T9 + 0.023612 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -84.2732 * tfactors.T943i + (1.0/3.0) * -2.98061 * tfactors.T923i
                                      + -0.531361 + (5.0/3.0) * 0.023612 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_V47_He4_to_Mn51_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V47_He4_to_Mn51_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V49_p_to_Cr50_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // V49 + p --> Cr50

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  35.2335 + -34.217 * tfactors.T913i + 1.28258 * tfactors.T913
                             + -1.5098 * tfactors.T9 + 0.142011 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -34.217 * tfactors.T943i + (1.0/3.0) * 1.28258 * tfactors.T923i
                                      + -1.5098 + (5.0/3.0) * 0.142011 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_V49_p_to_Cr50_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V49_p_to_Cr50_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V49_He4_to_Mn53_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // V49 + He4 --> Mn53

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  62.4555 + -84.3633 * tfactors.T913i + -5.28933 * tfactors.T913
                             + -1.3103 * tfactors.T9 + 0.187088 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -84.3633 * tfactors.T943i + (1.0/3.0) * -5.28933 * tfactors.T923i
                                      + -1.3103 + (5.0/3.0) * 0.187088 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_V49_He4_to_Mn53_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V49_He4_to_Mn53_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr46_He4_to_Fe50_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cr46 + He4 --> Fe50

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  56.3568 + -86.6499 * tfactors.T913i + -6.44312 * tfactors.T913
                             + 1.15026 * tfactors.T9 + -0.173961 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -86.6499 * tfactors.T943i + (1.0/3.0) * -6.44312 * tfactors.T923i
                                      + 1.15026 + (5.0/3.0) * -0.173961 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Cr46_He4_to_Fe50_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr46_He4_to_Fe50_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr48_p_to_Mn49_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cr48 + p --> Mn49

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  38.5336 + -35.1969 * tfactors.T913i + -3.36393 * tfactors.T913
                             + -1.55965 * tfactors.T9 + 0.237119 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -35.1969 * tfactors.T943i + (1.0/3.0) * -3.36393 * tfactors.T923i
                                      + -1.55965 + (5.0/3.0) * 0.237119 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Cr48_p_to_Mn49_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr48_p_to_Mn49_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr48_He4_to_Fe52_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cr48 + He4 --> Fe52

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  65.1754 + -86.7459 * tfactors.T913i + -9.79373 * tfactors.T913
                             + -0.772169 * tfactors.T9 + 0.155883 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -86.7459 * tfactors.T943i + (1.0/3.0) * -9.79373 * tfactors.T923i
                                      + -0.772169 + (5.0/3.0) * 0.155883 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Cr48_He4_to_Fe52_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr48_He4_to_Fe52_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr50_p_to_Mn51_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cr50 + p --> Mn51

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  35.8532 + -35.2065 * tfactors.T913i + 0.884621 * tfactors.T913
                             + -1.76421 * tfactors.T9 + 0.202494 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -35.2065 * tfactors.T943i + (1.0/3.0) * 0.884621 * tfactors.T923i
                                      + -1.76421 + (5.0/3.0) * 0.202494 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Cr50_p_to_Mn51_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr50_p_to_Mn51_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr50_He4_to_Fe54_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cr50 + He4 --> Fe54

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  59.7193 + -86.8349 * tfactors.T913i + 2.37343 * tfactors.T913
                             + -3.65677 * tfactors.T9 + 0.412128 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -86.8349 * tfactors.T943i + (1.0/3.0) * 2.37343 * tfactors.T923i
                                      + -3.65677 + (5.0/3.0) * 0.412128 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Cr50_He4_to_Fe54_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr50_He4_to_Fe54_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr52_p_to_Mn53_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cr52 + p --> Mn53

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  35.6954 + -35.2154 * tfactors.T913i + 0.73019 * tfactors.T913
                             + -1.26018 * tfactors.T9 + 0.121141 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -35.2154 * tfactors.T943i + (1.0/3.0) * 0.73019 * tfactors.T923i
                                      + -1.26018 + (5.0/3.0) * 0.121141 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Cr52_p_to_Mn53_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr52_p_to_Mn53_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr52_He4_to_Fe56_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cr52 + He4 --> Fe56

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  58.0146 + -86.9175 * tfactors.T913i + 5.655 * tfactors.T913
                             + -4.53625 * tfactors.T9 + 0.479806 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -86.9175 * tfactors.T943i + (1.0/3.0) * 5.655 * tfactors.T923i
                                      + -4.53625 + (5.0/3.0) * 0.479806 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Cr52_He4_to_Fe56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr52_He4_to_Fe56_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn49_p_to_Fe50_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Mn49 + p --> Fe50

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  36.0838 + -36.173 * tfactors.T913i + 1.28192 * tfactors.T913
                             + -3.04331 * tfactors.T9 + 0.393688 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -36.173 * tfactors.T943i + (1.0/3.0) * 1.28192 * tfactors.T923i
                                      + -3.04331 + (5.0/3.0) * 0.393688 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Mn49_p_to_Fe50_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn49_p_to_Fe50_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn49_He4_to_Co53_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Mn49 + He4 --> Co53

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  53.7894 + -89.1861 * tfactors.T913i + -0.631651 * tfactors.T913
                             + -0.826494 * tfactors.T9 + 0.0355619 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -89.1861 * tfactors.T943i + (1.0/3.0) * -0.631651 * tfactors.T923i
                                      + -0.826494 + (5.0/3.0) * 0.0355619 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Mn49_He4_to_Co53_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn49_He4_to_Co53_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn51_p_to_Fe52_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Mn51 + p --> Fe52

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  36.2596 + -36.1825 * tfactors.T913i + 0.873042 * tfactors.T913
                             + -2.89731 * tfactors.T9 + 0.364394 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -36.1825 * tfactors.T943i + (1.0/3.0) * 0.873042 * tfactors.T923i
                                      + -2.89731 + (5.0/3.0) * 0.364394 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Mn51_p_to_Fe52_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn51_p_to_Fe52_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn51_He4_to_Co55_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Mn51 + He4 --> Co55

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  65.9219 + -89.274 * tfactors.T913i + -10.4373 * tfactors.T913
                             + 1.00492 * tfactors.T9 + -0.125548 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -89.274 * tfactors.T943i + (1.0/3.0) * -10.4373 * tfactors.T923i
                                      + 1.00492 + (5.0/3.0) * -0.125548 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Mn51_He4_to_Co55_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn51_He4_to_Co55_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn53_p_to_Fe54_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Mn53 + p --> Fe54

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  35.586 + -36.1913 * tfactors.T913i + 1.44056 * tfactors.T913
                             + -1.50734 * tfactors.T9 + 0.142124 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -36.1913 * tfactors.T943i + (1.0/3.0) * 1.44056 * tfactors.T923i
                                      + -1.50734 + (5.0/3.0) * 0.142124 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Mn53_p_to_Fe54_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn53_p_to_Fe54_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn53_He4_to_Co57_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Mn53 + He4 --> Co57

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  62.0252 + -89.3558 * tfactors.T913i + 0.000882861 * tfactors.T913
                             + -2.79327 * tfactors.T9 + 0.309057 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -89.3558 * tfactors.T943i + (1.0/3.0) * 0.000882861 * tfactors.T923i
                                      + -2.79327 + (5.0/3.0) * 0.309057 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Mn53_He4_to_Co57_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn53_He4_to_Co57_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe52_p_to_Co53_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Fe52 + p --> Co53

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  39.2879 + -37.1457 * tfactors.T913i + -5.67417 * tfactors.T913
                             + -0.559644 * tfactors.T9 + 0.102768 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -37.1457 * tfactors.T943i + (1.0/3.0) * -5.67417 * tfactors.T923i
                                      + -0.559644 + (5.0/3.0) * 0.102768 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Fe52_p_to_Co53_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe52_p_to_Co53_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe52_He4_to_Ni56_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Fe52 + He4 --> Ni56

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  66.6417 + -91.6819 * tfactors.T913i + -9.51885 * tfactors.T913
                             + -0.533014 * tfactors.T9 + 0.0892607 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -91.6819 * tfactors.T943i + (1.0/3.0) * -9.51885 * tfactors.T923i
                                      + -0.533014 + (5.0/3.0) * 0.0892607 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Fe52_He4_to_Ni56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe52_He4_to_Ni56_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe54_p_to_Co55_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Fe54 + p --> Co55

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  36.2304 + -37.1544 * tfactors.T913i + 0.950364 * tfactors.T913
                             + -1.77529 * tfactors.T9 + 0.198562 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -37.1544 * tfactors.T943i + (1.0/3.0) * 0.950364 * tfactors.T923i
                                      + -1.77529 + (5.0/3.0) * 0.198562 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Fe54_p_to_Co55_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe54_p_to_Co55_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe54_He4_to_Ni58_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Fe54 + He4 --> Ni58

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  60.2478 + -91.7628 * tfactors.T913i + 4.23027 * tfactors.T913
                             + -3.31305 * tfactors.T9 + 0.271293 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -91.7628 * tfactors.T943i + (1.0/3.0) * 4.23027 * tfactors.T923i
                                      + -3.31305 + (5.0/3.0) * 0.271293 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Fe54_He4_to_Ni58_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe54_He4_to_Ni58_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe56_p_to_Co57_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Fe56 + p --> Co57

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  36.0665 + -37.1625 * tfactors.T913i + 1.06776 * tfactors.T913
                             + -1.31689 * tfactors.T9 + 0.122089 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -37.1625 * tfactors.T943i + (1.0/3.0) * 1.06776 * tfactors.T923i
                                      + -1.31689 + (5.0/3.0) * 0.122089 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Fe56_p_to_Co57_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe56_p_to_Co57_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co55_n_to_Co56_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Co55 + n --> Co56

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  17.065 + -1.86357 * tfactors.T913
                             + 0.616591 * tfactors.T9 + -0.0839313 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -1.86357 * tfactors.T923i
                                      + 0.616591 + (5.0/3.0) * -0.0839313 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Co55_n_to_Co56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co55_n_to_Co56_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co55_p_to_Ni56_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Co55 + p --> Ni56

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  37.3736 + -38.1053 * tfactors.T913i + -0.210947 * tfactors.T913
                             + -2.68377 * tfactors.T9 + 0.355814 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -38.1053 * tfactors.T943i + (1.0/3.0) * -0.210947 * tfactors.T923i
                                      + -2.68377 + (5.0/3.0) * 0.355814 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Co55_p_to_Ni56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co55_p_to_Ni56_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co55_He4_to_Cu59_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Co55 + He4 --> Cu59

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  65.2921 + -94.1404 * tfactors.T913i + -2.62786 * tfactors.T913
                             + -2.12066 * tfactors.T9 + 0.237999 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -94.1404 * tfactors.T943i + (1.0/3.0) * -2.62786 * tfactors.T923i
                                      + -2.12066 + (5.0/3.0) * 0.237999 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Co55_He4_to_Cu59_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co55_He4_to_Cu59_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co56_n_to_Co57_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Co56 + n --> Co57

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  17.3552 + -1.37855 * tfactors.T913
                             + 0.299896 * tfactors.T9 + -0.04382 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -1.37855 * tfactors.T923i
                                      + 0.299896 + (5.0/3.0) * -0.04382 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Co56_n_to_Co57_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co56_n_to_Co57_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co57_p_to_Ni58_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Co57 + p --> Ni58

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  36.0159 + -38.1133 * tfactors.T913i + 1.77414 * tfactors.T913
                             + -1.48268 * tfactors.T9 + 0.121073 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -38.1133 * tfactors.T943i + (1.0/3.0) * 1.77414 * tfactors.T923i
                                      + -1.48268 + (5.0/3.0) * 0.121073 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Co57_p_to_Ni58_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co57_p_to_Ni58_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni58_p_to_Cu59_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ni58 + p --> Cu59

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // si13n
        ln_set_rate =  36.6662 + -39.0526 * tfactors.T913i + 1.0436 * tfactors.T913
                             + -2.10834 * tfactors.T9 + 0.239559 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -39.0526 * tfactors.T943i + (1.0/3.0) * 1.0436 * tfactors.T923i
                                      + -2.10834 + (5.0/3.0) * 0.239559 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ni58_p_to_Cu59_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni58_p_to_Cu59_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_O16_C12_to_p_Al27_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // O16 + C12 --> p + Al27

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // cf88r
        ln_set_rate =  68.5253 + 0.205134 * tfactors.T9i + -119.242 * tfactors.T913i + 13.3667 * tfactors.T913
                             + 0.295425 * tfactors.T9 + -0.267288 * tfactors.T953 + -9.91729 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  -0.205134 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -119.242 * tfactors.T943i + (1.0/3.0) * 13.3667 * tfactors.T923i
                                      + 0.295425 + (5.0/3.0) * -0.267288 * tfactors.T923 + -9.91729 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_O16_C12_to_p_Al27_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_O16_C12_to_p_Al27_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_O16_C12_to_He4_Mg24_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // O16 + C12 --> He4 + Mg24

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // cf88r
        ln_set_rate =  48.5341 + 0.37204 * tfactors.T9i + -133.413 * tfactors.T913i + 50.1572 * tfactors.T913
                             + -3.15987 * tfactors.T9 + 0.0178251 * tfactors.T953 + -23.7027 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  -0.37204 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -133.413 * tfactors.T943i + (1.0/3.0) * 50.1572 * tfactors.T923i
                                      + -3.15987 + (5.0/3.0) * 0.0178251 * tfactors.T923 + -23.7027 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_O16_C12_to_He4_Mg24_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_O16_C12_to_He4_Mg24_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_O16_O16_to_p_P31_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // O16 + O16 --> p + P31

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // cf88r
        ln_set_rate =  85.2628 + 0.223453 * tfactors.T9i + -145.844 * tfactors.T913i + 8.72612 * tfactors.T913
                             + -0.554035 * tfactors.T9 + -0.137562 * tfactors.T953 + -6.88807 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  -0.223453 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -145.844 * tfactors.T943i + (1.0/3.0) * 8.72612 * tfactors.T923i
                                      + -0.554035 + (5.0/3.0) * -0.137562 * tfactors.T923 + -6.88807 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_O16_O16_to_p_P31_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_O16_O16_to_p_P31_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_O16_O16_to_He4_Si28_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // O16 + O16 --> He4 + Si28

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // cf88r
        ln_set_rate =  97.2435 + -0.268514 * tfactors.T9i + -119.324 * tfactors.T913i + -32.2497 * tfactors.T913
                             + 1.46214 * tfactors.T9 + -0.200893 * tfactors.T953 + 13.2148 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  0.268514 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -119.324 * tfactors.T943i + (1.0/3.0) * -32.2497 * tfactors.T923i
                                      + 1.46214 + (5.0/3.0) * -0.200893 * tfactors.T923 + 13.2148 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_O16_O16_to_He4_Si28_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_O16_O16_to_He4_Si28_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ne20_C12_to_p_P31_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ne20 + C12 --> p + P31

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // rolfr
        ln_set_rate =  -268.136 + -38.7624 * tfactors.T9i + 361.154 * tfactors.T913i + -92.643 * tfactors.T913
                             + -9.98738 * tfactors.T9 + 0.892737 * tfactors.T953 + 161.042 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  38.7624 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 361.154 * tfactors.T943i + (1.0/3.0) * -92.643 * tfactors.T923i
                                      + -9.98738 + (5.0/3.0) * 0.892737 * tfactors.T923 + 161.042 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ne20_C12_to_p_P31_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne20_C12_to_p_P31_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ne20_C12_to_He4_Si28_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ne20 + C12 --> He4 + Si28

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // rolfr
        ln_set_rate =  -308.905 + -47.2175 * tfactors.T9i + 514.197 * tfactors.T913i + -200.896 * tfactors.T913
                             + -6.42713 * tfactors.T9 + 0.758256 * tfactors.T953 + 236.359 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  47.2175 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 514.197 * tfactors.T943i + (1.0/3.0) * -200.896 * tfactors.T923i
                                      + -6.42713 + (5.0/3.0) * 0.758256 * tfactors.T923 + 236.359 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ne20_C12_to_He4_Si28_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne20_C12_to_He4_Si28_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mg23_n_to_He4_Ne20_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Mg23 + n --> He4 + Ne20

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  17.7448 + 1.83199 * tfactors.T913
                             + -0.290485 * tfactors.T9 + 0.0242929 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 1.83199 * tfactors.T923i
                                      + -0.290485 + (5.0/3.0) * 0.0242929 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Mg23_n_to_He4_Ne20_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg23_n_to_He4_Ne20_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mg23_n_to_C12_C12_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Mg23 + n --> C12 + C12

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // cf88r
        ln_set_rate =  -12.8267 + 11.4826 * tfactors.T913
                             + 1.82849 * tfactors.T9 + -0.34844 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 11.4826 * tfactors.T923i
                                      + 1.82849 + (5.0/3.0) * -0.34844 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Mg23_n_to_C12_C12_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg23_n_to_C12_C12_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Al27_p_to_He4_Mg24_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Al27 + p --> He4 + Mg24

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10r
        ln_set_rate =  -7.02789 + -4.2425 * tfactors.T9i + 18.0416 * tfactors.T913
                             + -1.54137 * tfactors.T9 + 0.0847506 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  4.2425 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 18.0416 * tfactors.T923i
                                      + -1.54137 + (5.0/3.0) * 0.0847506 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10r
        ln_set_rate =  -26.8683 + -0.963012 * tfactors.T9i + 5.18642 * tfactors.T913i + -34.7936 * tfactors.T913
                             + 168.225 * tfactors.T9 + -115.825 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  0.963012 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 5.18642 * tfactors.T943i + (1.0/3.0) * -34.7936 * tfactors.T923i
                                      + 168.225 + (5.0/3.0) * -115.825 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10n
        ln_set_rate =  29.4576 + -26.4162 * tfactors.T913i
                             + -2.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -26.4162 * tfactors.T943i
                                      + (5.0/3.0) * -2.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Al27_p_to_He4_Mg24_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Al27_p_to_He4_Mg24_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Al27_He4_to_p_Si30_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Al27 + He4 --> p + Si30

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  47.4856 + -56.5351 * tfactors.T913i + 1.60477 * tfactors.T913
                             + -1.40594 * tfactors.T9 + 0.127353 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -56.5351 * tfactors.T943i + (1.0/3.0) * 1.60477 * tfactors.T923i
                                      + -1.40594 + (5.0/3.0) * 0.127353 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Al27_He4_to_p_Si30_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Al27_He4_to_p_Si30_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si27_n_to_p_Al27_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Si27 + n --> p + Al27

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  17.84 + 2.10552 * tfactors.T913
                             + -0.318695 * tfactors.T9 + 0.0257496 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 2.10552 * tfactors.T923i
                                      + -0.318695 + (5.0/3.0) * 0.0257496 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Si27_n_to_p_Al27_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si27_n_to_p_Al27_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si27_n_to_He4_Mg24_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Si27 + n --> He4 + Mg24

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  18.0785 + 0.942981 * tfactors.T913
                             + -0.104624 * tfactors.T9 + 0.00723421 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 0.942981 * tfactors.T923i
                                      + -0.104624 + (5.0/3.0) * 0.00723421 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Si27_n_to_He4_Mg24_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si27_n_to_He4_Mg24_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si27_n_to_C12_O16_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Si27 + n --> C12 + O16

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // cf88r
        ln_set_rate =  -131.755 + 3.44391 * tfactors.T9i + -293.089 * tfactors.T913i + 414.404 * tfactors.T913
                             + -28.0562 * tfactors.T9 + 1.61807 * tfactors.T953 + -178.28 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  -3.44391 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -293.089 * tfactors.T943i + (1.0/3.0) * 414.404 * tfactors.T923i
                                      + -28.0562 + (5.0/3.0) * 1.61807 * tfactors.T923 + -178.28 * tfactors.T9i;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Si27_n_to_C12_O16_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si27_n_to_C12_O16_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P29_He4_to_p_S32_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // P29 + He4 --> p + S32

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  50.2824 + -62.3802 * tfactors.T913i + 0.459085 * tfactors.T913
                             + -0.870169 * tfactors.T9 + 0.0631143 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -62.3802 * tfactors.T943i + (1.0/3.0) * 0.459085 * tfactors.T923i
                                      + -0.870169 + (5.0/3.0) * 0.0631143 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_P29_He4_to_p_S32_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P29_He4_to_p_S32_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P31_p_to_He4_Si28_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // P31 + p --> He4 + Si28

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10r
        ln_set_rate =  -10.893 + -3.42575 * tfactors.T9i + 21.521 * tfactors.T913
                             + -1.90355 * tfactors.T9 + 0.092724 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  3.42575 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 21.521 * tfactors.T923i
                                      + -1.90355 + (5.0/3.0) * 0.092724 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10r
        ln_set_rate =  -12.919 + -1.87716 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  1.87716 * tfactors.T9i * tfactors.T9i
                                      + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10n
        ln_set_rate =  60.8829 + -31.932 * tfactors.T913i + -77.0334 * tfactors.T913
                             + -43.6847 * tfactors.T9 + -4.28955 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -31.932 * tfactors.T943i + (1.0/3.0) * -77.0334 * tfactors.T923i
                                      + -43.6847 + (5.0/3.0) * -4.28955 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_P31_p_to_He4_Si28_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P31_p_to_He4_Si28_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P31_He4_to_p_S34_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // P31 + He4 --> p + S34

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  49.0993 + -62.5433 * tfactors.T913i + 1.69339 * tfactors.T913
                             + -1.27872 * tfactors.T9 + 0.106616 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -62.5433 * tfactors.T943i + (1.0/3.0) * 1.69339 * tfactors.T923i
                                      + -1.27872 + (5.0/3.0) * 0.106616 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_P31_He4_to_p_S34_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P31_He4_to_p_S34_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S30_n_to_He4_Si27_reaclib(T& rate_eval, const tf_t& tfactors) {

        // S30 + n --> He4 + Si27

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  14.3423 + 1.55858 * tfactors.T913
                             + -0.0243677 * tfactors.T9 + 0.00835247 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 1.55858 * tfactors.T923i
                                      + -0.0243677 + (5.0/3.0) * 0.00835247 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_S30_n_to_He4_Si27_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S30_n_to_He4_Si27_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S30_He4_to_p_Cl33_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // S30 + He4 --> p + Cl33

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  50.7513 + -65.211 * tfactors.T913i + 0.403887 * tfactors.T913
                             + -0.808239 * tfactors.T9 + 0.0581461 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -65.211 * tfactors.T943i + (1.0/3.0) * 0.403887 * tfactors.T923i
                                      + -0.808239 + (5.0/3.0) * 0.0581461 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_S30_He4_to_p_Cl33_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S30_He4_to_p_Cl33_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl33_He4_to_p_Ar36_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cl33 + He4 --> p + Ar36

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  52.1588 + -68.1442 * tfactors.T913i + 0.291238 * tfactors.T913
                             + -0.791384 * tfactors.T9 + 0.0524823 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -68.1442 * tfactors.T943i + (1.0/3.0) * 0.291238 * tfactors.T923i
                                      + -0.791384 + (5.0/3.0) * 0.0524823 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Cl33_He4_to_p_Ar36_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl33_He4_to_p_Ar36_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl35_p_to_He4_S32_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cl35 + p --> He4 + S32

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10r
        ln_set_rate =  -1.01202 + -3.93495 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  3.93495 * tfactors.T9i * tfactors.T9i
                                      + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10r
        ln_set_rate =  -57.5294 + -0.532931 * tfactors.T9i + 25.5338 * tfactors.T913
                             + 6.45824 * tfactors.T9 + -0.950294 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  0.532931 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 25.5338 * tfactors.T923i
                                      + 6.45824 + (5.0/3.0) * -0.950294 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10n
        ln_set_rate =  32.12 + -30.9147 * tfactors.T913i + -1.2345 * tfactors.T913
                             + 22.5118 * tfactors.T9 + -33.0589 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -30.9147 * tfactors.T943i + (1.0/3.0) * -1.2345 * tfactors.T923i
                                      + 22.5118 + (5.0/3.0) * -33.0589 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10r
        ln_set_rate =  2.29121 + -6.00976 * tfactors.T9i + 5.33756 * tfactors.T913
                             + 1.64418 * tfactors.T9 + -0.246167 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  6.00976 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 5.33756 * tfactors.T923i
                                      + 1.64418 + (5.0/3.0) * -0.246167 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Cl35_p_to_He4_S32_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl35_p_to_He4_S32_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl35_He4_to_p_Ar38_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cl35 + He4 --> p + Ar38

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  51.1272 + -68.2848 * tfactors.T913i + 2.5993 * tfactors.T913
                             + -1.59144 * tfactors.T9 + 0.137745 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -68.2848 * tfactors.T943i + (1.0/3.0) * 2.5993 * tfactors.T923i
                                      + -1.59144 + (5.0/3.0) * 0.137745 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Cl35_He4_to_p_Ar38_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl35_He4_to_p_Ar38_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar34_He4_to_p_K37_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ar34 + He4 --> p + K37

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  53.1837 + -70.8666 * tfactors.T913i + 0.263069 * tfactors.T913
                             + -0.802985 * tfactors.T9 + 0.0550287 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -70.8666 * tfactors.T943i + (1.0/3.0) * 0.263069 * tfactors.T923i
                                      + -0.802985 + (5.0/3.0) * 0.0550287 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ar34_He4_to_p_K37_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar34_He4_to_p_K37_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K37_He4_to_p_Ca40_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // K37 + He4 --> p + Ca40

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  54.2966 + -73.6776 * tfactors.T913i + 0.493294 * tfactors.T913
                             + -0.902307 * tfactors.T9 + 0.0641391 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -73.6776 * tfactors.T943i + (1.0/3.0) * 0.493294 * tfactors.T923i
                                      + -0.902307 + (5.0/3.0) * 0.0641391 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_K37_He4_to_p_Ca40_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K37_He4_to_p_Ca40_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K39_p_to_He4_Ar36_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // K39 + p --> He4 + Ar36

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  20.5166 + -30.0732 * tfactors.T913i + 7.03263 * tfactors.T913
                             + -1.10085 * tfactors.T9 + 0.133768 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -30.0732 * tfactors.T943i + (1.0/3.0) * 7.03263 * tfactors.T923i
                                      + -1.10085 + (5.0/3.0) * 0.133768 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_K39_p_to_He4_Ar36_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K39_p_to_He4_Ar36_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca38_He4_to_p_Sc41_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ca38 + He4 --> p + Sc41

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  56.0864 + -76.3063 * tfactors.T913i + -0.741136 * tfactors.T913
                             + -0.511467 * tfactors.T9 + 0.0171615 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -76.3063 * tfactors.T943i + (1.0/3.0) * -0.741136 * tfactors.T923i
                                      + -0.511467 + (5.0/3.0) * 0.0171615 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ca38_He4_to_p_Sc41_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca38_He4_to_p_Sc41_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca42_p_to_He4_K39_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ca42 + p --> He4 + K39

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  -78.212 + -31.1381 * tfactors.T913i + 117.868 * tfactors.T913
                             + -27.0455 * tfactors.T9 + 2.75661 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -31.1381 * tfactors.T943i + (1.0/3.0) * 117.868 * tfactors.T923i
                                      + -27.0455 + (5.0/3.0) * 2.75661 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ca42_p_to_He4_K39_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca42_p_to_He4_K39_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc41_He4_to_p_Ti44_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Sc41 + He4 --> p + Ti44

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  57.4853 + -79.013 * tfactors.T913i + -0.989013 * tfactors.T913
                             + -0.463943 * tfactors.T9 + 0.00902542 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -79.013 * tfactors.T943i + (1.0/3.0) * -0.989013 * tfactors.T923i
                                      + -0.463943 + (5.0/3.0) * 0.00902542 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Sc41_He4_to_p_Ti44_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc41_He4_to_p_Ti44_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc43_p_to_He4_Ca40_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Sc43 + p --> He4 + Ca40

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  34.8559 + -32.1734 * tfactors.T913i + 0.0296879 * tfactors.T913
                             + -0.95232 * tfactors.T9 + 0.129022 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -32.1734 * tfactors.T943i + (1.0/3.0) * 0.0296879 * tfactors.T923i
                                      + -0.95232 + (5.0/3.0) * 0.129022 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Sc43_p_to_He4_Ca40_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc43_p_to_He4_Ca40_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc43_He4_to_p_Ti46_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Sc43 + He4 --> p + Ti46

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  57.3723 + -79.122 * tfactors.T913i + -0.495539 * tfactors.T913
                             + -0.618414 * tfactors.T9 + 0.0175873 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -79.122 * tfactors.T943i + (1.0/3.0) * -0.495539 * tfactors.T923i
                                      + -0.618414 + (5.0/3.0) * 0.0175873 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Sc43_He4_to_p_Ti46_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc43_He4_to_p_Ti46_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti42_He4_to_p_V45_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ti42 + He4 --> p + V45

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  58.4447 + -81.5598 * tfactors.T913i + -1.03231 * tfactors.T913
                             + -0.471149 * tfactors.T9 + 0.0112292 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -81.5598 * tfactors.T943i + (1.0/3.0) * -1.03231 * tfactors.T923i
                                      + -0.471149 + (5.0/3.0) * 0.0112292 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ti42_He4_to_p_V45_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti42_He4_to_p_V45_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti44_He4_to_p_V47_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ti44 + He4 --> p + V47

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // chw0r
        ln_set_rate =  -76.5154 + -10.7931 * tfactors.T9i + 70.2835 * tfactors.T913
                             + -7.99061 * tfactors.T9 + 0.486213 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  10.7931 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 70.2835 * tfactors.T923i
                                      + -7.99061 + (5.0/3.0) * 0.486213 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ti44_He4_to_p_V47_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti44_He4_to_p_V47_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V45_He4_to_p_Cr48_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // V45 + He4 --> p + Cr48

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  59.6809 + -84.1758 * tfactors.T913i + -1.12403 * tfactors.T913
                             + -0.469581 * tfactors.T9 + 0.00932056 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -84.1758 * tfactors.T943i + (1.0/3.0) * -1.12403 * tfactors.T923i
                                      + -0.469581 + (5.0/3.0) * 0.00932056 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_V45_He4_to_p_Cr48_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V45_He4_to_p_Cr48_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V47_He4_to_p_Cr50_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // V47 + He4 --> p + Cr50

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  59.9106 + -84.2732 * tfactors.T913i + -1.10027 * tfactors.T913
                             + -0.487084 * tfactors.T9 + 0.00453269 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -84.2732 * tfactors.T943i + (1.0/3.0) * -1.10027 * tfactors.T923i
                                      + -0.487084 + (5.0/3.0) * 0.00453269 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_V47_He4_to_p_Cr50_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V47_He4_to_p_Cr50_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V49_p_to_He4_Ti46_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // V49 + p --> He4 + Ti46

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  -20.0996 + -34.217 * tfactors.T913i + 47.2584 * tfactors.T913
                             + -9.18871 * tfactors.T9 + 0.87823 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -34.217 * tfactors.T943i + (1.0/3.0) * 47.2584 * tfactors.T923i
                                      + -9.18871 + (5.0/3.0) * 0.87823 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_V49_p_to_He4_Ti46_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V49_p_to_He4_Ti46_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V49_He4_to_p_Cr52_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // V49 + He4 --> p + Cr52

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  57.428 + -84.3633 * tfactors.T913i + 2.73419 * tfactors.T913
                             + -2.05655 * tfactors.T9 + 0.203834 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -84.3633 * tfactors.T943i + (1.0/3.0) * 2.73419 * tfactors.T923i
                                      + -2.05655 + (5.0/3.0) * 0.203834 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_V49_He4_to_p_Cr52_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V49_He4_to_p_Cr52_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr46_He4_to_p_Mn49_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cr46 + He4 --> p + Mn49

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  60.5689 + -86.6499 * tfactors.T913i + -1.08241 * tfactors.T913
                             + -0.516569 * tfactors.T9 + 0.0131721 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -86.6499 * tfactors.T943i + (1.0/3.0) * -1.08241 * tfactors.T923i
                                      + -0.516569 + (5.0/3.0) * 0.0131721 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Cr46_He4_to_p_Mn49_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr46_He4_to_p_Mn49_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr48_He4_to_p_Mn51_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cr48 + He4 --> p + Mn51

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  59.2276 + -86.7459 * tfactors.T913i + 1.05653 * tfactors.T913
                             + -1.15757 * tfactors.T9 + 0.0877546 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -86.7459 * tfactors.T943i + (1.0/3.0) * 1.05653 * tfactors.T923i
                                      + -1.15757 + (5.0/3.0) * 0.0877546 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Cr48_He4_to_p_Mn51_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr48_He4_to_p_Mn51_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn49_He4_to_p_Fe52_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Mn49 + He4 --> p + Fe52

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  62.437 + -89.1861 * tfactors.T913i + -1.31711 * tfactors.T913
                             + -0.507338 * tfactors.T9 + 0.0108137 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -89.1861 * tfactors.T943i + (1.0/3.0) * -1.31711 * tfactors.T923i
                                      + -0.507338 + (5.0/3.0) * 0.0108137 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Mn49_He4_to_p_Fe52_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn49_He4_to_p_Fe52_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn51_He4_to_p_Fe54_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Mn51 + He4 --> p + Fe54

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  62.2777 + -89.274 * tfactors.T913i + -0.862452 * tfactors.T913
                             + -0.635672 * tfactors.T9 + 0.0196464 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -89.274 * tfactors.T943i + (1.0/3.0) * -0.862452 * tfactors.T923i
                                      + -0.635672 + (5.0/3.0) * 0.0196464 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Mn51_He4_to_p_Fe54_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn51_He4_to_p_Fe54_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn53_p_to_He4_Cr50_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Mn53 + p --> He4 + Cr50

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  -68.3424 + -36.1913 * tfactors.T913i + 99.0857 * tfactors.T913
                             + -20.8943 * tfactors.T9 + 2.06209 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -36.1913 * tfactors.T943i + (1.0/3.0) * 99.0857 * tfactors.T923i
                                      + -20.8943 + (5.0/3.0) * 2.06209 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Mn53_p_to_He4_Cr50_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn53_p_to_He4_Cr50_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn53_He4_to_p_Fe56_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Mn53 + He4 --> p + Fe56

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  54.0176 + -89.3558 * tfactors.T913i + 9.25287 * tfactors.T913
                             + -3.70352 * tfactors.T9 + 0.372202 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -89.3558 * tfactors.T943i + (1.0/3.0) * 9.25287 * tfactors.T923i
                                      + -3.70352 + (5.0/3.0) * 0.372202 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Mn53_He4_to_p_Fe56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn53_He4_to_p_Fe56_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe50_He4_to_p_Co53_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Fe50 + He4 --> p + Co53

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  62.7305 + -91.5952 * tfactors.T913i + -1.33713 * tfactors.T913
                             + -0.465416 * tfactors.T9 + 0.00761392 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -91.5952 * tfactors.T943i + (1.0/3.0) * -1.33713 * tfactors.T923i
                                      + -0.465416 + (5.0/3.0) * 0.00761392 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Fe50_He4_to_p_Co53_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe50_He4_to_p_Co53_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe52_He4_to_p_Co55_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Fe52 + He4 --> p + Co55

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  62.2207 + -91.6819 * tfactors.T913i + -0.329235 * tfactors.T913
                             + -0.780924 * tfactors.T9 + 0.0425179 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -91.6819 * tfactors.T943i + (1.0/3.0) * -0.329235 * tfactors.T923i
                                      + -0.780924 + (5.0/3.0) * 0.0425179 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Fe52_He4_to_p_Co55_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe52_He4_to_p_Co55_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co53_He4_to_p_Ni56_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Co53 + He4 --> p + Ni56

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  64.2463 + -94.0605 * tfactors.T913i + -1.48203 * tfactors.T913
                             + -0.476442 * tfactors.T9 + 0.00631587 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -94.0605 * tfactors.T943i + (1.0/3.0) * -1.48203 * tfactors.T923i
                                      + -0.476442 + (5.0/3.0) * 0.00631587 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Co53_He4_to_p_Ni56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co53_He4_to_p_Ni56_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co55_He4_to_p_Ni58_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Co55 + He4 --> p + Ni58

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  60.2281 + -94.1404 * tfactors.T913i + 3.39179 * tfactors.T913
                             + -1.71062 * tfactors.T9 + 0.133003 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -94.1404 * tfactors.T943i + (1.0/3.0) * 3.39179 * tfactors.T923i
                                      + -1.71062 + (5.0/3.0) * 0.133003 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Co55_He4_to_p_Ni58_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co55_He4_to_p_Ni58_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co56_n_to_p_Fe56_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Co56 + n --> p + Fe56

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  20.4539 + -1.13331 * tfactors.T913
                             + 0.347185 * tfactors.T9 + -0.0328879 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -1.13331 * tfactors.T923i
                                      + 0.347185 + (5.0/3.0) * -0.0328879 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Co56_n_to_p_Fe56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co56_n_to_p_Fe56_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co56_n_to_He4_Mn53_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Co56 + n --> He4 + Mn53

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  12.5672 + -2.81854 * tfactors.T913
                             + 1.84188 * tfactors.T9 + -0.164509 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -2.81854 * tfactors.T923i
                                      + 1.84188 + (5.0/3.0) * -0.164509 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Co56_n_to_He4_Mn53_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co56_n_to_He4_Mn53_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co57_p_to_He4_Fe54_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Co57 + p --> He4 + Fe54

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  -2.1962 + -38.1133 * tfactors.T913i + 29.3541 * tfactors.T913
                             + -4.75966 * tfactors.T9 + 0.40418 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -38.1133 * tfactors.T943i + (1.0/3.0) * 29.3541 * tfactors.T923i
                                      + -4.75966 + (5.0/3.0) * 0.40418 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Co57_p_to_He4_Fe54_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co57_p_to_He4_Fe54_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni56_n_to_p_Co56_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ni56 + n --> p + Co56

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  15.5693 + 1.76846 * tfactors.T913
                             + 0.197992 * tfactors.T9 + -0.017494 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 1.76846 * tfactors.T923i
                                      + 0.197992 + (5.0/3.0) * -0.017494 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ni56_n_to_p_Co56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni56_n_to_p_Co56_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu59_n_to_He4_Co56_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cu59 + n --> He4 + Co56

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  11.0108 + -1.00882 * tfactors.T913
                             + 1.88993 * tfactors.T9 + -0.202319 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -1.00882 * tfactors.T923i
                                      + 1.88993 + (5.0/3.0) * -0.202319 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Cu59_n_to_He4_Co56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu59_n_to_He4_Co56_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu59_p_to_He4_Ni56_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cu59 + p --> He4 + Ni56

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  20.7688 + -39.9808 * tfactors.T913i + 13.6773 * tfactors.T913
                             + -3.76429 * tfactors.T9 + 0.438096 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -39.9808 * tfactors.T943i + (1.0/3.0) * 13.6773 * tfactors.T923i
                                      + -3.76429 + (5.0/3.0) * 0.438096 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Cu59_p_to_He4_Ni56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu59_p_to_He4_Ni56_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_He4_He4_He4_to_C12_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // 3 He4 --> C12

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // fy05r
        ln_set_rate =  -11.7884 + -1.02446 * tfactors.T9i + -23.57 * tfactors.T913i + 20.4886 * tfactors.T913
                             + -12.9882 * tfactors.T9 + -20.0 * tfactors.T953 + -2.16667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  1.02446 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -23.57 * tfactors.T943i + (1.0/3.0) * 20.4886 * tfactors.T923i
                                      + -12.9882 + (5.0/3.0) * -20.0 * tfactors.T923 + -2.16667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // fy05n
        ln_set_rate =  -0.971052 + -37.06 * tfactors.T913i + 29.3493 * tfactors.T913
                             + -115.507 * tfactors.T9 + -10.0 * tfactors.T953 + -1.33333 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -37.06 * tfactors.T943i + (1.0/3.0) * 29.3493 * tfactors.T923i
                                      + -115.507 + (5.0/3.0) * -10.0 * tfactors.T923 + -1.33333 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // fy05r
        ln_set_rate =  -24.3505 + -4.12656 * tfactors.T9i + -13.49 * tfactors.T913i + 21.4259 * tfactors.T913
                             + -1.34769 * tfactors.T9 + 0.0879816 * tfactors.T953 + -13.1653 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  4.12656 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -13.49 * tfactors.T943i + (1.0/3.0) * 21.4259 * tfactors.T923i
                                      + -1.34769 + (5.0/3.0) * 0.0879816 * tfactors.T923 + -13.1653 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_He4_He4_He4_to_C12_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He4_He4_He4_to_C12_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni56_n_to_Ni57_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ni56 + n --> Ni57

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  16.0765 + -1.19665 * tfactors.T913
                             + 0.507179 * tfactors.T9 + -0.074604 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -1.19665 * tfactors.T923i
                                      + 0.507179 + (5.0/3.0) * -0.074604 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ni56_n_to_Ni57_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni56_n_to_Ni57_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni57_n_to_Ni58_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ni57 + n --> Ni58

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  17.2731 + -1.90814 * tfactors.T913
                             + 0.493188 * tfactors.T9 + -0.0684633 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -1.90814 * tfactors.T923i
                                      + 0.493188 + (5.0/3.0) * -0.0684633 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ni57_n_to_Ni58_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni57_n_to_Ni58_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co53_n_to_Co54_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Co53 + n --> Co54

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  11.2037 + -1.14795 * tfactors.T913
                             + 0.413083 * tfactors.T9 + -0.0515627 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -1.14795 * tfactors.T923i
                                      + 0.413083 + (5.0/3.0) * -0.0515627 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Co53_n_to_Co54_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co53_n_to_Co54_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co54_n_to_Co55_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Co54 + n --> Co55

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  14.3761 + 0.433187 * tfactors.T913
                             + 0.0888642 * tfactors.T9 + -0.0315521 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 0.433187 * tfactors.T923i
                                      + 0.0888642 + (5.0/3.0) * -0.0315521 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Co54_n_to_Co55_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co54_n_to_Co55_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe50_n_to_Fe51_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Fe50 + n --> Fe51

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  13.2307 + -1.43737 * tfactors.T913
                             + 0.511909 * tfactors.T9 + -0.0531273 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -1.43737 * tfactors.T923i
                                      + 0.511909 + (5.0/3.0) * -0.0531273 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Fe50_n_to_Fe51_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe50_n_to_Fe51_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe51_n_to_Fe52_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Fe51 + n --> Fe52

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  11.0052 + -0.566649 * tfactors.T913
                             + 0.19001 * tfactors.T9 + -0.0219395 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -0.566649 * tfactors.T923i
                                      + 0.19001 + (5.0/3.0) * -0.0219395 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Fe51_n_to_Fe52_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe51_n_to_Fe52_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe52_n_to_Fe53_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Fe52 + n --> Fe53

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  15.8885 + -0.344319 * tfactors.T913
                             + 0.178277 * tfactors.T9 + -0.0334326 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -0.344319 * tfactors.T923i
                                      + 0.178277 + (5.0/3.0) * -0.0334326 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Fe52_n_to_Fe53_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe52_n_to_Fe53_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe53_n_to_Fe54_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Fe53 + n --> Fe54

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  16.4534 + -1.10421 * tfactors.T913
                             + 0.379905 * tfactors.T9 + -0.0581878 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -1.10421 * tfactors.T923i
                                      + 0.379905 + (5.0/3.0) * -0.0581878 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Fe53_n_to_Fe54_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe53_n_to_Fe54_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe54_n_to_Fe55_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Fe54 + n --> Fe55

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  -0.80864 + 0.0591716 * tfactors.T9i + -8.66617 * tfactors.T913i + 26.4472 * tfactors.T913
                             + -1.9222 * tfactors.T9 + 0.0986404 * tfactors.T953 + -9.78317 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  -0.0591716 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -8.66617 * tfactors.T943i + (1.0/3.0) * 26.4472 * tfactors.T923i
                                      + -1.9222 + (5.0/3.0) * 0.0986404 * tfactors.T923 + -9.78317 * tfactors.T9i;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Fe54_n_to_Fe55_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe54_n_to_Fe55_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe55_n_to_Fe56_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Fe55 + n --> Fe56

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  21.7202 + -0.0955677 * tfactors.T9i + 8.06062 * tfactors.T913i + -14.4809 * tfactors.T913
                             + 0.94252 * tfactors.T9 + -0.0776007 * tfactors.T953 + 6.47093 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  0.0955677 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 8.06062 * tfactors.T943i + (1.0/3.0) * -14.4809 * tfactors.T923i
                                      + 0.94252 + (5.0/3.0) * -0.0776007 * tfactors.T923 + 6.47093 * tfactors.T9i;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Fe55_n_to_Fe56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe55_n_to_Fe56_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn49_n_to_Mn50_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Mn49 + n --> Mn50

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  11.269 + -0.759512 * tfactors.T913
                             + 0.246336 * tfactors.T9 + -0.0261234 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -0.759512 * tfactors.T923i
                                      + 0.246336 + (5.0/3.0) * -0.0261234 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Mn49_n_to_Mn50_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn49_n_to_Mn50_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn50_n_to_Mn51_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Mn50 + n --> Mn51

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  15.4434 + -0.538879 * tfactors.T913
                             + 0.284528 * tfactors.T9 + -0.0459849 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -0.538879 * tfactors.T923i
                                      + 0.284528 + (5.0/3.0) * -0.0459849 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Mn50_n_to_Mn51_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn50_n_to_Mn51_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn51_n_to_Mn52_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Mn51 + n --> Mn52

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  17.2055 + -1.83611 * tfactors.T913
                             + 0.460384 * tfactors.T9 + -0.0584947 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -1.83611 * tfactors.T923i
                                      + 0.460384 + (5.0/3.0) * -0.0584947 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Mn51_n_to_Mn52_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn51_n_to_Mn52_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn52_n_to_Mn53_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Mn52 + n --> Mn53

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  17.0941 + -0.91084 * tfactors.T913
                             + 0.162511 * tfactors.T9 + -0.0279619 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -0.91084 * tfactors.T923i
                                      + 0.162511 + (5.0/3.0) * -0.0279619 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Mn52_n_to_Mn53_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn52_n_to_Mn53_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr46_n_to_Cr47_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cr46 + n --> Cr47

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  15.4822 + -1.5945 * tfactors.T913
                             + 0.39528 * tfactors.T9 + -0.0372417 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -1.5945 * tfactors.T923i
                                      + 0.39528 + (5.0/3.0) * -0.0372417 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Cr46_n_to_Cr47_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr46_n_to_Cr47_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr47_n_to_Cr48_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cr47 + n --> Cr48

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  12.5822 + -1.36332 * tfactors.T913
                             + 0.288027 * tfactors.T9 + -0.0333017 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -1.36332 * tfactors.T923i
                                      + 0.288027 + (5.0/3.0) * -0.0333017 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Cr47_n_to_Cr48_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr47_n_to_Cr48_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr48_n_to_Cr49_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cr48 + n --> Cr49

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  15.8544 + -0.000452357 * tfactors.T913
                             + 0.0412838 * tfactors.T9 + -0.0160865 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -0.000452357 * tfactors.T923i
                                      + 0.0412838 + (5.0/3.0) * -0.0160865 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Cr48_n_to_Cr49_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr48_n_to_Cr49_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr49_n_to_Cr50_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cr49 + n --> Cr50

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  17.8355 + -1.98927 * tfactors.T913
                             + 0.465124 * tfactors.T9 + -0.0582249 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -1.98927 * tfactors.T923i
                                      + 0.465124 + (5.0/3.0) * -0.0582249 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Cr49_n_to_Cr50_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr49_n_to_Cr50_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr50_n_to_Cr51_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cr50 + n --> Cr51

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  15.1066 + -0.0225613 * tfactors.T9i + 0.815037 * tfactors.T913i + -0.708732 * tfactors.T913
                             + 0.386295 * tfactors.T9 + -0.0615661 * tfactors.T953 + 0.0454627 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  0.0225613 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 0.815037 * tfactors.T943i + (1.0/3.0) * -0.708732 * tfactors.T923i
                                      + 0.386295 + (5.0/3.0) * -0.0615661 * tfactors.T923 + 0.0454627 * tfactors.T9i;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Cr50_n_to_Cr51_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr50_n_to_Cr51_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr51_n_to_Cr52_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cr51 + n --> Cr52

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  20.1253 + -0.0860167 * tfactors.T9i + 6.95563 * tfactors.T913i + -11.4485 * tfactors.T913
                             + 0.736956 * tfactors.T9 + -0.0695323 * tfactors.T953 + 5.31244 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  0.0860167 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 6.95563 * tfactors.T943i + (1.0/3.0) * -11.4485 * tfactors.T923i
                                      + 0.736956 + (5.0/3.0) * -0.0695323 * tfactors.T923 + 5.31244 * tfactors.T9i;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Cr51_n_to_Cr52_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr51_n_to_Cr52_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V45_n_to_V46_reaclib(T& rate_eval, const tf_t& tfactors) {

        // V45 + n --> V46

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  12.5642 + -1.62632 * tfactors.T913
                             + 0.391829 * tfactors.T9 + -0.0403093 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -1.62632 * tfactors.T923i
                                      + 0.391829 + (5.0/3.0) * -0.0403093 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_V45_n_to_V46_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V45_n_to_V46_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V46_n_to_V47_reaclib(T& rate_eval, const tf_t& tfactors) {

        // V46 + n --> V47

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  15.7405 + -0.692996 * tfactors.T913
                             + 0.316873 * tfactors.T9 + -0.0417235 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -0.692996 * tfactors.T923i
                                      + 0.316873 + (5.0/3.0) * -0.0417235 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_V46_n_to_V47_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V46_n_to_V47_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V47_n_to_V48_reaclib(T& rate_eval, const tf_t& tfactors) {

        // V47 + n --> V48

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  17.5513 + -2.40159 * tfactors.T913
                             + 0.594573 * tfactors.T9 + -0.0682896 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -2.40159 * tfactors.T923i
                                      + 0.594573 + (5.0/3.0) * -0.0682896 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_V47_n_to_V48_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V47_n_to_V48_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V48_n_to_V49_reaclib(T& rate_eval, const tf_t& tfactors) {

        // V48 + n --> V49

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  17.8835 + -0.917026 * tfactors.T913
                             + -0.109162 * tfactors.T9 + 0.0127488 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -0.917026 * tfactors.T923i
                                      + -0.109162 + (5.0/3.0) * 0.0127488 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_V48_n_to_V49_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V48_n_to_V49_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti42_n_to_Ti43_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ti42 + n --> Ti43

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  13.6337 + -1.00856 * tfactors.T913
                             + 0.342738 * tfactors.T9 + -0.0367066 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -1.00856 * tfactors.T923i
                                      + 0.342738 + (5.0/3.0) * -0.0367066 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ti42_n_to_Ti43_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti42_n_to_Ti43_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti43_n_to_Ti44_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ti43 + n --> Ti44

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  11.5528 + -1.09388 * tfactors.T913
                             + 0.248209 * tfactors.T9 + -0.02653 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -1.09388 * tfactors.T923i
                                      + 0.248209 + (5.0/3.0) * -0.02653 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ti43_n_to_Ti44_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti43_n_to_Ti44_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti44_n_to_Ti45_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ti44 + n --> Ti45

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  15.2506 + 1.01203 * tfactors.T913
                             + -0.201174 * tfactors.T9 + 0.00360954 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 1.01203 * tfactors.T923i
                                      + -0.201174 + (5.0/3.0) * 0.00360954 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ti44_n_to_Ti45_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti44_n_to_Ti45_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti45_n_to_Ti46_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ti45 + n --> Ti46

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  17.133 + -1.55633 * tfactors.T913
                             + 0.300783 * tfactors.T9 + -0.0373696 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -1.55633 * tfactors.T923i
                                      + 0.300783 + (5.0/3.0) * -0.0373696 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ti45_n_to_Ti46_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti45_n_to_Ti46_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc41_n_to_Sc42_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Sc41 + n --> Sc42

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  10.6014 + -1.01509 * tfactors.T913
                             + 0.362316 * tfactors.T9 + -0.0412181 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -1.01509 * tfactors.T923i
                                      + 0.362316 + (5.0/3.0) * -0.0412181 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Sc41_n_to_Sc42_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc41_n_to_Sc42_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc42_n_to_Sc43_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Sc42 + n --> Sc43

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  14.9443 + -0.921604 * tfactors.T913
                             + 0.333782 * tfactors.T9 + -0.0404241 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -0.921604 * tfactors.T923i
                                      + 0.333782 + (5.0/3.0) * -0.0404241 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Sc42_n_to_Sc43_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc42_n_to_Sc43_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca38_n_to_Ca39_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ca38 + n --> Ca39

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  12.0675 + -0.418748 * tfactors.T913
                             + 0.243656 * tfactors.T9 + -0.0268305 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -0.418748 * tfactors.T923i
                                      + 0.243656 + (5.0/3.0) * -0.0268305 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ca38_n_to_Ca39_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca38_n_to_Ca39_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca39_n_to_Ca40_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ca39 + n --> Ca40

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  9.33485 + 0.436603 * tfactors.T913
                             + 0.098368 * tfactors.T9 + -0.0226294 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 0.436603 * tfactors.T923i
                                      + 0.098368 + (5.0/3.0) * -0.0226294 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ca39_n_to_Ca40_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca39_n_to_Ca40_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca40_n_to_Ca41_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ca40 + n --> Ca41

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  3.04601 + 0.0855375 * tfactors.T9i + -9.18636 * tfactors.T913i + 21.0095 * tfactors.T913
                             + -1.33326 * tfactors.T9 + 0.0769347 * tfactors.T953 + -8.66245 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  -0.0855375 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -9.18636 * tfactors.T943i + (1.0/3.0) * 21.0095 * tfactors.T923i
                                      + -1.33326 + (5.0/3.0) * 0.0769347 * tfactors.T923 + -8.66245 * tfactors.T9i;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ca40_n_to_Ca41_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca40_n_to_Ca41_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca41_n_to_Ca42_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ca41 + n --> Ca42

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  20.4461 + -0.0409481 * tfactors.T9i + 3.7037 * tfactors.T913i + -9.99246 * tfactors.T913
                             + 1.05894 * tfactors.T9 + -0.0978564 * tfactors.T953 + 3.43604 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  0.0409481 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 3.7037 * tfactors.T943i + (1.0/3.0) * -9.99246 * tfactors.T923i
                                      + 1.05894 + (5.0/3.0) * -0.0978564 * tfactors.T923 + 3.43604 * tfactors.T9i;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ca41_n_to_Ca42_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca41_n_to_Ca42_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar34_n_to_Ar35_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ar34 + n --> Ar35

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  10.6999 + 1.14834 * tfactors.T913
                             + -0.142306 * tfactors.T9 + 0.00918793 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 1.14834 * tfactors.T923i
                                      + -0.142306 + (5.0/3.0) * 0.00918793 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ar34_n_to_Ar35_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar34_n_to_Ar35_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar35_n_to_Ar36_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ar35 + n --> Ar36

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  8.76542 + 0.542023 * tfactors.T913
                             + 0.0223868 * tfactors.T9 + -0.0120062 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 0.542023 * tfactors.T923i
                                      + 0.0223868 + (5.0/3.0) * -0.0120062 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ar35_n_to_Ar36_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar35_n_to_Ar36_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar36_n_to_Ar37_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ar36 + n --> Ar37

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  12.0149 + 0.0317044 * tfactors.T9i + -3.1764 * tfactors.T913i + 5.13191 * tfactors.T913
                             + -0.00639688 * tfactors.T9 + -0.0292833 * tfactors.T953 + -2.74683 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  -0.0317044 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -3.1764 * tfactors.T943i + (1.0/3.0) * 5.13191 * tfactors.T923i
                                      + -0.00639688 + (5.0/3.0) * -0.0292833 * tfactors.T923 + -2.74683 * tfactors.T9i;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ar36_n_to_Ar37_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar36_n_to_Ar37_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar37_n_to_Ar38_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ar37 + n --> Ar38

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  14.7933 + -0.825362 * tfactors.T913
                             + 0.336634 * tfactors.T9 + -0.0509617 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -0.825362 * tfactors.T923i
                                      + 0.336634 + (5.0/3.0) * -0.0509617 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ar37_n_to_Ar38_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar37_n_to_Ar38_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si28_n_to_Si29_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Si28 + n --> Si29

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ka02r
        ln_set_rate =  6.9158 + -0.38 * tfactors.T9i + 7.68863 * tfactors.T913
                             + -1.7991 * tfactors.T9 + -1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  0.38 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 7.68863 * tfactors.T923i
                                      + -1.7991 + -1.5 * tfactors.T9i;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // ka02 
        ln_set_rate =  8.77553;
        amrex::ignore_unused(tfactors);

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 = 0.0;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Si28_n_to_Si29_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si28_n_to_Si29_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si29_n_to_Si30_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Si29 + n --> Si30

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ka02r
        ln_set_rate =  9.60115 + -0.179366 * tfactors.T9i + 5.50678 * tfactors.T913
                             + -2.85656 * tfactors.T9 + -1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  0.179366 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 5.50678 * tfactors.T923i
                                      + -2.85656 + -1.5 * tfactors.T9i;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // ka02n
        ln_set_rate =  11.8023 + 0.650904 * tfactors.T913;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 0.650904 * tfactors.T923i;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Si29_n_to_Si30_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si29_n_to_Si30_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S30_n_to_S31_reaclib(T& rate_eval, const tf_t& tfactors) {

        // S30 + n --> S31

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  10.3285 + 1.62298 * tfactors.T913
                             + -0.278802 * tfactors.T9 + 0.0210647 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 1.62298 * tfactors.T923i
                                      + -0.278802 + (5.0/3.0) * 0.0210647 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_S30_n_to_S31_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S30_n_to_S31_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S31_n_to_S32_reaclib(T& rate_eval, const tf_t& tfactors) {

        // S31 + n --> S32

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  7.56582 + 1.71463 * tfactors.T913
                             + -0.221804 * tfactors.T9 + 0.00880104 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 1.71463 * tfactors.T923i
                                      + -0.221804 + (5.0/3.0) * 0.00880104 * tfactors.T923;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_S31_n_to_S32_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S31_n_to_S32_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S32_n_to_S33_reaclib(T& rate_eval, const tf_t& tfactors) {

        // S32 + n --> S33

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  12.4466 + 0.198828 * tfactors.T9i + -15.0178 * tfactors.T913i + 16.3567 * tfactors.T913
                             + -0.436839 * tfactors.T9 + -0.00574462 * tfactors.T953 + -9.78034 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  -0.198828 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -15.0178 * tfactors.T943i + (1.0/3.0) * 16.3567 * tfactors.T923i
                                      + -0.436839 + (5.0/3.0) * -0.00574462 * tfactors.T923 + -9.78034 * tfactors.T9i;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_S32_n_to_S33_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S32_n_to_S33_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S33_n_to_S34_reaclib(T& rate_eval, const tf_t& tfactors) {

        // S33 + n --> S34

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  -4.17932 + 0.24641 * tfactors.T9i + -23.1351 * tfactors.T913i + 43.3623 * tfactors.T913
                             + -2.74733 * tfactors.T9 + 0.159045 * tfactors.T953 + -19.6332 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  -0.24641 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -23.1351 * tfactors.T943i + (1.0/3.0) * 43.3623 * tfactors.T923i
                                      + -2.74733 + (5.0/3.0) * 0.159045 * tfactors.T923 + -19.6332 * tfactors.T9i;

        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_S33_n_to_S34_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S33_n_to_S34_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_C12_C12_to_p_Na23_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // C12 + C12 --> p + Na23

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // cf88r
        ln_set_rate =  60.9649 + -84.165 * tfactors.T913i + -1.4191 * tfactors.T913
                             + -0.114619 * tfactors.T9 + -0.070307 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -84.165 * tfactors.T943i + (1.0/3.0) * -1.4191 * tfactors.T923i
                                      + -0.114619 + (5.0/3.0) * -0.070307 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_C12_C12_to_p_Na23_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_C12_C12_to_p_Na23_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Na23_p_to_He4_Ne20_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Na23 + p --> He4 + Ne20

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10r
        ln_set_rate =  -6.58736 + -2.31577 * tfactors.T9i + 19.7297 * tfactors.T913
                             + -2.20987 * tfactors.T9 + 0.153374 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  2.31577 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 19.7297 * tfactors.T923i
                                      + -2.20987 + (5.0/3.0) * 0.153374 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10r
        ln_set_rate =  0.0178295 + -1.86103 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  1.86103 * tfactors.T9i * tfactors.T9i
                                      + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10n
        ln_set_rate =  18.9756 + -20.0024 * tfactors.T913i + 11.5988 * tfactors.T913
                             + -1.37398 * tfactors.T9 + -1.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -20.0024 * tfactors.T943i + (1.0/3.0) * 11.5988 * tfactors.T923i
                                      + -1.37398 + (5.0/3.0) * -1.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Na23_p_to_He4_Ne20_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Na23_p_to_He4_Ne20_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Na23_p_to_Mg24_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Na23 + p --> Mg24

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10r
        ln_set_rate =  9.0594 + -3.28029 * tfactors.T9i + -0.360588 * tfactors.T913
                             + 1.4187 * tfactors.T9 + -0.184061 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  3.28029 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.360588 * tfactors.T923i
                                      + 1.4187 + (5.0/3.0) * -0.184061 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10r
        ln_set_rate =  -5.02585 + -1.61219 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  1.61219 * tfactors.T9i * tfactors.T9i
                                      + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10n
        ln_set_rate =  18.9075 + -20.6428 * tfactors.T913i + 1.52954 * tfactors.T913
                             + 2.7487 * tfactors.T9 + -1.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -20.6428 * tfactors.T943i + (1.0/3.0) * 1.52954 * tfactors.T923i
                                      + 2.7487 + (5.0/3.0) * -1.0 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Na23_p_to_Mg24_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Na23_p_to_Mg24_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_C12_C12_to_He4_Ne20_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // C12 + C12 --> He4 + Ne20

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // cf88r
        ln_set_rate =  61.2863 + -84.165 * tfactors.T913i + -1.56627 * tfactors.T913
                             + -0.0736084 * tfactors.T9 + -0.072797 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -84.165 * tfactors.T943i + (1.0/3.0) * -1.56627 * tfactors.T923i
                                      + -0.0736084 + (5.0/3.0) * -0.072797 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_C12_C12_to_He4_Ne20_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_C12_C12_to_He4_Ne20_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ne20_He4_to_Mg24_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ne20 + He4 --> Mg24

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10r
        ln_set_rate =  -8.79827 + -12.7809 * tfactors.T9i + 16.9229 * tfactors.T913
                             + -2.57325 * tfactors.T9 + 0.208997 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  12.7809 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 16.9229 * tfactors.T923i
                                      + -2.57325 + (5.0/3.0) * 0.208997 * tfactors.T923 + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10r
        ln_set_rate =  1.98307 + -9.22026 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  9.22026 * tfactors.T9i * tfactors.T9i
                                      + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10r
        ln_set_rate =  -38.7055 + -2.50605 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  2.50605 * tfactors.T9i * tfactors.T9i
                                      + -1.5 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        // il10n
        ln_set_rate =  24.5058 + -46.2525 * tfactors.T913i + 5.58901 * tfactors.T913
                             + 7.61843 * tfactors.T9 + -3.683 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -46.2525 * tfactors.T943i + (1.0/3.0) * 5.58901 * tfactors.T923i
                                      + 7.61843 + (5.0/3.0) * -3.683 * tfactors.T923 + -0.666667 * tfactors.T9i;

            dln_set_rate_dT9 += dlog_scor_dT * 1.0e9_rt;
        }

        // avoid underflows by zeroing rates in [0.0, 1.e-100]
        ln_set_rate = std::max(ln_set_rate, -230.0_rt);
        set_rate = std::exp(ln_set_rate);
        rate += set_rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            if (ln_set_rate > -230.0_rt) {
                drate_dT += set_rate * dln_set_rate_dT9 * 1.0e-9;
            }
        }

        rate_eval.screened_rates(k_Ne20_He4_to_Mg24_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne20_He4_to_Mg24_reaclib) = drate_dT;
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

        {
            // C12_He4_to_O16_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_C12);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_C12);
            }
#endif
            rate_C12_He4_to_O16_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // O16_He4_to_Ne20_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_O16);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_O16);
            }
#endif
            rate_O16_He4_to_Ne20_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mg23_n_to_Mg24_reaclib
            rate_Mg23_n_to_Mg24_reaclib(rate_eval, tfactors);
        }

        {
            // Mg23_He4_to_Si27_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mg23);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mg23);
            }
#endif
            rate_Mg23_He4_to_Si27_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mg24_He4_to_Si28_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mg24);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mg24);
            }
#endif
            rate_Mg24_He4_to_Si28_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Al27_p_to_Si28_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Al27);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Al27);
            }
#endif
            rate_Al27_p_to_Si28_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Al27_He4_to_P31_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Al27);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Al27);
            }
#endif
            rate_Al27_He4_to_P31_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Si27_n_to_Si28_reaclib
            rate_Si27_n_to_Si28_reaclib(rate_eval, tfactors);
        }

        {
            // Si28_p_to_P29_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Si28);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Si28);
            }
#endif
            rate_Si28_p_to_P29_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Si28_He4_to_S32_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si28);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si28);
            }
#endif
            rate_Si28_He4_to_S32_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Si30_p_to_P31_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Si30);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Si30);
            }
#endif
            rate_Si30_p_to_P31_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Si30_He4_to_S34_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Si30);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Si30);
            }
#endif
            rate_Si30_He4_to_S34_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // P29_p_to_S30_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_P29);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_P29);
            }
#endif
            rate_P29_p_to_S30_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // P29_He4_to_Cl33_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_P29);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_P29);
            }
#endif
            rate_P29_He4_to_Cl33_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // P31_p_to_S32_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_P31);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_P31);
            }
#endif
            rate_P31_p_to_S32_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // P31_He4_to_Cl35_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_P31);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_P31);
            }
#endif
            rate_P31_He4_to_Cl35_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // S30_He4_to_Ar34_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_S30);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_S30);
            }
#endif
            rate_S30_He4_to_Ar34_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // S32_p_to_Cl33_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_S32);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_S32);
            }
#endif
            rate_S32_p_to_Cl33_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // S32_He4_to_Ar36_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_S32);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_S32);
            }
#endif
            rate_S32_He4_to_Ar36_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // S34_p_to_Cl35_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_S34);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_S34);
            }
#endif
            rate_S34_p_to_Cl35_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // S34_He4_to_Ar38_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_S34);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_S34);
            }
#endif
            rate_S34_He4_to_Ar38_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cl33_p_to_Ar34_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cl33);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cl33);
            }
#endif
            rate_Cl33_p_to_Ar34_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cl33_He4_to_K37_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cl33);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cl33);
            }
#endif
            rate_Cl33_He4_to_K37_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cl35_p_to_Ar36_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cl35);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cl35);
            }
#endif
            rate_Cl35_p_to_Ar36_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cl35_He4_to_K39_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cl35);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cl35);
            }
#endif
            rate_Cl35_He4_to_K39_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ar34_He4_to_Ca38_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ar34);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ar34);
            }
#endif
            rate_Ar34_He4_to_Ca38_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ar36_p_to_K37_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ar36);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ar36);
            }
#endif
            rate_Ar36_p_to_K37_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ar36_He4_to_Ca40_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ar36);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ar36);
            }
#endif
            rate_Ar36_He4_to_Ca40_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ar38_p_to_K39_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ar38);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ar38);
            }
#endif
            rate_Ar38_p_to_K39_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ar38_He4_to_Ca42_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ar38);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ar38);
            }
#endif
            rate_Ar38_He4_to_Ca42_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // K37_p_to_Ca38_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_K37);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_K37);
            }
#endif
            rate_K37_p_to_Ca38_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // K37_He4_to_Sc41_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_K37);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_K37);
            }
#endif
            rate_K37_He4_to_Sc41_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // K39_p_to_Ca40_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_K39);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_K39);
            }
#endif
            rate_K39_p_to_Ca40_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // K39_He4_to_Sc43_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_K39);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_K39);
            }
#endif
            rate_K39_He4_to_Sc43_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca38_He4_to_Ti42_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca38);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca38);
            }
#endif
            rate_Ca38_He4_to_Ti42_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca40_p_to_Sc41_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ca40);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca40);
            }
#endif
            rate_Ca40_p_to_Sc41_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca40_He4_to_Ti44_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca40);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca40);
            }
#endif
            rate_Ca40_He4_to_Ti44_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca42_p_to_Sc43_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ca42);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca42);
            }
#endif
            rate_Ca42_p_to_Sc43_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca42_He4_to_Ti46_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca42);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca42);
            }
#endif
            rate_Ca42_He4_to_Ti46_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc41_p_to_Ti42_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Sc41);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc41);
            }
#endif
            rate_Sc41_p_to_Ti42_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc41_He4_to_V45_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc41);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc41);
            }
#endif
            rate_Sc41_He4_to_V45_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc43_p_to_Ti44_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Sc43);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc43);
            }
#endif
            rate_Sc43_p_to_Ti44_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc43_He4_to_V47_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc43);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc43);
            }
#endif
            rate_Sc43_He4_to_V47_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti42_He4_to_Cr46_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti42);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti42);
            }
#endif
            rate_Ti42_He4_to_Cr46_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti44_p_to_V45_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti44);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti44);
            }
#endif
            rate_Ti44_p_to_V45_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti44_He4_to_Cr48_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti44);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti44);
            }
#endif
            rate_Ti44_He4_to_Cr48_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti46_p_to_V47_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ti46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ti46);
            }
#endif
            rate_Ti46_p_to_V47_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti46_He4_to_Cr50_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti46);
            }
#endif
            rate_Ti46_He4_to_Cr50_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V45_p_to_Cr46_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_V45);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V45);
            }
#endif
            rate_V45_p_to_Cr46_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V45_He4_to_Mn49_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_V45);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V45);
            }
#endif
            rate_V45_He4_to_Mn49_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V47_p_to_Cr48_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_V47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V47);
            }
#endif
            rate_V47_p_to_Cr48_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V47_He4_to_Mn51_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_V47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V47);
            }
#endif
            rate_V47_He4_to_Mn51_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V49_p_to_Cr50_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_V49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V49);
            }
#endif
            rate_V49_p_to_Cr50_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V49_He4_to_Mn53_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_V49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V49);
            }
#endif
            rate_V49_He4_to_Mn53_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr46_He4_to_Fe50_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr46);
            }
#endif
            rate_Cr46_He4_to_Fe50_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr48_p_to_Mn49_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cr48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr48);
            }
#endif
            rate_Cr48_p_to_Mn49_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr48_He4_to_Fe52_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr48);
            }
#endif
            rate_Cr48_He4_to_Fe52_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr50_p_to_Mn51_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cr50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr50);
            }
#endif
            rate_Cr50_p_to_Mn51_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr50_He4_to_Fe54_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr50);
            }
#endif
            rate_Cr50_He4_to_Fe54_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr52_p_to_Mn53_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cr52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cr52);
            }
#endif
            rate_Cr52_p_to_Mn53_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr52_He4_to_Fe56_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr52);
            }
#endif
            rate_Cr52_He4_to_Fe56_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn49_p_to_Fe50_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mn49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn49);
            }
#endif
            rate_Mn49_p_to_Fe50_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn49_He4_to_Co53_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn49);
            }
#endif
            rate_Mn49_He4_to_Co53_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn51_p_to_Fe52_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mn51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn51);
            }
#endif
            rate_Mn51_p_to_Fe52_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn51_He4_to_Co55_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn51);
            }
#endif
            rate_Mn51_He4_to_Co55_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn53_p_to_Fe54_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mn53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn53);
            }
#endif
            rate_Mn53_p_to_Fe54_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn53_He4_to_Co57_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn53);
            }
#endif
            rate_Mn53_He4_to_Co57_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Fe52_p_to_Co53_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe52);
            }
#endif
            rate_Fe52_p_to_Co53_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Fe52_He4_to_Ni56_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Fe52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe52);
            }
#endif
            rate_Fe52_He4_to_Ni56_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Fe54_p_to_Co55_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe54);
            }
#endif
            rate_Fe54_p_to_Co55_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Fe54_He4_to_Ni58_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Fe54);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe54);
            }
#endif
            rate_Fe54_He4_to_Ni58_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Fe56_p_to_Co57_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Fe56);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Fe56);
            }
#endif
            rate_Fe56_p_to_Co57_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Co55_n_to_Co56_reaclib
            rate_Co55_n_to_Co56_reaclib(rate_eval, tfactors);
        }

        {
            // Co55_p_to_Ni56_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Co55);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co55);
            }
#endif
            rate_Co55_p_to_Ni56_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Co55_He4_to_Cu59_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Co55);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co55);
            }
#endif
            rate_Co55_He4_to_Cu59_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Co56_n_to_Co57_reaclib
            rate_Co56_n_to_Co57_reaclib(rate_eval, tfactors);
        }

        {
            // Co57_p_to_Ni58_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Co57);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co57);
            }
#endif
            rate_Co57_p_to_Ni58_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ni58_p_to_Cu59_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ni58);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ni58);
            }
#endif
            rate_Ni58_p_to_Cu59_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // O16_C12_to_p_Al27_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_C12_O16);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_O16);
            }
#endif
            rate_O16_C12_to_p_Al27_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // O16_C12_to_He4_Mg24_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_C12_O16);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_O16);
            }
#endif
            rate_O16_C12_to_He4_Mg24_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            rate_O16_O16_to_p_P31_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            rate_O16_O16_to_He4_Si28_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ne20_C12_to_p_P31_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_C12_Ne20);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_Ne20);
            }
#endif
            rate_Ne20_C12_to_p_P31_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ne20_C12_to_He4_Si28_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_C12_Ne20);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_Ne20);
            }
#endif
            rate_Ne20_C12_to_He4_Si28_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mg23_n_to_He4_Ne20_reaclib
            rate_Mg23_n_to_He4_Ne20_reaclib(rate_eval, tfactors);
        }

        {
            // Mg23_n_to_C12_C12_reaclib
            rate_Mg23_n_to_C12_C12_reaclib(rate_eval, tfactors);
        }

        {
            // Al27_p_to_He4_Mg24_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Al27);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Al27);
            }
#endif
            rate_Al27_p_to_He4_Mg24_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Al27_He4_to_p_Si30_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Al27);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Al27);
            }
#endif
            rate_Al27_He4_to_p_Si30_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Si27_n_to_p_Al27_reaclib
            rate_Si27_n_to_p_Al27_reaclib(rate_eval, tfactors);
        }

        {
            // Si27_n_to_He4_Mg24_reaclib
            rate_Si27_n_to_He4_Mg24_reaclib(rate_eval, tfactors);
        }

        {
            // Si27_n_to_C12_O16_reaclib
            rate_Si27_n_to_C12_O16_reaclib(rate_eval, tfactors);
        }

        {
            // P29_He4_to_p_S32_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_P29);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_P29);
            }
#endif
            rate_P29_He4_to_p_S32_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // P31_p_to_He4_Si28_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_P31);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_P31);
            }
#endif
            rate_P31_p_to_He4_Si28_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // P31_He4_to_p_S34_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_P31);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_P31);
            }
#endif
            rate_P31_He4_to_p_S34_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // S30_n_to_He4_Si27_reaclib
            rate_S30_n_to_He4_Si27_reaclib(rate_eval, tfactors);
        }

        {
            // S30_He4_to_p_Cl33_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_S30);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_S30);
            }
#endif
            rate_S30_He4_to_p_Cl33_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cl33_He4_to_p_Ar36_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cl33);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cl33);
            }
#endif
            rate_Cl33_He4_to_p_Ar36_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cl35_p_to_He4_S32_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cl35);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cl35);
            }
#endif
            rate_Cl35_p_to_He4_S32_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cl35_He4_to_p_Ar38_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cl35);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cl35);
            }
#endif
            rate_Cl35_He4_to_p_Ar38_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ar34_He4_to_p_K37_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ar34);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ar34);
            }
#endif
            rate_Ar34_He4_to_p_K37_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // K37_He4_to_p_Ca40_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_K37);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_K37);
            }
#endif
            rate_K37_He4_to_p_Ca40_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // K39_p_to_He4_Ar36_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_K39);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_K39);
            }
#endif
            rate_K39_p_to_He4_Ar36_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca38_He4_to_p_Sc41_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ca38);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ca38);
            }
#endif
            rate_Ca38_He4_to_p_Sc41_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ca42_p_to_He4_K39_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Ca42);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Ca42);
            }
#endif
            rate_Ca42_p_to_He4_K39_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc41_He4_to_p_Ti44_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc41);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc41);
            }
#endif
            rate_Sc41_He4_to_p_Ti44_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc43_p_to_He4_Ca40_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Sc43);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Sc43);
            }
#endif
            rate_Sc43_p_to_He4_Ca40_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Sc43_He4_to_p_Ti46_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Sc43);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Sc43);
            }
#endif
            rate_Sc43_He4_to_p_Ti46_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti42_He4_to_p_V45_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti42);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti42);
            }
#endif
            rate_Ti42_He4_to_p_V45_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ti44_He4_to_p_V47_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ti44);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ti44);
            }
#endif
            rate_Ti44_He4_to_p_V47_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V45_He4_to_p_Cr48_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_V45);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V45);
            }
#endif
            rate_V45_He4_to_p_Cr48_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V47_He4_to_p_Cr50_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_V47);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V47);
            }
#endif
            rate_V47_He4_to_p_Cr50_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V49_p_to_He4_Ti46_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_V49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_V49);
            }
#endif
            rate_V49_p_to_He4_Ti46_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // V49_He4_to_p_Cr52_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_V49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_V49);
            }
#endif
            rate_V49_He4_to_p_Cr52_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr46_He4_to_p_Mn49_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr46);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr46);
            }
#endif
            rate_Cr46_He4_to_p_Mn49_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Cr48_He4_to_p_Mn51_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Cr48);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Cr48);
            }
#endif
            rate_Cr48_He4_to_p_Mn51_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn49_He4_to_p_Fe52_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn49);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn49);
            }
#endif
            rate_Mn49_He4_to_p_Fe52_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn51_He4_to_p_Fe54_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn51);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn51);
            }
#endif
            rate_Mn51_He4_to_p_Fe54_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn53_p_to_He4_Cr50_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mn53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mn53);
            }
#endif
            rate_Mn53_p_to_He4_Cr50_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mn53_He4_to_p_Fe56_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Mn53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Mn53);
            }
#endif
            rate_Mn53_He4_to_p_Fe56_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Fe50_He4_to_p_Co53_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Fe50);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe50);
            }
#endif
            rate_Fe50_He4_to_p_Co53_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Fe52_He4_to_p_Co55_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Fe52);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Fe52);
            }
#endif
            rate_Fe52_He4_to_p_Co55_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Co53_He4_to_p_Ni56_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Co53);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co53);
            }
#endif
            rate_Co53_He4_to_p_Ni56_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Co55_He4_to_p_Ni58_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Co55);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Co55);
            }
#endif
            rate_Co55_He4_to_p_Ni58_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Co56_n_to_p_Fe56_reaclib
            rate_Co56_n_to_p_Fe56_reaclib(rate_eval, tfactors);
        }

        {
            // Co56_n_to_He4_Mn53_reaclib
            rate_Co56_n_to_He4_Mn53_reaclib(rate_eval, tfactors);
        }

        {
            // Co57_p_to_He4_Fe54_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Co57);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Co57);
            }
#endif
            rate_Co57_p_to_He4_Fe54_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ni56_n_to_p_Co56_reaclib
            rate_Ni56_n_to_p_Co56_reaclib(rate_eval, tfactors);
        }

        {
            // Cu59_n_to_He4_Co56_reaclib
            rate_Cu59_n_to_He4_Co56_reaclib(rate_eval, tfactors);
        }

        {
            // Cu59_p_to_He4_Ni56_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Cu59);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Cu59);
            }
#endif
            rate_Cu59_p_to_He4_Ni56_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            rate_He4_He4_He4_to_C12_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ni56_n_to_Ni57_reaclib
            rate_Ni56_n_to_Ni57_reaclib(rate_eval, tfactors);
        }

        {
            // Ni57_n_to_Ni58_reaclib
            rate_Ni57_n_to_Ni58_reaclib(rate_eval, tfactors);
        }

        {
            // Co53_n_to_Co54_reaclib
            rate_Co53_n_to_Co54_reaclib(rate_eval, tfactors);
        }

        {
            // Co54_n_to_Co55_reaclib
            rate_Co54_n_to_Co55_reaclib(rate_eval, tfactors);
        }

        {
            // Fe50_n_to_Fe51_reaclib
            rate_Fe50_n_to_Fe51_reaclib(rate_eval, tfactors);
        }

        {
            // Fe51_n_to_Fe52_reaclib
            rate_Fe51_n_to_Fe52_reaclib(rate_eval, tfactors);
        }

        {
            // Fe52_n_to_Fe53_reaclib
            rate_Fe52_n_to_Fe53_reaclib(rate_eval, tfactors);
        }

        {
            // Fe53_n_to_Fe54_reaclib
            rate_Fe53_n_to_Fe54_reaclib(rate_eval, tfactors);
        }

        {
            // Fe54_n_to_Fe55_reaclib
            rate_Fe54_n_to_Fe55_reaclib(rate_eval, tfactors);
        }

        {
            // Fe55_n_to_Fe56_reaclib
            rate_Fe55_n_to_Fe56_reaclib(rate_eval, tfactors);
        }

        {
            // Mn49_n_to_Mn50_reaclib
            rate_Mn49_n_to_Mn50_reaclib(rate_eval, tfactors);
        }

        {
            // Mn50_n_to_Mn51_reaclib
            rate_Mn50_n_to_Mn51_reaclib(rate_eval, tfactors);
        }

        {
            // Mn51_n_to_Mn52_reaclib
            rate_Mn51_n_to_Mn52_reaclib(rate_eval, tfactors);
        }

        {
            // Mn52_n_to_Mn53_reaclib
            rate_Mn52_n_to_Mn53_reaclib(rate_eval, tfactors);
        }

        {
            // Cr46_n_to_Cr47_reaclib
            rate_Cr46_n_to_Cr47_reaclib(rate_eval, tfactors);
        }

        {
            // Cr47_n_to_Cr48_reaclib
            rate_Cr47_n_to_Cr48_reaclib(rate_eval, tfactors);
        }

        {
            // Cr48_n_to_Cr49_reaclib
            rate_Cr48_n_to_Cr49_reaclib(rate_eval, tfactors);
        }

        {
            // Cr49_n_to_Cr50_reaclib
            rate_Cr49_n_to_Cr50_reaclib(rate_eval, tfactors);
        }

        {
            // Cr50_n_to_Cr51_reaclib
            rate_Cr50_n_to_Cr51_reaclib(rate_eval, tfactors);
        }

        {
            // Cr51_n_to_Cr52_reaclib
            rate_Cr51_n_to_Cr52_reaclib(rate_eval, tfactors);
        }

        {
            // V45_n_to_V46_reaclib
            rate_V45_n_to_V46_reaclib(rate_eval, tfactors);
        }

        {
            // V46_n_to_V47_reaclib
            rate_V46_n_to_V47_reaclib(rate_eval, tfactors);
        }

        {
            // V47_n_to_V48_reaclib
            rate_V47_n_to_V48_reaclib(rate_eval, tfactors);
        }

        {
            // V48_n_to_V49_reaclib
            rate_V48_n_to_V49_reaclib(rate_eval, tfactors);
        }

        {
            // Ti42_n_to_Ti43_reaclib
            rate_Ti42_n_to_Ti43_reaclib(rate_eval, tfactors);
        }

        {
            // Ti43_n_to_Ti44_reaclib
            rate_Ti43_n_to_Ti44_reaclib(rate_eval, tfactors);
        }

        {
            // Ti44_n_to_Ti45_reaclib
            rate_Ti44_n_to_Ti45_reaclib(rate_eval, tfactors);
        }

        {
            // Ti45_n_to_Ti46_reaclib
            rate_Ti45_n_to_Ti46_reaclib(rate_eval, tfactors);
        }

        {
            // Sc41_n_to_Sc42_reaclib
            rate_Sc41_n_to_Sc42_reaclib(rate_eval, tfactors);
        }

        {
            // Sc42_n_to_Sc43_reaclib
            rate_Sc42_n_to_Sc43_reaclib(rate_eval, tfactors);
        }

        {
            // Ca38_n_to_Ca39_reaclib
            rate_Ca38_n_to_Ca39_reaclib(rate_eval, tfactors);
        }

        {
            // Ca39_n_to_Ca40_reaclib
            rate_Ca39_n_to_Ca40_reaclib(rate_eval, tfactors);
        }

        {
            // Ca40_n_to_Ca41_reaclib
            rate_Ca40_n_to_Ca41_reaclib(rate_eval, tfactors);
        }

        {
            // Ca41_n_to_Ca42_reaclib
            rate_Ca41_n_to_Ca42_reaclib(rate_eval, tfactors);
        }

        {
            // Ar34_n_to_Ar35_reaclib
            rate_Ar34_n_to_Ar35_reaclib(rate_eval, tfactors);
        }

        {
            // Ar35_n_to_Ar36_reaclib
            rate_Ar35_n_to_Ar36_reaclib(rate_eval, tfactors);
        }

        {
            // Ar36_n_to_Ar37_reaclib
            rate_Ar36_n_to_Ar37_reaclib(rate_eval, tfactors);
        }

        {
            // Ar37_n_to_Ar38_reaclib
            rate_Ar37_n_to_Ar38_reaclib(rate_eval, tfactors);
        }

        {
            // Si28_n_to_Si29_reaclib
            rate_Si28_n_to_Si29_reaclib(rate_eval, tfactors);
        }

        {
            // Si29_n_to_Si30_reaclib
            rate_Si29_n_to_Si30_reaclib(rate_eval, tfactors);
        }

        {
            // S30_n_to_S31_reaclib
            rate_S30_n_to_S31_reaclib(rate_eval, tfactors);
        }

        {
            // S31_n_to_S32_reaclib
            rate_S31_n_to_S32_reaclib(rate_eval, tfactors);
        }

        {
            // S32_n_to_S33_reaclib
            rate_S32_n_to_S33_reaclib(rate_eval, tfactors);
        }

        {
            // S33_n_to_S34_reaclib
            rate_S33_n_to_S34_reaclib(rate_eval, tfactors);
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
            rate_C12_C12_to_p_Na23_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Na23_p_to_He4_Ne20_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Na23);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Na23);
            }
#endif
            rate_Na23_p_to_He4_Ne20_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Na23_p_to_Mg24_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Na23);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Na23);
            }
#endif
            rate_Na23_p_to_Mg24_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
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
            rate_C12_C12_to_He4_Ne20_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Ne20_He4_to_Mg24_reaclib
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Ne20);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Ne20);
            }
#endif
            rate_Ne20_He4_to_Mg24_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }


    }
}
#endif
