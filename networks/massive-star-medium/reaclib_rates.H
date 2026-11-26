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
    void rate_O16_to_He4_C12_reaclib(T& rate_eval, const tf_t& tfactors) {

        // O16 --> He4 + C12

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // nac2 
        ln_set_rate =  94.3131 + -84.503 * tfactors.T9i + 58.9128 * tfactors.T913i + -148.273 * tfactors.T913
                             + 9.08324 * tfactors.T9 + -0.541041 * tfactors.T953 + 71.8554 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  84.503 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 58.9128 * tfactors.T943i + (1.0/3.0) * -148.273 * tfactors.T923i
                                      + 9.08324 + (5.0/3.0) * -0.541041 * tfactors.T923 + 71.8554 * tfactors.T9i;

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
        ln_set_rate =  279.295 + -84.9515 * tfactors.T9i + 103.411 * tfactors.T913i + -420.567 * tfactors.T913
                             + 64.0874 * tfactors.T9 + -12.4624 * tfactors.T953 + 138.803 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  84.9515 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 103.411 * tfactors.T943i + (1.0/3.0) * -420.567 * tfactors.T923i
                                      + 64.0874 + (5.0/3.0) * -12.4624 * tfactors.T923 + 138.803 * tfactors.T9i;

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

        rate_eval.screened_rates(k_O16_to_He4_C12_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_O16_to_He4_C12_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ne20_to_He4_O16_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ne20 --> He4 + O16

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // co10r
        ln_set_rate =  28.6431 + -65.246 * tfactors.T9i;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  65.246 * tfactors.T9i * tfactors.T9i;

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
        ln_set_rate =  48.6604 + -54.8875 * tfactors.T9i + -39.7262 * tfactors.T913i + -0.210799 * tfactors.T913
                             + 0.442879 * tfactors.T9 + -0.0797753 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  54.8875 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -39.7262 * tfactors.T943i + (1.0/3.0) * -0.210799 * tfactors.T923i
                                      + 0.442879 + (5.0/3.0) * -0.0797753 * tfactors.T923 + 0.833333 * tfactors.T9i;

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
        ln_set_rate =  34.2658 + -67.6518 * tfactors.T9i + -3.65925 * tfactors.T913
                             + 0.714224 * tfactors.T9 + -0.00107508 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  67.6518 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -3.65925 * tfactors.T923i
                                      + 0.714224 + (5.0/3.0) * -0.00107508 * tfactors.T923;

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

        rate_eval.screened_rates(k_Ne20_to_He4_O16_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne20_to_He4_O16_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mg24_to_n_Mg23_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Mg24 --> n + Mg23

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  32.0344 + -191.835 * tfactors.T9i + 2.66964 * tfactors.T913
                             + -0.448904 * tfactors.T9 + 0.0326505 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  191.835 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 2.66964 * tfactors.T923i
                                      + -0.448904 + (5.0/3.0) * 0.0326505 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Mg24_to_n_Mg23_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg24_to_n_Mg23_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si27_to_He4_Mg23_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Si27 --> He4 + Mg23

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  68.1469 + -108.333 * tfactors.T9i + -53.203 * tfactors.T913i + -4.6318 * tfactors.T913
                             + -0.130951 * tfactors.T9 + 0.014691 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  108.333 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -53.203 * tfactors.T943i + (1.0/3.0) * -4.6318 * tfactors.T923i
                                      + -0.130951 + (5.0/3.0) * 0.014691 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Si27_to_He4_Mg23_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si27_to_He4_Mg23_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si28_to_n_Si27_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Si28 --> n + Si27

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  34.5835 + -199.363 * tfactors.T9i + 0.559804 * tfactors.T913
                             + -0.0233808 * tfactors.T9 + -0.00121313 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  199.363 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 0.559804 * tfactors.T923i
                                      + -0.0233808 + (5.0/3.0) * -0.00121313 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Si28_to_n_Si27_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si28_to_n_Si27_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si28_to_p_Al27_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Si28 --> p + Al27

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10n
        ln_set_rate =  46.5494 + -134.445 * tfactors.T9i + -23.2205 * tfactors.T913i
                             + -2.0 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  134.445 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -23.2205 * tfactors.T943i
                                      + (5.0/3.0) * -2.0 * tfactors.T923 + 0.833333 * tfactors.T9i;

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
        ln_set_rate =  111.466 + -134.832 * tfactors.T9i + -26.8327 * tfactors.T913i + -116.137 * tfactors.T913
                             + 0.00950567 * tfactors.T9 + 0.00999755 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  134.832 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -26.8327 * tfactors.T943i + (1.0/3.0) * -116.137 * tfactors.T923i
                                      + 0.00950567 + (5.0/3.0) * 0.00999755 * tfactors.T923;

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
        ln_set_rate =  11.7765 + -136.349 * tfactors.T9i + 23.8634 * tfactors.T913
                             + -3.70135 * tfactors.T9 + 0.28964 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  136.349 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 23.8634 * tfactors.T923i
                                      + -3.70135 + (5.0/3.0) * 0.28964 * tfactors.T923;

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

        rate_eval.screened_rates(k_Si28_to_p_Al27_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si28_to_p_Al27_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si28_to_He4_Mg24_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Si28 --> He4 + Mg24

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // st08r
        ln_set_rate =  32.9006 + -131.488 * tfactors.T9i;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  131.488 * tfactors.T9i * tfactors.T9i;

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
        ln_set_rate =  -25.6886 + -128.693 * tfactors.T9i + 21.3721 * tfactors.T913i + 37.7649 * tfactors.T913
                             + -4.10635 * tfactors.T9 + 0.249618 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  128.693 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 21.3721 * tfactors.T943i + (1.0/3.0) * 37.7649 * tfactors.T923i
                                      + -4.10635 + (5.0/3.0) * 0.249618 * tfactors.T923;

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

        rate_eval.screened_rates(k_Si28_to_He4_Mg24_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si28_to_He4_Mg24_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si29_to_n_Si28_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Si29 --> n + Si28

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ka02 
        ln_set_rate =  31.7355 + -98.3365 * tfactors.T9i
                             + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  98.3365 * tfactors.T9i * tfactors.T9i
                                      + 1.5 * tfactors.T9i;

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

        // ka02r
        ln_set_rate =  29.8758 + -98.7165 * tfactors.T9i + 7.68863 * tfactors.T913
                             + -1.7991 * tfactors.T9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  98.7165 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 7.68863 * tfactors.T923i
                                      + -1.7991;

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

        rate_eval.screened_rates(k_Si29_to_n_Si28_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si29_to_n_Si28_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si30_to_n_Si29_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Si30 --> n + Si29

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ka02n
        ln_set_rate =  36.1504 + -123.112 * tfactors.T9i + 0.650904 * tfactors.T913
                             + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  123.112 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 0.650904 * tfactors.T923i
                                      + 1.5 * tfactors.T9i;

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

        // ka02r
        ln_set_rate =  33.9492 + -123.292 * tfactors.T9i + 5.50678 * tfactors.T913
                             + -2.85656 * tfactors.T9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  123.292 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 5.50678 * tfactors.T923i
                                      + -2.85656;

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

        rate_eval.screened_rates(k_Si30_to_n_Si29_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si30_to_n_Si29_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si31_to_n_Si30_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Si31 --> n + Si30

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  13.0759 + -76.3737 * tfactors.T9i + -12.3502 * tfactors.T913i + 34.6486 * tfactors.T913
                             + -1.89409 * tfactors.T9 + 0.0781979 * tfactors.T953 + -12.9601 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  76.3737 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -12.3502 * tfactors.T943i + (1.0/3.0) * 34.6486 * tfactors.T923i
                                      + -1.89409 + (5.0/3.0) * 0.0781979 * tfactors.T923 + -12.9601 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Si31_to_n_Si30_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si31_to_n_Si30_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P29_to_p_Si28_reaclib(T& rate_eval, const tf_t& tfactors) {

        // P29 --> p + Si28

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10n
        ln_set_rate =  39.1379 + -31.8984 * tfactors.T9i + -23.8173 * tfactors.T913i + 7.08203 * tfactors.T913
                             + -1.44753 * tfactors.T9 + 0.0804296 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  31.8984 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -23.8173 * tfactors.T943i + (1.0/3.0) * 7.08203 * tfactors.T923i
                                      + -1.44753 + (5.0/3.0) * 0.0804296 * tfactors.T923 + 0.833333 * tfactors.T9i;

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
        ln_set_rate =  28.6997 + -36.0408 * tfactors.T9i;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  36.0408 * tfactors.T9i * tfactors.T9i;

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

        rate_eval.screened_rates(k_P29_to_p_Si28_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P29_to_p_Si28_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P31_to_p_Si30_reaclib(T& rate_eval, const tf_t& tfactors) {

        // P31 --> p + Si30

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // de20r
        ln_set_rate =  18.7213 + -85.9282 * tfactors.T9i
                             + 6.49034 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  85.9282 * tfactors.T9i * tfactors.T9i
                                      + 6.49034 * tfactors.T9i;

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
        ln_set_rate =  4.04359 + -85.6217 * tfactors.T9i
                             + 2.80331 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  85.6217 * tfactors.T9i * tfactors.T9i
                                      + 2.80331 * tfactors.T9i;

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
        ln_set_rate =  -1115.38 + -180.553 * tfactors.T9i
                             + -895.258 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  180.553 * tfactors.T9i * tfactors.T9i
                                      + -895.258 * tfactors.T9i;

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
        ln_set_rate =  32.9288 + -90.2661 * tfactors.T9i
                             + -0.070816 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  90.2661 * tfactors.T9i * tfactors.T9i
                                      + -0.070816 * tfactors.T9i;

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
        ln_set_rate =  -11.8961 + -85.2694 * tfactors.T9i
                             + -0.128387 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  85.2694 * tfactors.T9i * tfactors.T9i
                                      + -0.128387 * tfactors.T9i;

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
        ln_set_rate =  35.3329 + -91.3175 * tfactors.T9i
                             + 0.3809 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  91.3175 * tfactors.T9i * tfactors.T9i
                                      + 0.3809 * tfactors.T9i;

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
        ln_set_rate =  31.761 + -89.8588 * tfactors.T9i
                             + 2.7883 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  89.8588 * tfactors.T9i * tfactors.T9i
                                      + 2.7883 * tfactors.T9i;

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
        ln_set_rate =  -311.303 + -85.8097 * tfactors.T9i
                             + -77.047 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  85.8097 * tfactors.T9i * tfactors.T9i
                                      + -77.047 * tfactors.T9i;

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

        rate_eval.screened_rates(k_P31_to_p_Si30_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P31_to_p_Si30_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P31_to_He4_Al27_reaclib(T& rate_eval, const tf_t& tfactors) {

        // P31 --> He4 + Al27

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  73.2168 + -112.206 * tfactors.T9i + -56.5351 * tfactors.T913i + -0.896208 * tfactors.T913
                             + -1.72024 * tfactors.T9 + 0.185409 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  112.206 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -56.5351 * tfactors.T943i + (1.0/3.0) * -0.896208 * tfactors.T923i
                                      + -1.72024 + (5.0/3.0) * 0.185409 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_P31_to_He4_Al27_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P31_to_He4_Al27_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S32_to_p_P31_reaclib(T& rate_eval, const tf_t& tfactors) {

        // S32 --> p + P31

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10r
        ln_set_rate =  25.1729 + -106.637 * tfactors.T9i + 8.09341 * tfactors.T913
                             + -0.615971 * tfactors.T9 + 0.031159 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  106.637 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 8.09341 * tfactors.T923i
                                      + -0.615971 + (5.0/3.0) * 0.031159 * tfactors.T923;

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
        ln_set_rate =  21.6829 + -105.119 * tfactors.T9i;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  105.119 * tfactors.T9i * tfactors.T9i;

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
        ln_set_rate =  43.6109 + -102.86 * tfactors.T9i + -25.3278 * tfactors.T913i + 6.4931 * tfactors.T913
                             + -9.27513 * tfactors.T9 + -0.610439 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  102.86 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -25.3278 * tfactors.T943i + (1.0/3.0) * 6.4931 * tfactors.T923i
                                      + -9.27513 + (5.0/3.0) * -0.610439 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_S32_to_p_P31_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S32_to_p_P31_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S32_to_He4_Si28_reaclib(T& rate_eval, const tf_t& tfactors) {

        // S32 --> He4 + Si28

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  72.813 + -80.626 * tfactors.T9i + -59.4896 * tfactors.T913i + 4.47205 * tfactors.T913
                             + -4.78989 * tfactors.T9 + 0.557201 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  80.626 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -59.4896 * tfactors.T943i + (1.0/3.0) * 4.47205 * tfactors.T923i
                                      + -4.78989 + (5.0/3.0) * 0.557201 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_S32_to_He4_Si28_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S32_to_He4_Si28_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S33_to_n_S32_reaclib(T& rate_eval, const tf_t& tfactors) {

        // S33 --> n + S32

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  34.7199 + -100.079 * tfactors.T9i + -15.0178 * tfactors.T913i + 16.3567 * tfactors.T913
                             + -0.436839 * tfactors.T9 + -0.00574462 * tfactors.T953 + -8.28034 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  100.079 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -15.0178 * tfactors.T943i + (1.0/3.0) * 16.3567 * tfactors.T923i
                                      + -0.436839 + (5.0/3.0) * -0.00574462 * tfactors.T923 + -8.28034 * tfactors.T9i;

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

        rate_eval.screened_rates(k_S33_to_n_S32_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S33_to_n_S32_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S33_to_He4_Si29_reaclib(T& rate_eval, const tf_t& tfactors) {

        // S33 --> He4 + Si29

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  73.7708 + -82.576 * tfactors.T9i + -59.5755 * tfactors.T913i + 1.06274 * tfactors.T913
                             + -3.07529 * tfactors.T9 + 0.372011 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  82.576 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -59.5755 * tfactors.T943i + (1.0/3.0) * 1.06274 * tfactors.T923i
                                      + -3.07529 + (5.0/3.0) * 0.372011 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_S33_to_He4_Si29_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S33_to_He4_Si29_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S34_to_n_S33_reaclib(T& rate_eval, const tf_t& tfactors) {

        // S34 --> n + S33

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  20.868 + -132.246 * tfactors.T9i + -23.1351 * tfactors.T913i + 43.3623 * tfactors.T913
                             + -2.74733 * tfactors.T9 + 0.159045 * tfactors.T953 + -18.1332 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  132.246 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -23.1351 * tfactors.T943i + (1.0/3.0) * 43.3623 * tfactors.T923i
                                      + -2.74733 + (5.0/3.0) * 0.159045 * tfactors.T923 + -18.1332 * tfactors.T9i;

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

        rate_eval.screened_rates(k_S34_to_n_S33_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S34_to_n_S33_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S34_to_p_P33_reaclib(T& rate_eval, const tf_t& tfactors) {

        // S34 --> p + P33

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // rath 
        ln_set_rate =  198.208 + -130.135 * tfactors.T9i + 161.691 * tfactors.T913i + -350.015 * tfactors.T913
                             + 24.3597 * tfactors.T9 + -1.56396 * tfactors.T953 + 151.463 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  130.135 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 161.691 * tfactors.T943i + (1.0/3.0) * -350.015 * tfactors.T923i
                                      + 24.3597 + (5.0/3.0) * -1.56396 * tfactors.T923 + 151.463 * tfactors.T9i;

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

        rate_eval.screened_rates(k_S34_to_p_P33_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S34_to_p_P33_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S34_to_He4_Si30_reaclib(T& rate_eval, const tf_t& tfactors) {

        // S34 --> He4 + Si30

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  72.4051 + -91.9515 * tfactors.T9i + -59.6559 * tfactors.T913i + 3.70141 * tfactors.T913
                             + -3.12537 * tfactors.T9 + 0.307626 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  91.9515 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -59.6559 * tfactors.T943i + (1.0/3.0) * 3.70141 * tfactors.T923i
                                      + -3.12537 + (5.0/3.0) * 0.307626 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_S34_to_He4_Si30_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S34_to_He4_Si30_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S35_to_n_S34_reaclib(T& rate_eval, const tf_t& tfactors) {

        // S35 --> n + S34

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  63.6891 + -81.2902 * tfactors.T9i + 25.2182 * tfactors.T913i + -59.2704 * tfactors.T913
                             + 4.22256 * tfactors.T9 + -0.274247 * tfactors.T953 + 26.032 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  81.2902 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 25.2182 * tfactors.T943i + (1.0/3.0) * -59.2704 * tfactors.T923i
                                      + 4.22256 + (5.0/3.0) * -0.274247 * tfactors.T923 + 26.032 * tfactors.T9i;

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

        rate_eval.screened_rates(k_S35_to_n_S34_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S35_to_n_S34_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S35_to_He4_Si31_reaclib(T& rate_eval, const tf_t& tfactors) {

        // S35 --> He4 + Si31

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // rath 
        ln_set_rate =  -48.3481 + -98.5134 * tfactors.T9i + -98.5176 * tfactors.T913i + 171.207 * tfactors.T913
                             + -14.3839 * tfactors.T9 + 0.960708 * tfactors.T953 + -60.4936 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  98.5134 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -98.5176 * tfactors.T943i + (1.0/3.0) * 171.207 * tfactors.T923i
                                      + -14.3839 + (5.0/3.0) * 0.960708 * tfactors.T923 + -60.4936 * tfactors.T9i;

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

        rate_eval.screened_rates(k_S35_to_He4_Si31_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S35_to_He4_Si31_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl33_to_p_S32_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cl33 --> p + S32

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10r
        ln_set_rate =  91.6078 + -29.4248 * tfactors.T9i + -33.7204 * tfactors.T913i + -32.7355 * tfactors.T913
                             + 3.92526 * tfactors.T9 + -0.250479 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  29.4248 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -33.7204 * tfactors.T943i + (1.0/3.0) * -32.7355 * tfactors.T923i
                                      + 3.92526 + (5.0/3.0) * -0.250479 * tfactors.T923;

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
        ln_set_rate =  74.7432 + -26.4211 * tfactors.T9i + -29.7741 * tfactors.T913i + -87.4473 * tfactors.T913
                             + 182.189 * tfactors.T9 + -128.625 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  26.4211 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -29.7741 * tfactors.T943i + (1.0/3.0) * -87.4473 * tfactors.T923i
                                      + 182.189 + (5.0/3.0) * -128.625 * tfactors.T923 + 0.833333 * tfactors.T9i;

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
        ln_set_rate =  -4.96497 + -27.2952 * tfactors.T9i;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  27.2952 * tfactors.T9i * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cl33_to_p_S32_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl33_to_p_S32_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl33_to_He4_P29_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cl33 --> He4 + P29

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  66.203 + -75.1475 * tfactors.T9i + -62.3802 * tfactors.T913i + 0.593062 * tfactors.T913
                             + -1.14226 * tfactors.T9 + 0.0934776 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  75.1475 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -62.3802 * tfactors.T943i + (1.0/3.0) * 0.593062 * tfactors.T923i
                                      + -1.14226 + (5.0/3.0) * 0.0934776 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cl33_to_He4_P29_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl33_to_He4_P29_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl35_to_p_S34_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cl35 --> p + S34

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // se19r
        ln_set_rate =  22.7104 + -77.1013 * tfactors.T9i;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  77.1013 * tfactors.T9i * tfactors.T9i;

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
        ln_set_rate =  -37.7488 + -101.622 * tfactors.T9i + 454.53 * tfactors.T913i + -375.557 * tfactors.T913
                             + 12.6533 * tfactors.T9 + -0.408677 * tfactors.T953 + 246.618 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  101.622 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 454.53 * tfactors.T943i + (1.0/3.0) * -375.557 * tfactors.T923i
                                      + 12.6533 + (5.0/3.0) * -0.408677 * tfactors.T923 + 246.618 * tfactors.T9i;

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
        ln_set_rate =  14.5754 + -75.9594 * tfactors.T9i;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  75.9594 * tfactors.T9i * tfactors.T9i;

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
        ln_set_rate =  30.5303 + -79.6891 * tfactors.T9i;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  79.6891 * tfactors.T9i * tfactors.T9i;

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
        ln_set_rate =  7.39781 + -75.3374 * tfactors.T9i;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  75.3374 * tfactors.T9i * tfactors.T9i;

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
        ln_set_rate =  29.2773 + -79.0392 * tfactors.T9i;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  79.0392 * tfactors.T9i * tfactors.T9i;

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
        ln_set_rate =  4.33927 + -75.1053 * tfactors.T9i;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  75.1053 * tfactors.T9i * tfactors.T9i;

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
        ln_set_rate =  28.0819 + -78.5751 * tfactors.T9i;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  78.5751 * tfactors.T9i * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cl35_to_p_S34_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl35_to_p_S34_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl35_to_He4_P31_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cl35 --> He4 + P31

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  74.6679 + -81.2033 * tfactors.T9i + -62.5433 * tfactors.T913i + -2.95026 * tfactors.T913
                             + -0.89652 * tfactors.T9 + 0.0774126 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  81.2033 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -62.5433 * tfactors.T943i + (1.0/3.0) * -2.95026 * tfactors.T923i
                                      + -0.89652 + (5.0/3.0) * 0.0774126 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cl35_to_He4_P31_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl35_to_He4_P31_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl37_to_He4_P33_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cl37 --> He4 + P33

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // rath 
        ln_set_rate =  75.154 + -94.4688 * tfactors.T9i + 21.2001 * tfactors.T913i + -84.0086 * tfactors.T913
                             + 1.39426 * tfactors.T9 + 0.0896178 * tfactors.T953 + 49.3178 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  94.4688 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 21.2001 * tfactors.T943i + (1.0/3.0) * -84.0086 * tfactors.T923i
                                      + 1.39426 + (5.0/3.0) * 0.0896178 * tfactors.T923 + 49.3178 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cl37_to_He4_P33_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl37_to_He4_P33_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar36_to_p_Cl35_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ar36 --> p + Cl35

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10r
        ln_set_rate =  17.2028 + -102.37 * tfactors.T9i + 18.0179 * tfactors.T913
                             + -2.86304 * tfactors.T9 + 0.250854 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  102.37 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 18.0179 * tfactors.T923i
                                      + -2.86304 + (5.0/3.0) * 0.250854 * tfactors.T923;

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
        ln_set_rate =  16.0169 + -100.729 * tfactors.T9i;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  100.729 * tfactors.T9i * tfactors.T9i;

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
        ln_set_rate =  -17.4751 + -99.2838 * tfactors.T9i;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  99.2838 * tfactors.T9i * tfactors.T9i;

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
        ln_set_rate =  60.7366 + -98.7191 * tfactors.T9i + -27.8971 * tfactors.T913i + -16.2304 * tfactors.T913
                             + 35.255 * tfactors.T9 + -25.8411 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  98.7191 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -27.8971 * tfactors.T943i + (1.0/3.0) * -16.2304 * tfactors.T923i
                                      + 35.255 + (5.0/3.0) * -25.8411 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ar36_to_p_Cl35_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar36_to_p_Cl35_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar36_to_He4_S32_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ar36 --> He4 + S32

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  73.8164 + -77.0627 * tfactors.T9i + -65.3709 * tfactors.T913i + 5.68294 * tfactors.T913
                             + -5.00388 * tfactors.T9 + 0.571407 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  77.0627 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -65.3709 * tfactors.T943i + (1.0/3.0) * 5.68294 * tfactors.T923i
                                      + -5.00388 + (5.0/3.0) * 0.571407 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ar36_to_He4_S32_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar36_to_He4_S32_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar37_to_n_Ar36_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ar37 --> n + Ar36

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  34.2933 + -101.941 * tfactors.T9i + -3.1764 * tfactors.T913i + 5.13191 * tfactors.T913
                             + -0.00639688 * tfactors.T9 + -0.0292833 * tfactors.T953 + -1.24683 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  101.941 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -3.1764 * tfactors.T943i + (1.0/3.0) * 5.13191 * tfactors.T923i
                                      + -0.00639688 + (5.0/3.0) * -0.0292833 * tfactors.T923 + -1.24683 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ar37_to_n_Ar36_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar37_to_n_Ar36_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar37_to_He4_S33_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ar37 --> He4 + S33

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  74.852 + -78.7549 * tfactors.T9i + -65.4446 * tfactors.T913i + 3.59607 * tfactors.T913
                             + -3.40501 * tfactors.T9 + 0.363961 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  78.7549 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -65.4446 * tfactors.T943i + (1.0/3.0) * 3.59607 * tfactors.T923i
                                      + -3.40501 + (5.0/3.0) * 0.363961 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ar37_to_He4_S33_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar37_to_He4_S33_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar38_to_n_Ar37_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ar38 --> n + Ar37

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  39.8454 + -137.376 * tfactors.T9i + -0.825362 * tfactors.T913
                             + 0.336634 * tfactors.T9 + -0.0509617 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  137.376 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.825362 * tfactors.T923i
                                      + 0.336634 + (5.0/3.0) * -0.0509617 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ar38_to_n_Ar37_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar38_to_n_Ar37_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar38_to_p_Cl37_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ar38 --> p + Cl37

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  59.0616 + -118.853 * tfactors.T9i + -27.9113 * tfactors.T913i + 0.282028 * tfactors.T913
                             + -1.80122 * tfactors.T9 + 0.21751 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  118.853 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -27.9113 * tfactors.T943i + (1.0/3.0) * 0.282028 * tfactors.T923i
                                      + -1.80122 + (5.0/3.0) * 0.21751 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ar38_to_p_Cl37_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar38_to_p_Cl37_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar38_to_He4_S34_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ar38 --> He4 + S34

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // cn18 
        ln_set_rate =  253.654 + -101.617 * tfactors.T9i + 635.332 * tfactors.T913i + -902.05 * tfactors.T913
                             + 42.5623 * tfactors.T9 + -2.08434 * tfactors.T953 + 478.85 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  101.617 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 635.332 * tfactors.T943i + (1.0/3.0) * -902.05 * tfactors.T923i
                                      + 42.5623 + (5.0/3.0) * -2.08434 * tfactors.T923 + 478.85 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ar38_to_He4_S34_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar38_to_He4_S34_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K37_to_p_Ar36_reaclib(T& rate_eval, const tf_t& tfactors) {

        // K37 --> p + Ar36

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10r
        ln_set_rate =  26.976 + -25.1829 * tfactors.T9i;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  25.1829 * tfactors.T9i * tfactors.T9i;

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
        ln_set_rate =  25.4611 + -30.5381 * tfactors.T9i + 7.50962 * tfactors.T913
                             + -1.68675 * tfactors.T9 + 0.189725 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  30.5381 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 7.50962 * tfactors.T923i
                                      + -1.68675 + (5.0/3.0) * 0.189725 * tfactors.T923;

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
        ln_set_rate =  43.7706 + -21.5578 * tfactors.T9i + -28.9682 * tfactors.T913i
                             + -1.0 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  21.5578 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -28.9682 * tfactors.T943i
                                      + (5.0/3.0) * -1.0 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_K37_to_p_Ar36_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K37_to_p_Ar36_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K37_to_He4_Cl33_reaclib(T& rate_eval, const tf_t& tfactors) {

        // K37 --> He4 + Cl33

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  71.6063 + -72.1995 * tfactors.T9i + -68.1442 * tfactors.T913i + -2.72746 * tfactors.T913
                             + -0.317537 * tfactors.T9 + 0.0162782 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  72.1995 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -68.1442 * tfactors.T943i + (1.0/3.0) * -2.72746 * tfactors.T923i
                                      + -0.317537 + (5.0/3.0) * 0.0162782 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_K37_to_He4_Cl33_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K37_to_He4_Cl33_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K39_to_p_Ar38_reaclib(T& rate_eval, const tf_t& tfactors) {

        // K39 --> p + Ar38

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  57.5639 + -74.0533 * tfactors.T9i + -29.0021 * tfactors.T913i + -0.525968 * tfactors.T913
                             + -1.94216 * tfactors.T9 + 0.267346 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  74.0533 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -29.0021 * tfactors.T943i + (1.0/3.0) * -0.525968 * tfactors.T923i
                                      + -1.94216 + (5.0/3.0) * 0.267346 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_K39_to_p_Ar38_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K39_to_p_Ar38_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K39_to_He4_Cl35_reaclib(T& rate_eval, const tf_t& tfactors) {

        // K39 --> He4 + Cl35

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  77.6477 + -83.7658 * tfactors.T9i + -68.2848 * tfactors.T913i + 0.0178545 * tfactors.T913
                             + -2.06783 * tfactors.T9 + 0.199659 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  83.7658 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -68.2848 * tfactors.T943i + (1.0/3.0) * 0.0178545 * tfactors.T923i
                                      + -2.06783 + (5.0/3.0) * 0.199659 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_K39_to_He4_Cl35_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K39_to_He4_Cl35_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca40_to_p_K39_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ca40 --> p + K39

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // lo18r
        ln_set_rate =  2786.44 + -101.871 * tfactors.T9i + 802.18 * tfactors.T913i + -4010.27 * tfactors.T913
                             + 1137.69 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  101.871 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 802.18 * tfactors.T943i + (1.0/3.0) * -4010.27 * tfactors.T923i
                                      + 1137.69 * tfactors.T9i;

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
        ln_set_rate =  613.153 + -109.213 * tfactors.T9i + 641.844 * tfactors.T913i + -1248.49 * tfactors.T913
                             + 566.426 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  109.213 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 641.844 * tfactors.T943i + (1.0/3.0) * -1248.49 * tfactors.T923i
                                      + 566.426 * tfactors.T9i;

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
        ln_set_rate =  127.306 + -98.3134 * tfactors.T9i + 41.1723 * tfactors.T913i + -149.299 * tfactors.T913
                             + 10.5229 * tfactors.T9 + -0.68208 * tfactors.T953 + 60.7367 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  98.3134 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 41.1723 * tfactors.T943i + (1.0/3.0) * -149.299 * tfactors.T923i
                                      + 10.5229 + (5.0/3.0) * -0.68208 * tfactors.T923 + 60.7367 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ca40_to_p_K39_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca40_to_p_K39_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca40_to_He4_Ar36_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ca40 --> He4 + Ar36

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  77.2826 + -81.6916 * tfactors.T9i + -71.0046 * tfactors.T913i + 4.0656 * tfactors.T913
                             + -5.26509 * tfactors.T9 + 0.683546 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  81.6916 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -71.0046 * tfactors.T943i + (1.0/3.0) * 4.0656 * tfactors.T923i
                                      + -5.26509 + (5.0/3.0) * 0.683546 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ca40_to_He4_Ar36_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca40_to_He4_Ar36_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca42_to_He4_Ar38_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ca42 --> He4 + Ar38

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  77.5813 + -72.6145 * tfactors.T9i + -71.1296 * tfactors.T913i + 2.75299 * tfactors.T913
                             + -2.43406 * tfactors.T9 + 0.199511 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  72.6145 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -71.1296 * tfactors.T943i + (1.0/3.0) * 2.75299 * tfactors.T923i
                                      + -2.43406 + (5.0/3.0) * 0.199511 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ca42_to_He4_Ar38_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca42_to_He4_Ar38_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc41_to_p_Ca40_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Sc41 --> p + Ca40

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10r
        ln_set_rate =  72.2762 + -46.1048 * tfactors.T9i + -33.7287 * tfactors.T913
                             + 5.22449 * tfactors.T9 + -0.357104 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  46.1048 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -33.7287 * tfactors.T923i
                                      + 5.22449 + (5.0/3.0) * -0.357104 * tfactors.T923;

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
        ln_set_rate =  41.5969 + -12.5925 * tfactors.T9i + -31.1711 * tfactors.T913i
                             + -1.0 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  12.5925 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -31.1711 * tfactors.T943i
                                      + (5.0/3.0) * -1.0 * tfactors.T923 + 0.833333 * tfactors.T9i;

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
        ln_set_rate =  27.2189 + -19.9148 * tfactors.T9i;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  19.9148 * tfactors.T9i * tfactors.T9i;

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

        rate_eval.screened_rates(k_Sc41_to_p_Ca40_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc41_to_p_Ca40_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc41_to_He4_K37_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Sc41 --> He4 + K37

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  74.8459 + -72.7267 * tfactors.T9i + -73.6776 * tfactors.T913i + -5.3733 * tfactors.T913
                             + 0.685776 * tfactors.T9 + -0.123735 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  72.7267 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -73.6776 * tfactors.T943i + (1.0/3.0) * -5.3733 * tfactors.T923i
                                      + 0.685776 + (5.0/3.0) * -0.123735 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Sc41_to_He4_K37_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc41_to_He4_K37_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc43_to_p_Ca42_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Sc43 --> p + Ca42

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  57.4102 + -57.2082 * tfactors.T9i + -31.1381 * tfactors.T913i + -0.151667 * tfactors.T913
                             + -1.76926 * tfactors.T9 + 0.219569 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  57.2082 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -31.1381 * tfactors.T943i + (1.0/3.0) * -0.151667 * tfactors.T923i
                                      + -1.76926 + (5.0/3.0) * 0.219569 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Sc43_to_p_Ca42_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc43_to_p_Ca42_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc43_to_He4_K39_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Sc43 --> He4 + K39

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  78.3727 + -55.7693 * tfactors.T9i + -73.8006 * tfactors.T913i + 1.87885 * tfactors.T913
                             + -2.75862 * tfactors.T9 + 0.279678 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  55.7693 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -73.8006 * tfactors.T943i + (1.0/3.0) * 1.87885 * tfactors.T923i
                                      + -2.75862 + (5.0/3.0) * 0.279678 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Sc43_to_He4_K39_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc43_to_He4_K39_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti44_to_p_Sc43_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ti44 --> p + Sc43

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  62.5939 + -100.373 * tfactors.T9i + -32.1734 * tfactors.T913i + -1.77078 * tfactors.T913
                             + -2.21706 * tfactors.T9 + 0.298499 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  100.373 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -32.1734 * tfactors.T943i + (1.0/3.0) * -1.77078 * tfactors.T923i
                                      + -2.21706 + (5.0/3.0) * 0.298499 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ti44_to_p_Sc43_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti44_to_p_Sc43_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti44_to_He4_Ca40_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ti44 --> He4 + Ca40

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // chw0 
        ln_set_rate =  78.6991 + -59.4974 * tfactors.T9i + -76.4273 * tfactors.T913i + 3.87451 * tfactors.T913
                             + -3.61477 * tfactors.T9 + 0.367451 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  59.4974 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -76.4273 * tfactors.T943i + (1.0/3.0) * 3.87451 * tfactors.T923i
                                      + -3.61477 + (5.0/3.0) * 0.367451 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ti44_to_He4_Ca40_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti44_to_He4_Ca40_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti46_to_p_Sc45_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ti46 --> p + Sc45

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  60.5905 + -120.044 * tfactors.T9i + -32.1843 * tfactors.T913i + 1.38642 * tfactors.T913
                             + -1.64811 * tfactors.T9 + 0.157323 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  120.044 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -32.1843 * tfactors.T943i + (1.0/3.0) * 1.38642 * tfactors.T923i
                                      + -1.64811 + (5.0/3.0) * 0.157323 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ti46_to_p_Sc45_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti46_to_p_Sc45_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti46_to_He4_Ca42_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ti46 --> He4 + Ca42

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  79.7886 + -92.8507 * tfactors.T9i + -76.5379 * tfactors.T913i + 2.80099 * tfactors.T913
                             + -2.34788 * tfactors.T9 + 0.150014 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  92.8507 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -76.5379 * tfactors.T943i + (1.0/3.0) * 2.80099 * tfactors.T923i
                                      + -2.34788 + (5.0/3.0) * 0.150014 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ti46_to_He4_Ca42_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti46_to_He4_Ca42_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V47_to_p_Ti46_reaclib(T& rate_eval, const tf_t& tfactors) {

        // V47 --> p + Ti46

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  58.0555 + -59.9674 * tfactors.T9i + -33.2034 * tfactors.T913i + 0.505619 * tfactors.T913
                             + -1.73654 * tfactors.T9 + 0.207342 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  59.9674 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -33.2034 * tfactors.T943i + (1.0/3.0) * 0.505619 * tfactors.T923i
                                      + -1.73654 + (5.0/3.0) * 0.207342 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_V47_to_p_Ti46_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V47_to_p_Ti46_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V47_to_He4_Sc43_reaclib(T& rate_eval, const tf_t& tfactors) {

        // V47 --> He4 + Sc43

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  84.6713 + -95.6099 * tfactors.T9i + -79.122 * tfactors.T913i + -7.07006 * tfactors.T913
                             + 0.424183 * tfactors.T9 + -0.0665231 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  95.6099 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -79.122 * tfactors.T943i + (1.0/3.0) * -7.07006 * tfactors.T923i
                                      + 0.424183 + (5.0/3.0) * -0.0665231 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_V47_to_He4_Sc43_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V47_to_He4_Sc43_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V49_to_p_Ti48_reaclib(T& rate_eval, const tf_t& tfactors) {

        // V49 --> p + Ti48

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  56.8861 + -78.4254 * tfactors.T9i + -33.2133 * tfactors.T913i + 0.564373 * tfactors.T913
                             + -1.18789 * tfactors.T9 + 0.110892 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  78.4254 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -33.2133 * tfactors.T943i + (1.0/3.0) * 0.564373 * tfactors.T923i
                                      + -1.18789 + (5.0/3.0) * 0.110892 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_V49_to_p_Ti48_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V49_to_p_Ti48_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V49_to_He4_Sc45_reaclib(T& rate_eval, const tf_t& tfactors) {

        // V49 --> He4 + Sc45

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  85.107 + -108.085 * tfactors.T9i + -79.222 * tfactors.T913i + -5.31022 * tfactors.T913
                             + -1.37323 * tfactors.T9 + 0.21679 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  108.085 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -79.222 * tfactors.T943i + (1.0/3.0) * -5.31022 * tfactors.T923i
                                      + -1.37323 + (5.0/3.0) * 0.21679 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_V49_to_He4_Sc45_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V49_to_He4_Sc45_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr48_to_p_V47_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cr48 --> p + V47

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // nfisn
        ln_set_rate =  48.892 + -93.8243 * tfactors.T9i + -45.9868 * tfactors.T913i + 13.6822 * tfactors.T913
                             + -0.376902 * tfactors.T9 + -0.0194875 * tfactors.T953 + -6.92325 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  93.8243 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -45.9868 * tfactors.T943i + (1.0/3.0) * 13.6822 * tfactors.T923i
                                      + -0.376902 + (5.0/3.0) * -0.0194875 * tfactors.T923 + -6.92325 * tfactors.T9i;

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
        ln_set_rate =  65.6231 + -94.5854 * tfactors.T9i + -110.655 * tfactors.T913i + 83.0232 * tfactors.T913
                             + -19.7762 * tfactors.T9 + 3.03961 * tfactors.T953 + -47.9742 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  94.5854 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -110.655 * tfactors.T943i + (1.0/3.0) * 83.0232 * tfactors.T923i
                                      + -19.7762 + (5.0/3.0) * 3.03961 * tfactors.T923 + -47.9742 * tfactors.T9i;

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
        ln_set_rate =  67.7403 + -100.13 * tfactors.T9i + -34.0548 * tfactors.T913i + -3.41973 * tfactors.T913
                             + 1.16501 * tfactors.T9 + -0.105543 * tfactors.T953 + -6.20886 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  100.13 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -34.0548 * tfactors.T943i + (1.0/3.0) * -3.41973 * tfactors.T923i
                                      + 1.16501 + (5.0/3.0) * -0.105543 * tfactors.T923 + -6.20886 * tfactors.T9i;

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
        ln_set_rate =  536.523 + -99.3659 * tfactors.T9i + 317.171 * tfactors.T913i + -911.679 * tfactors.T913
                             + 94.4245 * tfactors.T9 + -10.1973 * tfactors.T953 + 332.227 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  99.3659 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 317.171 * tfactors.T943i + (1.0/3.0) * -911.679 * tfactors.T923i
                                      + 94.4245 + (5.0/3.0) * -10.1973 * tfactors.T923 + 332.227 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cr48_to_p_V47_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr48_to_p_V47_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr48_to_He4_Ti44_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cr48 --> He4 + Ti44

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  89.7573 + -89.3041 * tfactors.T9i + -81.667 * tfactors.T913i + -10.6333 * tfactors.T913
                             + -0.672613 * tfactors.T9 + 0.161209 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  89.3041 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -81.667 * tfactors.T943i + (1.0/3.0) * -10.6333 * tfactors.T923i
                                      + -0.672613 + (5.0/3.0) * 0.161209 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cr48_to_He4_Ti44_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr48_to_He4_Ti44_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr50_to_p_V49_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cr50 --> p + V49

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  60.9884 + -111.305 * tfactors.T9i + -34.217 * tfactors.T913i + 1.28258 * tfactors.T913
                             + -1.5098 * tfactors.T9 + 0.142011 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  111.305 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -34.217 * tfactors.T943i + (1.0/3.0) * 1.28258 * tfactors.T923i
                                      + -1.5098 + (5.0/3.0) * 0.142011 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cr50_to_p_V49_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr50_to_p_V49_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr50_to_He4_Ti46_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cr50 --> He4 + Ti46

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  82.0042 + -99.346 * tfactors.T9i + -81.7658 * tfactors.T913i + 3.48637 * tfactors.T913
                             + -3.9188 * tfactors.T9 + 0.440356 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  99.346 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -81.7658 * tfactors.T943i + (1.0/3.0) * 3.48637 * tfactors.T923i
                                      + -3.9188 + (5.0/3.0) * 0.440356 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cr50_to_He4_Ti46_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr50_to_He4_Ti46_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr52_to_p_V51_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cr52 --> p + V51

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  59.9066 + -121.899 * tfactors.T9i + -34.2261 * tfactors.T913i + 3.30454 * tfactors.T913
                             + -2.61654 * tfactors.T9 + 0.283602 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  121.899 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -34.2261 * tfactors.T943i + (1.0/3.0) * 3.30454 * tfactors.T923i
                                      + -2.61654 + (5.0/3.0) * 0.283602 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cr52_to_p_V51_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr52_to_p_V51_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr52_to_He4_Ti48_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cr52 --> He4 + Ti48

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  81.6507 + -108.55 * tfactors.T9i + -81.8569 * tfactors.T913i + 4.92305 * tfactors.T913
                             + -5.04112 * tfactors.T9 + 0.6175 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  108.55 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -81.8569 * tfactors.T943i + (1.0/3.0) * 4.92305 * tfactors.T923i
                                      + -5.04112 + (5.0/3.0) * 0.6175 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cr52_to_He4_Ti48_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr52_to_He4_Ti48_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn51_to_p_Cr50_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Mn51 --> p + Cr50

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  57.7375 + -61.1651 * tfactors.T9i + -35.2065 * tfactors.T913i + 0.884621 * tfactors.T913
                             + -1.76421 * tfactors.T9 + 0.202494 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  61.1651 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -35.2065 * tfactors.T943i + (1.0/3.0) * 0.884621 * tfactors.T923i
                                      + -1.76421 + (5.0/3.0) * 0.202494 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Mn51_to_p_Cr50_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn51_to_p_Cr50_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn51_to_He4_V47_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Mn51 --> He4 + V47

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  81.4259 + -100.544 * tfactors.T9i + -84.2732 * tfactors.T913i + -2.98061 * tfactors.T913
                             + -0.531361 * tfactors.T9 + 0.023612 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  100.544 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -84.2732 * tfactors.T943i + (1.0/3.0) * -2.98061 * tfactors.T923i
                                      + -0.531361 + (5.0/3.0) * 0.023612 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Mn51_to_He4_V47_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn51_to_He4_V47_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn53_to_p_Cr52_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Mn53 --> p + Cr52

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  57.2932 + -76.1248 * tfactors.T9i + -35.2154 * tfactors.T913i + 0.73019 * tfactors.T913
                             + -1.26018 * tfactors.T9 + 0.121141 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  76.1248 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -35.2154 * tfactors.T943i + (1.0/3.0) * 0.73019 * tfactors.T923i
                                      + -1.26018 + (5.0/3.0) * 0.121141 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Mn53_to_p_Cr52_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn53_to_p_Cr52_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn53_to_He4_V49_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Mn53 --> He4 + V49

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  87.4299 + -106.249 * tfactors.T9i + -84.3633 * tfactors.T913i + -5.28933 * tfactors.T913
                             + -1.3103 * tfactors.T9 + 0.187088 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  106.249 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -84.3633 * tfactors.T943i + (1.0/3.0) * -5.28933 * tfactors.T923i
                                      + -1.3103 + (5.0/3.0) * 0.187088 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Mn53_to_He4_V49_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn53_to_He4_V49_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn55_to_p_Cr54_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Mn55 --> p + Cr54

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  56.8219 + -93.6136 * tfactors.T9i + -35.2237 * tfactors.T913i + 2.60768 * tfactors.T913
                             + -1.7864 * tfactors.T9 + 0.134189 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  93.6136 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -35.2237 * tfactors.T943i + (1.0/3.0) * 2.60768 * tfactors.T923i
                                      + -1.7864 + (5.0/3.0) * 0.134189 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Mn55_to_p_Cr54_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn55_to_p_Cr54_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn55_to_He4_V51_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Mn55 --> He4 + V51

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  86.8796 + -92.0714 * tfactors.T9i + -84.4467 * tfactors.T913i + -1.29542 * tfactors.T913
                             + -3.56909 * tfactors.T9 + 0.513926 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  92.0714 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -84.4467 * tfactors.T943i + (1.0/3.0) * -1.29542 * tfactors.T923i
                                      + -3.56909 + (5.0/3.0) * 0.513926 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Mn55_to_He4_V51_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn55_to_He4_V51_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe52_to_p_Mn51_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Fe52 --> p + Mn51

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  61.728 + -85.6325 * tfactors.T9i + -36.1825 * tfactors.T913i + 0.873042 * tfactors.T913
                             + -2.89731 * tfactors.T9 + 0.364394 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  85.6325 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -36.1825 * tfactors.T943i + (1.0/3.0) * 0.873042 * tfactors.T923i
                                      + -2.89731 + (5.0/3.0) * 0.364394 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Fe52_to_p_Mn51_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe52_to_p_Mn51_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe52_to_He4_Cr48_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Fe52 --> He4 + Cr48

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  90.1474 + -92.109 * tfactors.T9i + -86.7459 * tfactors.T913i + -9.79373 * tfactors.T913
                             + -0.772169 * tfactors.T9 + 0.155883 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  92.109 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -86.7459 * tfactors.T943i + (1.0/3.0) * -9.79373 * tfactors.T923i
                                      + -0.772169 + (5.0/3.0) * 0.155883 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Fe52_to_He4_Cr48_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe52_to_He4_Cr48_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe54_to_p_Mn53_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Fe54 --> p + Mn53

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  61.3432 + -102.741 * tfactors.T9i + -36.1913 * tfactors.T913i + 1.44056 * tfactors.T913
                             + -1.50734 * tfactors.T9 + 0.142124 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  102.741 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -36.1913 * tfactors.T943i + (1.0/3.0) * 1.44056 * tfactors.T923i
                                      + -1.50734 + (5.0/3.0) * 0.142124 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Fe54_to_p_Mn53_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe54_to_p_Mn53_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe54_to_He4_Cr50_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Fe54 --> He4 + Cr50

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  84.6959 + -97.6852 * tfactors.T9i + -86.8349 * tfactors.T913i + 2.37343 * tfactors.T913
                             + -3.65677 * tfactors.T9 + 0.412128 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  97.6852 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -86.8349 * tfactors.T943i + (1.0/3.0) * 2.37343 * tfactors.T923i
                                      + -3.65677 + (5.0/3.0) * 0.412128 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Fe54_to_He4_Cr50_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe54_to_He4_Cr50_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe56_to_p_Mn55_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Fe56 --> p + Mn55

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  61.6593 + -118.177 * tfactors.T9i + -36.1995 * tfactors.T913i + 1.15346 * tfactors.T913
                             + -1.94437 * tfactors.T9 + 0.218429 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  118.177 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -36.1995 * tfactors.T943i + (1.0/3.0) * 1.15346 * tfactors.T923i
                                      + -1.94437 + (5.0/3.0) * 0.218429 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Fe56_to_p_Mn55_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe56_to_p_Mn55_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe56_to_He4_Cr52_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Fe56 --> He4 + Cr52

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  82.9955 + -88.349 * tfactors.T9i + -86.9175 * tfactors.T913i + 5.655 * tfactors.T913
                             + -4.53625 * tfactors.T9 + 0.479806 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  88.349 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -86.9175 * tfactors.T943i + (1.0/3.0) * 5.655 * tfactors.T923i
                                      + -4.53625 + (5.0/3.0) * 0.479806 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Fe56_to_He4_Cr52_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe56_to_He4_Cr52_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe58_to_He4_Cr54_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Fe58 --> He4 + Cr54

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  88.2775 + -88.7256 * tfactors.T9i + -86.9945 * tfactors.T913i + -1.05674 * tfactors.T913
                             + -4.13443 * tfactors.T9 + 0.599141 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  88.7256 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -86.9945 * tfactors.T943i + (1.0/3.0) * -1.05674 * tfactors.T923i
                                      + -4.13443 + (5.0/3.0) * 0.599141 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Fe58_to_He4_Cr54_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe58_to_He4_Cr54_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co55_to_p_Fe54_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Co55 --> p + Fe54

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  57.8292 + -58.766 * tfactors.T9i + -37.1544 * tfactors.T913i + 0.950364 * tfactors.T913
                             + -1.77529 * tfactors.T9 + 0.198562 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  58.766 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -37.1544 * tfactors.T943i + (1.0/3.0) * 0.950364 * tfactors.T923i
                                      + -1.77529 + (5.0/3.0) * 0.198562 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Co55_to_p_Fe54_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co55_to_p_Fe54_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co55_to_He4_Mn51_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Co55 --> He4 + Mn51

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  90.613 + -95.2861 * tfactors.T9i + -89.274 * tfactors.T913i + -10.4373 * tfactors.T913
                             + 1.00492 * tfactors.T9 + -0.125548 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  95.2861 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -89.274 * tfactors.T943i + (1.0/3.0) * -10.4373 * tfactors.T923i
                                      + 1.00492 + (5.0/3.0) * -0.125548 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Co55_to_He4_Mn51_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co55_to_He4_Mn51_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co56_to_n_Co55_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Co56 --> n + Co55

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  40.626 + -117.009 * tfactors.T9i + -1.86357 * tfactors.T913
                             + 0.616591 * tfactors.T9 + -0.0839313 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  117.009 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.86357 * tfactors.T923i
                                      + 0.616591 + (5.0/3.0) * -0.0839313 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Co56_to_n_Co55_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co56_to_n_Co55_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co57_to_n_Co56_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Co57 --> n + Co56

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  41.1522 + -132.015 * tfactors.T9i + -1.37855 * tfactors.T913
                             + 0.299896 * tfactors.T9 + -0.04382 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  132.015 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.37855 * tfactors.T923i
                                      + 0.299896 + (5.0/3.0) * -0.04382 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Co57_to_n_Co56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co57_to_n_Co56_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co57_to_p_Fe56_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Co57 --> p + Fe56

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  57.6663 + -69.9498 * tfactors.T9i + -37.1625 * tfactors.T913i + 1.06776 * tfactors.T913
                             + -1.31689 * tfactors.T9 + 0.122089 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  69.9498 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -37.1625 * tfactors.T943i + (1.0/3.0) * 1.06776 * tfactors.T923i
                                      + -1.31689 + (5.0/3.0) * 0.122089 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Co57_to_p_Fe56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co57_to_p_Fe56_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co57_to_He4_Mn53_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Co57 --> He4 + Mn53

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  87.0081 + -82.174 * tfactors.T9i + -89.3558 * tfactors.T913i + 0.000882861 * tfactors.T913
                             + -2.79327 * tfactors.T9 + 0.309057 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  82.174 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -89.3558 * tfactors.T943i + (1.0/3.0) * 0.000882861 * tfactors.T923i
                                      + -2.79327 + (5.0/3.0) * 0.309057 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Co57_to_He4_Mn53_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co57_to_He4_Mn53_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co59_to_p_Fe58_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Co59 --> p + Fe58

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  57.8189 + -85.4551 * tfactors.T9i + -37.17 * tfactors.T913i + 1.63285 * tfactors.T913
                             + -1.16636 * tfactors.T9 + 0.0508619 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  85.4551 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -37.17 * tfactors.T943i + (1.0/3.0) * 1.63285 * tfactors.T923i
                                      + -1.16636 + (5.0/3.0) * 0.0508619 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Co59_to_p_Fe58_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co59_to_p_Fe58_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co59_to_He4_Mn55_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Co59 --> He4 + Mn55

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  84.7829 + -80.5671 * tfactors.T9i + -89.432 * tfactors.T913i + 5.33428 * tfactors.T913
                             + -5.11461 * tfactors.T9 + 0.613808 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  80.5671 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -89.432 * tfactors.T943i + (1.0/3.0) * 5.33428 * tfactors.T923i
                                      + -5.11461 + (5.0/3.0) * 0.613808 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Co59_to_He4_Mn55_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co59_to_He4_Mn55_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni56_to_p_Co55_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ni56 --> p + Co55

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  63.1318 + -83.1473 * tfactors.T9i + -38.1053 * tfactors.T913i + -0.210947 * tfactors.T913
                             + -2.68377 * tfactors.T9 + 0.355814 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  83.1473 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -38.1053 * tfactors.T943i + (1.0/3.0) * -0.210947 * tfactors.T923i
                                      + -2.68377 + (5.0/3.0) * 0.355814 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ni56_to_p_Co55_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni56_to_p_Co55_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni56_to_He4_Fe52_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ni56 --> He4 + Fe52

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  91.6226 + -92.801 * tfactors.T9i + -91.6819 * tfactors.T913i + -9.51885 * tfactors.T913
                             + -0.533014 * tfactors.T9 + 0.0892607 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  92.801 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -91.6819 * tfactors.T943i + (1.0/3.0) * -9.51885 * tfactors.T923i
                                      + -0.533014 + (5.0/3.0) * 0.0892607 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ni56_to_He4_Fe52_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni56_to_He4_Fe52_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni58_to_p_Co57_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ni58 --> p + Co57

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  61.775 + -94.8373 * tfactors.T9i + -38.1133 * tfactors.T913i + 1.77414 * tfactors.T913
                             + -1.48268 * tfactors.T9 + 0.121073 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  94.8373 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -38.1133 * tfactors.T943i + (1.0/3.0) * 1.77414 * tfactors.T923i
                                      + -1.48268 + (5.0/3.0) * 0.121073 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ni58_to_p_Co57_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni58_to_p_Co57_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni58_to_He4_Fe54_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ni58 --> He4 + Fe54

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  85.2327 + -74.2705 * tfactors.T9i + -91.7628 * tfactors.T913i + 4.23027 * tfactors.T913
                             + -3.31305 * tfactors.T9 + 0.271293 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  74.2705 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -91.7628 * tfactors.T943i + (1.0/3.0) * 4.23027 * tfactors.T923i
                                      + -3.31305 + (5.0/3.0) * 0.271293 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ni58_to_He4_Fe54_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni58_to_He4_Fe54_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni60_to_p_Co59_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ni60 --> p + Co59

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  61.717 + -110.622 * tfactors.T9i + -38.1208 * tfactors.T913i + 2.42609 * tfactors.T913
                             + -1.5885 * tfactors.T9 + 0.110313 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  110.622 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -38.1208 * tfactors.T943i + (1.0/3.0) * 2.42609 * tfactors.T923i
                                      + -1.5885 + (5.0/3.0) * 0.110313 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ni60_to_p_Co59_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni60_to_p_Co59_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni60_to_He4_Fe56_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ni60 --> He4 + Fe56

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  85.6559 + -73.0113 * tfactors.T9i + -91.8383 * tfactors.T913i + 3.63237 * tfactors.T913
                             + -2.73367 * tfactors.T9 + 0.179684 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  73.0113 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -91.8383 * tfactors.T943i + (1.0/3.0) * 3.63237 * tfactors.T923i
                                      + -2.73367 + (5.0/3.0) * 0.179684 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ni60_to_He4_Fe56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni60_to_He4_Fe56_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu59_to_p_Ni58_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cu59 --> p + Ni58

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // si13n
        ln_set_rate =  58.96 + -39.67 * tfactors.T9i + -39.0526 * tfactors.T913i + 1.0436 * tfactors.T913
                             + -2.10834 * tfactors.T9 + 0.239559 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  39.67 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -39.0526 * tfactors.T943i + (1.0/3.0) * 1.0436 * tfactors.T923i
                                      + -2.10834 + (5.0/3.0) * 0.239559 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cu59_to_p_Ni58_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu59_to_p_Ni58_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu59_to_He4_Co55_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cu59 --> He4 + Co55

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  90.972 + -55.1745 * tfactors.T9i + -94.1404 * tfactors.T913i + -2.62786 * tfactors.T913
                             + -2.12066 * tfactors.T9 + 0.237999 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  55.1745 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -94.1404 * tfactors.T943i + (1.0/3.0) * -2.62786 * tfactors.T923i
                                      + -2.12066 + (5.0/3.0) * 0.237999 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cu59_to_He4_Co55_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu59_to_He4_Co55_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Zn60_to_p_Cu59_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Zn60 --> p + Cu59

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  62.2929 + -59.4094 * tfactors.T9i + -39.9808 * tfactors.T913i + 1.17311 * tfactors.T913
                             + -2.90486 * tfactors.T9 + 0.339644 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  59.4094 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -39.9808 * tfactors.T943i + (1.0/3.0) * 1.17311 * tfactors.T923i
                                      + -2.90486 + (5.0/3.0) * 0.339644 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Zn60_to_p_Cu59_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn60_to_p_Cu59_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Zn60_to_He4_Ni56_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Zn60 --> He4 + Ni56

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  86.0619 + -31.4367 * tfactors.T9i + -96.4898 * tfactors.T913i + 6.47209 * tfactors.T913
                             + -5.2029 * tfactors.T9 + 0.533391 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  31.4367 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -96.4898 * tfactors.T943i + (1.0/3.0) * 6.47209 * tfactors.T923i
                                      + -5.2029 + (5.0/3.0) * 0.533391 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Zn60_to_He4_Ni56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn60_to_He4_Ni56_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_C12_to_3He4_reaclib(T& rate_eval, const tf_t& tfactors) {

        // C12 --> 3 He4

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // fy05r
        ln_set_rate =  22.394 + -88.5493 * tfactors.T9i + -13.49 * tfactors.T913i + 21.4259 * tfactors.T913
                             + -1.34769 * tfactors.T9 + 0.0879816 * tfactors.T953 + -10.1653 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  88.5493 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -13.49 * tfactors.T943i + (1.0/3.0) * 21.4259 * tfactors.T923i
                                      + -1.34769 + (5.0/3.0) * 0.0879816 * tfactors.T923 + -10.1653 * tfactors.T9i;

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
        ln_set_rate =  34.9561 + -85.4472 * tfactors.T9i + -23.57 * tfactors.T913i + 20.4886 * tfactors.T913
                             + -12.9882 * tfactors.T9 + -20.0 * tfactors.T953 + 0.83333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  85.4472 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -23.57 * tfactors.T943i + (1.0/3.0) * 20.4886 * tfactors.T923i
                                      + -12.9882 + (5.0/3.0) * -20.0 * tfactors.T923 + 0.83333 * tfactors.T9i;

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
        ln_set_rate =  45.7734 + -84.4227 * tfactors.T9i + -37.06 * tfactors.T913i + 29.3493 * tfactors.T913
                             + -115.507 * tfactors.T9 + -10.0 * tfactors.T953 + 1.66667 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  84.4227 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -37.06 * tfactors.T943i + (1.0/3.0) * 29.3493 * tfactors.T923i
                                      + -115.507 + (5.0/3.0) * -10.0 * tfactors.T923 + 1.66667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_C12_to_3He4_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_C12_to_3He4_reaclib) = drate_dT;
        }
    }

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
    void rate_Si29_He4_to_S33_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Si29 + He4 --> S33

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  49.5657 + -59.5755 * tfactors.T913i + 1.06274 * tfactors.T913
                             + -3.07529 * tfactors.T9 + 0.372011 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -59.5755 * tfactors.T943i + (1.0/3.0) * 1.06274 * tfactors.T923i
                                      + -3.07529 + (5.0/3.0) * 0.372011 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Si29_He4_to_S33_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si29_He4_to_S33_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si30_n_to_Si31_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Si30 + n --> Si31

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  -9.19438 + 0.0685821 * tfactors.T9i + -12.3502 * tfactors.T913i + 34.6486 * tfactors.T913
                             + -1.89409 * tfactors.T9 + 0.0781979 * tfactors.T953 + -14.4601 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  -0.0685821 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -12.3502 * tfactors.T943i + (1.0/3.0) * 34.6486 * tfactors.T923i
                                      + -1.89409 + (5.0/3.0) * 0.0781979 * tfactors.T923 + -14.4601 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Si30_n_to_Si31_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si30_n_to_Si31_reaclib) = drate_dT;
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
    void rate_Si31_He4_to_S35_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Si31 + He4 --> S35

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // rath 
        ln_set_rate =  -73.2581 + -1.94072 * tfactors.T9i + -98.5176 * tfactors.T913i + 171.207 * tfactors.T913
                             + -14.3839 * tfactors.T9 + 0.960708 * tfactors.T953 + -61.9936 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  1.94072 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -98.5176 * tfactors.T943i + (1.0/3.0) * 171.207 * tfactors.T923i
                                      + -14.3839 + (5.0/3.0) * 0.960708 * tfactors.T923 + -61.9936 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Si31_He4_to_S35_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si31_He4_to_S35_reaclib) = drate_dT;
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
    void rate_P33_p_to_S34_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // P33 + p --> S34

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // rath 
        ln_set_rate =  173.854 + -3.84298 * tfactors.T9i + 161.691 * tfactors.T913i + -350.015 * tfactors.T913
                             + 24.3597 * tfactors.T9 + -1.56396 * tfactors.T953 + 149.963 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  3.84298 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 161.691 * tfactors.T943i + (1.0/3.0) * -350.015 * tfactors.T923i
                                      + 24.3597 + (5.0/3.0) * -1.56396 * tfactors.T923 + 149.963 * tfactors.T9i;

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

        rate_eval.screened_rates(k_P33_p_to_S34_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P33_p_to_S34_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P33_He4_to_Cl37_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // P33 + He4 --> Cl37

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // rath 
        ln_set_rate =  50.9267 + -3.38505 * tfactors.T9i + 21.2001 * tfactors.T913i + -84.0086 * tfactors.T913
                             + 1.39426 * tfactors.T9 + 0.0896178 * tfactors.T953 + 47.8178 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  3.38505 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 21.2001 * tfactors.T943i + (1.0/3.0) * -84.0086 * tfactors.T923i
                                      + 1.39426 + (5.0/3.0) * 0.0896178 * tfactors.T923 + 47.8178 * tfactors.T9i;

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

        rate_eval.screened_rates(k_P33_He4_to_Cl37_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P33_He4_to_Cl37_reaclib) = drate_dT;
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
    void rate_S33_He4_to_Ar37_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // S33 + He4 --> Ar37

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  49.9315 + -65.4446 * tfactors.T913i + 3.59607 * tfactors.T913
                             + -3.40501 * tfactors.T9 + 0.363961 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -65.4446 * tfactors.T943i + (1.0/3.0) * 3.59607 * tfactors.T923i
                                      + -3.40501 + (5.0/3.0) * 0.363961 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_S33_He4_to_Ar37_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S33_He4_to_Ar37_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S34_n_to_S35_reaclib(T& rate_eval, const tf_t& tfactors) {

        // S34 + n --> S35

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  41.4131 + -0.229262 * tfactors.T9i + 25.2182 * tfactors.T913i + -59.2704 * tfactors.T913
                             + 4.22256 * tfactors.T9 + -0.274247 * tfactors.T953 + 24.532 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  0.229262 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 25.2182 * tfactors.T943i + (1.0/3.0) * -59.2704 * tfactors.T923i
                                      + 4.22256 + (5.0/3.0) * -0.274247 * tfactors.T923 + 24.532 * tfactors.T9i;

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

        rate_eval.screened_rates(k_S34_n_to_S35_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S34_n_to_S35_reaclib) = drate_dT;
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
    void rate_Cl37_p_to_Ar38_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cl37 + p --> Ar38

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  34.0095 + -27.9113 * tfactors.T913i + 0.282028 * tfactors.T913
                             + -1.80122 * tfactors.T9 + 0.21751 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -27.9113 * tfactors.T943i + (1.0/3.0) * 0.282028 * tfactors.T923i
                                      + -1.80122 + (5.0/3.0) * 0.21751 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cl37_p_to_Ar38_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl37_p_to_Ar38_reaclib) = drate_dT;
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
    void rate_Sc45_p_to_Ti46_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Sc45 + p --> Ti46

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  34.8383 + -32.1843 * tfactors.T913i + 1.38642 * tfactors.T913
                             + -1.64811 * tfactors.T9 + 0.157323 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -32.1843 * tfactors.T943i + (1.0/3.0) * 1.38642 * tfactors.T923i
                                      + -1.64811 + (5.0/3.0) * 0.157323 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Sc45_p_to_Ti46_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc45_p_to_Ti46_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc45_He4_to_V49_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Sc45 + He4 --> V49

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  60.1427 + -79.222 * tfactors.T913i + -5.31022 * tfactors.T913
                             + -1.37323 * tfactors.T9 + 0.21679 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -79.222 * tfactors.T943i + (1.0/3.0) * -5.31022 * tfactors.T923i
                                      + -1.37323 + (5.0/3.0) * 0.21679 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Sc45_He4_to_V49_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc45_He4_to_V49_reaclib) = drate_dT;
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
    void rate_Ti48_p_to_V49_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ti48 + p --> V49

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  35.2907 + -33.2133 * tfactors.T913i + 0.564373 * tfactors.T913
                             + -1.18789 * tfactors.T9 + 0.110892 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -33.2133 * tfactors.T943i + (1.0/3.0) * 0.564373 * tfactors.T923i
                                      + -1.18789 + (5.0/3.0) * 0.110892 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ti48_p_to_V49_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti48_p_to_V49_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti48_He4_to_Cr52_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ti48 + He4 --> Cr52

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  56.6787 + -81.8569 * tfactors.T913i + 4.92305 * tfactors.T913
                             + -5.04112 * tfactors.T9 + 0.6175 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -81.8569 * tfactors.T943i + (1.0/3.0) * 4.92305 * tfactors.T923i
                                      + -5.04112 + (5.0/3.0) * 0.6175 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ti48_He4_to_Cr52_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti48_He4_to_Cr52_reaclib) = drate_dT;
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
    void rate_V51_p_to_Cr52_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // V51 + p --> Cr52

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  34.1505 + -34.2261 * tfactors.T913i + 3.30454 * tfactors.T913
                             + -2.61654 * tfactors.T9 + 0.283602 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -34.2261 * tfactors.T943i + (1.0/3.0) * 3.30454 * tfactors.T923i
                                      + -2.61654 + (5.0/3.0) * 0.283602 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_V51_p_to_Cr52_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V51_p_to_Cr52_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V51_He4_to_Mn55_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // V51 + He4 --> Mn55

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  61.6131 + -84.4467 * tfactors.T913i + -1.29542 * tfactors.T913
                             + -3.56909 * tfactors.T9 + 0.513926 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -84.4467 * tfactors.T943i + (1.0/3.0) * -1.29542 * tfactors.T923i
                                      + -3.56909 + (5.0/3.0) * 0.513926 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_V51_He4_to_Mn55_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V51_He4_to_Mn55_reaclib) = drate_dT;
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
    void rate_Cr54_p_to_Mn55_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cr54 + p --> Mn55

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  34.9354 + -35.2237 * tfactors.T913i + 2.60768 * tfactors.T913
                             + -1.7864 * tfactors.T9 + 0.134189 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -35.2237 * tfactors.T943i + (1.0/3.0) * 2.60768 * tfactors.T923i
                                      + -1.7864 + (5.0/3.0) * 0.134189 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cr54_p_to_Mn55_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr54_p_to_Mn55_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr54_He4_to_Fe58_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cr54 + He4 --> Fe58

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  63.2926 + -86.9945 * tfactors.T913i + -1.05674 * tfactors.T913
                             + -4.13443 * tfactors.T9 + 0.599141 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -86.9945 * tfactors.T943i + (1.0/3.0) * -1.05674 * tfactors.T923i
                                      + -4.13443 + (5.0/3.0) * 0.599141 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cr54_He4_to_Fe58_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr54_He4_to_Fe58_reaclib) = drate_dT;
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
    void rate_Mn55_p_to_Fe56_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Mn55 + p --> Fe56

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  36.1888 + -36.1995 * tfactors.T913i + 1.15346 * tfactors.T913
                             + -1.94437 * tfactors.T9 + 0.218429 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -36.1995 * tfactors.T943i + (1.0/3.0) * 1.15346 * tfactors.T923i
                                      + -1.94437 + (5.0/3.0) * 0.218429 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Mn55_p_to_Fe56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn55_p_to_Fe56_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn55_He4_to_Co59_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Mn55 + He4 --> Co59

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  60.0838 + -89.432 * tfactors.T913i + 5.33428 * tfactors.T913
                             + -5.11461 * tfactors.T9 + 0.613808 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -89.432 * tfactors.T943i + (1.0/3.0) * 5.33428 * tfactors.T923i
                                      + -5.11461 + (5.0/3.0) * 0.613808 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Mn55_He4_to_Co59_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn55_He4_to_Co59_reaclib) = drate_dT;
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
    void rate_Fe56_He4_to_Ni60_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Fe56 + He4 --> Ni60

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  60.6673 + -91.8383 * tfactors.T913i + 3.63237 * tfactors.T913
                             + -2.73367 * tfactors.T9 + 0.179684 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -91.8383 * tfactors.T943i + (1.0/3.0) * 3.63237 * tfactors.T923i
                                      + -2.73367 + (5.0/3.0) * 0.179684 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Fe56_He4_to_Ni60_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe56_He4_to_Ni60_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe58_p_to_Co59_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Fe58 + p --> Co59

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  36.2182 + -37.17 * tfactors.T913i + 1.63285 * tfactors.T913
                             + -1.16636 * tfactors.T9 + 0.0508619 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -37.17 * tfactors.T943i + (1.0/3.0) * 1.63285 * tfactors.T923i
                                      + -1.16636 + (5.0/3.0) * 0.0508619 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Fe58_p_to_Co59_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe58_p_to_Co59_reaclib) = drate_dT;
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
    void rate_Co59_p_to_Ni60_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Co59 + p --> Ni60

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  35.957 + -38.1208 * tfactors.T913i + 2.42609 * tfactors.T913
                             + -1.5885 * tfactors.T9 + 0.110313 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -38.1208 * tfactors.T943i + (1.0/3.0) * 2.42609 * tfactors.T923i
                                      + -1.5885 + (5.0/3.0) * 0.110313 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Co59_p_to_Ni60_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co59_p_to_Ni60_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni56_He4_to_Zn60_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ni56 + He4 --> Zn60

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  61.0733 + -96.4898 * tfactors.T913i + 6.47209 * tfactors.T913
                             + -5.2029 * tfactors.T9 + 0.533391 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -96.4898 * tfactors.T943i + (1.0/3.0) * 6.47209 * tfactors.T923i
                                      + -5.2029 + (5.0/3.0) * 0.533391 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ni56_He4_to_Zn60_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni56_He4_to_Zn60_reaclib) = drate_dT;
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
    void rate_Cu59_p_to_Zn60_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cu59 + p --> Zn60

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  37.226 + -39.9808 * tfactors.T913i + 1.17311 * tfactors.T913
                             + -2.90486 * tfactors.T9 + 0.339644 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -39.9808 * tfactors.T943i + (1.0/3.0) * 1.17311 * tfactors.T923i
                                      + -2.90486 + (5.0/3.0) * 0.339644 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cu59_p_to_Zn60_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu59_p_to_Zn60_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_C12_C12_to_n_Mg23_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // C12 + C12 --> n + Mg23

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // cf88r
        ln_set_rate =  -12.8056 + -30.1498 * tfactors.T9i + 11.4826 * tfactors.T913
                             + 1.82849 * tfactors.T9 + -0.34844 * tfactors.T953;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  30.1498 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 11.4826 * tfactors.T923i
                                      + 1.82849 + (5.0/3.0) * -0.34844 * tfactors.T923;

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

        rate_eval.screened_rates(k_C12_C12_to_n_Mg23_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_C12_C12_to_n_Mg23_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_O16_C12_to_n_Si27_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // O16 + C12 --> n + Si27

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // cf88r
        ln_set_rate =  -132.213 + -1.46479 * tfactors.T9i + -293.089 * tfactors.T913i + 414.404 * tfactors.T913
                             + -28.0562 * tfactors.T9 + 1.61807 * tfactors.T953 + -178.28 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  1.46479 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -293.089 * tfactors.T943i + (1.0/3.0) * 414.404 * tfactors.T923i
                                      + -28.0562 + (5.0/3.0) * 1.61807 * tfactors.T923 + -178.28 * tfactors.T9i;

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

        rate_eval.screened_rates(k_O16_C12_to_n_Si27_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_O16_C12_to_n_Si27_reaclib) = drate_dT;
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
    void rate_Ne20_He4_to_n_Mg23_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ne20 + He4 --> n + Mg23

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  17.9544 + -83.7215 * tfactors.T9i + 1.83199 * tfactors.T913
                             + -0.290485 * tfactors.T9 + 0.0242929 * tfactors.T953;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  83.7215 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 1.83199 * tfactors.T923i
                                      + -0.290485 + (5.0/3.0) * 0.0242929 * tfactors.T923;

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

        rate_eval.screened_rates(k_Ne20_He4_to_n_Mg23_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne20_He4_to_n_Mg23_reaclib) = drate_dT;
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
    void rate_Mg24_He4_to_n_Si27_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Mg24 + He4 --> n + Si27

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  18.6606 + -83.5022 * tfactors.T9i + 0.942981 * tfactors.T913
                             + -0.104624 * tfactors.T9 + 0.00723421 * tfactors.T953;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  83.5022 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 0.942981 * tfactors.T923i
                                      + -0.104624 + (5.0/3.0) * 0.00723421 * tfactors.T923;

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

        rate_eval.screened_rates(k_Mg24_He4_to_n_Si27_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg24_He4_to_n_Si27_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mg24_He4_to_p_Al27_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Mg24 + He4 --> p + Al27

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10r
        ln_set_rate =  -26.2862 + -19.5422 * tfactors.T9i + 5.18642 * tfactors.T913i + -34.7936 * tfactors.T913
                             + 168.225 * tfactors.T9 + -115.825 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  19.5422 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 5.18642 * tfactors.T943i + (1.0/3.0) * -34.7936 * tfactors.T923i
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

        // il10r
        ln_set_rate =  -6.44575 + -22.8216 * tfactors.T9i + 18.0416 * tfactors.T913
                             + -1.54137 * tfactors.T9 + 0.0847506 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  22.8216 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 18.0416 * tfactors.T923i
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

        // il10n
        ln_set_rate =  30.0397 + -18.5791 * tfactors.T9i + -26.4162 * tfactors.T913i
                             + -2.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  18.5791 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -26.4162 * tfactors.T943i
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

        rate_eval.screened_rates(k_Mg24_He4_to_p_Al27_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg24_He4_to_p_Al27_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mg24_He4_to_C12_O16_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Mg24 + He4 --> C12 + O16

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // cf88r
        ln_set_rate =  49.5738 + -78.202 * tfactors.T9i + -133.413 * tfactors.T913i + 50.1572 * tfactors.T913
                             + -3.15987 * tfactors.T9 + 0.0178251 * tfactors.T953 + -23.7027 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  78.202 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -133.413 * tfactors.T943i + (1.0/3.0) * 50.1572 * tfactors.T923i
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

        rate_eval.screened_rates(k_Mg24_He4_to_C12_O16_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg24_He4_to_C12_O16_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Al27_p_to_n_Si27_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Al27 + p --> n + Si27

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  17.84 + -64.9237 * tfactors.T9i + 2.10552 * tfactors.T913
                             + -0.318695 * tfactors.T9 + 0.0257496 * tfactors.T953;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  64.9237 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 2.10552 * tfactors.T923i
                                      + -0.318695 + (5.0/3.0) * 0.0257496 * tfactors.T923;

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

        rate_eval.screened_rates(k_Al27_p_to_n_Si27_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Al27_p_to_n_Si27_reaclib) = drate_dT;
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
    void rate_Al27_p_to_C12_O16_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Al27 + p --> C12 + O16

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // cf88r
        ln_set_rate =  68.9829 + -59.8017 * tfactors.T9i + -119.242 * tfactors.T913i + 13.3667 * tfactors.T913
                             + 0.295425 * tfactors.T9 + -0.267288 * tfactors.T953 + -9.91729 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  59.8017 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -119.242 * tfactors.T943i + (1.0/3.0) * 13.3667 * tfactors.T923i
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

        rate_eval.screened_rates(k_Al27_p_to_C12_O16_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Al27_p_to_C12_O16_reaclib) = drate_dT;
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
    void rate_Si28_He4_to_p_P31_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Si28 + He4 --> p + P31

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10r
        ln_set_rate =  -13.4595 + -24.112 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  24.112 * tfactors.T9i * tfactors.T9i
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
        ln_set_rate =  -11.4335 + -25.6606 * tfactors.T9i + 21.521 * tfactors.T913
                             + -1.90355 * tfactors.T9 + 0.092724 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  25.6606 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 21.521 * tfactors.T923i
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

        // il10n
        ln_set_rate =  60.3424 + -22.2348 * tfactors.T9i + -31.932 * tfactors.T913i + -77.0334 * tfactors.T913
                             + -43.6847 * tfactors.T9 + -4.28955 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  22.2348 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -31.932 * tfactors.T943i + (1.0/3.0) * -77.0334 * tfactors.T923i
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

        rate_eval.screened_rates(k_Si28_He4_to_p_P31_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si28_He4_to_p_P31_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si28_He4_to_C12_Ne20_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Si28 + He4 --> C12 + Ne20

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // rolfr
        ln_set_rate =  -307.762 + -186.722 * tfactors.T9i + 514.197 * tfactors.T913i + -200.896 * tfactors.T913
                             + -6.42713 * tfactors.T9 + 0.758256 * tfactors.T953 + 236.359 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  186.722 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 514.197 * tfactors.T943i + (1.0/3.0) * -200.896 * tfactors.T923i
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

        rate_eval.screened_rates(k_Si28_He4_to_C12_Ne20_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si28_He4_to_C12_Ne20_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si28_He4_to_O16_O16_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Si28 + He4 --> O16 + O16

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // cf88r
        ln_set_rate =  97.7904 + -111.595 * tfactors.T9i + -119.324 * tfactors.T913i + -32.2497 * tfactors.T913
                             + 1.46214 * tfactors.T9 + -0.200893 * tfactors.T953 + 13.2148 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  111.595 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -119.324 * tfactors.T943i + (1.0/3.0) * -32.2497 * tfactors.T923i
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

        rate_eval.screened_rates(k_Si28_He4_to_O16_O16_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si28_He4_to_O16_O16_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si29_p_to_n_P29_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Si29 + p --> n + P29

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  18.5351 + -66.4331 * tfactors.T9i + 1.48018 * tfactors.T913
                             + -0.177129 * tfactors.T9 + 0.0127163 * tfactors.T953;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  66.4331 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 1.48018 * tfactors.T923i
                                      + -0.177129 + (5.0/3.0) * 0.0127163 * tfactors.T923;

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

        rate_eval.screened_rates(k_Si29_p_to_n_P29_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si29_p_to_n_P29_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si29_He4_to_n_S32_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Si29 + He4 --> n + S32

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  -2.87932 + -17.7056 * tfactors.T9i + 9.48125 * tfactors.T913
                             + 0.4472 * tfactors.T9 + -0.119237 * tfactors.T953;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  17.7056 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 9.48125 * tfactors.T923i
                                      + 0.4472 + (5.0/3.0) * -0.119237 * tfactors.T923;

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

        rate_eval.screened_rates(k_Si29_He4_to_n_S32_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si29_He4_to_n_S32_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si30_p_to_He4_Al27_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Si30 + p --> He4 + Al27

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  50.5056 + -27.5284 * tfactors.T9i + -56.5351 * tfactors.T913i + 1.60477 * tfactors.T913
                             + -1.40594 * tfactors.T9 + 0.127353 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  27.5284 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -56.5351 * tfactors.T943i + (1.0/3.0) * 1.60477 * tfactors.T923i
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

        rate_eval.screened_rates(k_Si30_p_to_He4_Al27_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si30_p_to_He4_Al27_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si30_He4_to_n_S33_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Si30 + He4 --> n + S33

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  17.3228 + -40.5385 * tfactors.T9i + -0.725571 * tfactors.T913
                             + 0.430731 * tfactors.T9 + -0.0323663 * tfactors.T953;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  40.5385 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.725571 * tfactors.T923i
                                      + 0.430731 + (5.0/3.0) * -0.0323663 * tfactors.T923;

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

        rate_eval.screened_rates(k_Si30_He4_to_n_S33_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si30_He4_to_n_S33_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si30_He4_to_p_P33_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Si30 + He4 --> p + P33

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // rath 
        ln_set_rate =  81.0152 + -37.7526 * tfactors.T9i + 110.297 * tfactors.T913i + -193.863 * tfactors.T913
                             + 11.4977 * tfactors.T9 + -0.687235 * tfactors.T953 + 94.819 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  37.7526 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 110.297 * tfactors.T943i + (1.0/3.0) * -193.863 * tfactors.T923i
                                      + 11.4977 + (5.0/3.0) * -0.687235 * tfactors.T923 + 94.819 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Si30_He4_to_p_P33_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si30_He4_to_p_P33_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si31_p_to_n_P31_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Si31 + p --> n + P31

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // rath 
        ln_set_rate =  80.5195 + -3.3832 * tfactors.T9i + 119.352 * tfactors.T913i + -198.709 * tfactors.T913
                             + 10.5994 * tfactors.T9 + -0.583506 * tfactors.T953 + 100.693 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  3.3832 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 119.352 * tfactors.T943i + (1.0/3.0) * -198.709 * tfactors.T923i
                                      + 10.5994 + (5.0/3.0) * -0.583506 * tfactors.T923 + 100.693 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Si31_p_to_n_P31_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si31_p_to_n_P31_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si31_He4_to_n_S34_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Si31 + He4 --> n + S34

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // rath 
        ln_set_rate =  -5.10649 + -2.7438 * tfactors.T9i + -22.4045 * tfactors.T913i + 23.4647 * tfactors.T913
                             + -4.73922 * tfactors.T9 + 0.371237 * tfactors.T953 + 4.62837 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  2.7438 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -22.4045 * tfactors.T943i + (1.0/3.0) * 23.4647 * tfactors.T923i
                                      + -4.73922 + (5.0/3.0) * 0.371237 * tfactors.T923 + 4.62837 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Si31_He4_to_n_S34_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si31_He4_to_n_S34_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P29_n_to_p_Si29_reaclib(T& rate_eval, const tf_t& tfactors) {

        // P29 + n --> p + Si29

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  18.5351 + 1.48018 * tfactors.T913
                             + -0.177129 * tfactors.T9 + 0.0127163 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 1.48018 * tfactors.T923i
                                      + -0.177129 + (5.0/3.0) * 0.0127163 * tfactors.T923;

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

        rate_eval.screened_rates(k_P29_n_to_p_Si29_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P29_n_to_p_Si29_reaclib) = drate_dT;
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
    void rate_P31_n_to_p_Si31_reaclib(T& rate_eval, const tf_t& tfactors) {

        // P31 + n --> p + Si31

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // rath 
        ln_set_rate =  81.2126 + -11.6108 * tfactors.T9i + 119.352 * tfactors.T913i + -198.709 * tfactors.T913
                             + 10.5994 * tfactors.T9 + -0.583506 * tfactors.T953 + 100.693 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  11.6108 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 119.352 * tfactors.T943i + (1.0/3.0) * -198.709 * tfactors.T923i
                                      + 10.5994 + (5.0/3.0) * -0.583506 * tfactors.T923 + 100.693 * tfactors.T9i;

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

        rate_eval.screened_rates(k_P31_n_to_p_Si31_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P31_n_to_p_Si31_reaclib) = drate_dT;
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
    void rate_P31_p_to_C12_Ne20_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // P31 + p --> C12 + Ne20

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // rolfr
        ln_set_rate =  -266.452 + -156.019 * tfactors.T9i + 361.154 * tfactors.T913i + -92.643 * tfactors.T913
                             + -9.98738 * tfactors.T9 + 0.892737 * tfactors.T953 + 161.042 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  156.019 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 361.154 * tfactors.T943i + (1.0/3.0) * -92.643 * tfactors.T923i
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

        rate_eval.screened_rates(k_P31_p_to_C12_Ne20_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P31_p_to_C12_Ne20_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P31_p_to_O16_O16_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // P31 + p --> O16 + O16

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // cf88r
        ln_set_rate =  86.3501 + -88.8797 * tfactors.T9i + -145.844 * tfactors.T913i + 8.72612 * tfactors.T913
                             + -0.554035 * tfactors.T9 + -0.137562 * tfactors.T953 + -6.88807 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  88.8797 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -145.844 * tfactors.T943i + (1.0/3.0) * 8.72612 * tfactors.T923i
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

        rate_eval.screened_rates(k_P31_p_to_O16_O16_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P31_p_to_O16_O16_reaclib) = drate_dT;
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
    void rate_P33_p_to_n_S33_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // P33 + p --> n + S33

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  5.79458 + -6.19492 * tfactors.T9i + 8.02978 * tfactors.T913
                             + -0.479244 * tfactors.T9 + 0.0134889 * tfactors.T953;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  6.19492 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 8.02978 * tfactors.T923i
                                      + -0.479244 + (5.0/3.0) * 0.0134889 * tfactors.T923;

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

        rate_eval.screened_rates(k_P33_p_to_n_S33_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P33_p_to_n_S33_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P33_p_to_He4_Si30_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // P33 + p --> He4 + Si30

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // rath 
        ln_set_rate =  81.5654 + -3.41491 * tfactors.T9i + 110.297 * tfactors.T913i + -193.863 * tfactors.T913
                             + 11.4977 * tfactors.T9 + -0.687235 * tfactors.T953 + 94.819 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  3.41491 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 110.297 * tfactors.T943i + (1.0/3.0) * -193.863 * tfactors.T923i
                                      + 11.4977 + (5.0/3.0) * -0.687235 * tfactors.T923 + 94.819 * tfactors.T9i;

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

        rate_eval.screened_rates(k_P33_p_to_He4_Si30_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P33_p_to_He4_Si30_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S32_n_to_He4_Si29_reaclib(T& rate_eval, const tf_t& tfactors) {

        // S32 + n --> He4 + Si29

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  -0.947538 + 9.48125 * tfactors.T913
                             + 0.4472 * tfactors.T9 + -0.119237 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 9.48125 * tfactors.T923i
                                      + 0.4472 + (5.0/3.0) * -0.119237 * tfactors.T923;

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

        rate_eval.screened_rates(k_S32_n_to_He4_Si29_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S32_n_to_He4_Si29_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S32_p_to_He4_P29_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // S32 + p --> He4 + P29

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  52.2142 + -48.7275 * tfactors.T9i + -62.3802 * tfactors.T913i + 0.459085 * tfactors.T913
                             + -0.870169 * tfactors.T9 + 0.0631143 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  48.7275 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -62.3802 * tfactors.T943i + (1.0/3.0) * 0.459085 * tfactors.T923i
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

        rate_eval.screened_rates(k_S32_p_to_He4_P29_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S32_p_to_He4_P29_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S32_He4_to_p_Cl35_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // S32 + He4 --> p + Cl35

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10r
        ln_set_rate =  -57.395 + -22.1894 * tfactors.T9i + 25.5338 * tfactors.T913
                             + 6.45824 * tfactors.T9 + -0.950294 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  22.1894 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 25.5338 * tfactors.T923i
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
        ln_set_rate =  32.2544 + -21.6564 * tfactors.T9i + -30.9147 * tfactors.T913i + -1.2345 * tfactors.T913
                             + 22.5118 * tfactors.T9 + -33.0589 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  21.6564 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -30.9147 * tfactors.T943i + (1.0/3.0) * -1.2345 * tfactors.T923i
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
        ln_set_rate =  2.42563 + -27.6662 * tfactors.T9i + 5.33756 * tfactors.T913
                             + 1.64418 * tfactors.T9 + -0.246167 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  27.6662 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 5.33756 * tfactors.T923i
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

        // il10r
        ln_set_rate =  -0.877602 + -25.5914 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  25.5914 * tfactors.T9i * tfactors.T9i
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

        rate_eval.screened_rates(k_S32_He4_to_p_Cl35_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S32_He4_to_p_Cl35_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S33_n_to_p_P33_reaclib(T& rate_eval, const tf_t& tfactors) {

        // S33 + n --> p + P33

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  5.10143 + 8.02978 * tfactors.T913
                             + -0.479244 * tfactors.T9 + 0.0134889 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 8.02978 * tfactors.T923i
                                      + -0.479244 + (5.0/3.0) * 0.0134889 * tfactors.T923;

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

        rate_eval.screened_rates(k_S33_n_to_p_P33_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S33_n_to_p_P33_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S33_n_to_He4_Si30_reaclib(T& rate_eval, const tf_t& tfactors) {

        // S33 + n --> He4 + Si30

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  17.1798 + -0.725571 * tfactors.T913
                             + 0.430731 * tfactors.T9 + -0.0323663 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -0.725571 * tfactors.T923i
                                      + 0.430731 + (5.0/3.0) * -0.0323663 * tfactors.T923;

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

        rate_eval.screened_rates(k_S33_n_to_He4_Si30_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S33_n_to_He4_Si30_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S33_p_to_n_Cl33_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // S33 + p --> n + Cl33

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  18.5772 + -73.8616 * tfactors.T9i + 1.39633 * tfactors.T913
                             + -0.136457 * tfactors.T9 + 0.00585594 * tfactors.T953;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  73.8616 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 1.39633 * tfactors.T923i
                                      + -0.136457 + (5.0/3.0) * 0.00585594 * tfactors.T923;

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

        rate_eval.screened_rates(k_S33_p_to_n_Cl33_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S33_p_to_n_Cl33_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S33_He4_to_n_Ar36_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // S33 + He4 --> n + Ar36

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  0.398258 + -23.219 * tfactors.T9i + 5.44923 * tfactors.T913
                             + 1.17359 * tfactors.T9 + -0.175811 * tfactors.T953;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  23.219 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 5.44923 * tfactors.T923i
                                      + 1.17359 + (5.0/3.0) * -0.175811 * tfactors.T923;

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

        rate_eval.screened_rates(k_S33_He4_to_n_Ar36_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S33_He4_to_n_Ar36_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S34_n_to_He4_Si31_reaclib(T& rate_eval, const tf_t& tfactors) {

        // S34 + n --> He4 + Si31

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // rath 
        ln_set_rate =  -2.47246 + -18.2474 * tfactors.T9i + -22.4045 * tfactors.T913i + 23.4647 * tfactors.T913
                             + -4.73922 * tfactors.T9 + 0.371237 * tfactors.T953 + 4.62837 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  18.2474 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -22.4045 * tfactors.T943i + (1.0/3.0) * 23.4647 * tfactors.T923i
                                      + -4.73922 + (5.0/3.0) * 0.371237 * tfactors.T923 + 4.62837 * tfactors.T9i;

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

        rate_eval.screened_rates(k_S34_n_to_He4_Si31_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S34_n_to_He4_Si31_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S34_p_to_He4_P31_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // S34 + p --> He4 + P31

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  51.0402 + -7.27426 * tfactors.T9i + -62.5433 * tfactors.T913i + 1.69339 * tfactors.T913
                             + -1.27872 * tfactors.T9 + 0.106616 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  7.27426 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -62.5433 * tfactors.T943i + (1.0/3.0) * 1.69339 * tfactors.T923i
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

        rate_eval.screened_rates(k_S34_p_to_He4_P31_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S34_p_to_He4_P31_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S34_He4_to_n_Ar37_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // S34 + He4 --> n + Ar37

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // SM93 
        ln_set_rate =  18.5358 + -53.7288 * tfactors.T9i + -1.69894 * tfactors.T913
                             + 0.371187 * tfactors.T9 + -0.0457242 * tfactors.T953 + 1.0 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  53.7288 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.69894 * tfactors.T923i
                                      + 0.371187 + (5.0/3.0) * -0.0457242 * tfactors.T923 + 1.0 * tfactors.T9i;

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

        // SM93 
        ln_set_rate =  18.5358 + -53.7288 * tfactors.T9i + 0.522368 * tfactors.T913
                             + -1.03944 * tfactors.T9 + 0.142735 * tfactors.T953;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  53.7288 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 0.522368 * tfactors.T923i
                                      + -1.03944 + (5.0/3.0) * 0.142735 * tfactors.T923;

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

        rate_eval.screened_rates(k_S34_He4_to_n_Ar37_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S34_He4_to_n_Ar37_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S34_He4_to_p_Cl37_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // S34 + He4 --> p + Cl37

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // SM93 
        ln_set_rate =  80.8229 + -35.2081 * tfactors.T9i + -47.5185 * tfactors.T913i + -32.046 * tfactors.T913
                             + 4.43077 * tfactors.T9 + -0.368628 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  35.2081 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -47.5185 * tfactors.T943i + (1.0/3.0) * -32.046 * tfactors.T923i
                                      + 4.43077 + (5.0/3.0) * -0.368628 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_S34_He4_to_p_Cl37_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S34_He4_to_p_Cl37_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S35_p_to_n_Cl35_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // S35 + p --> n + Cl35

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  7.51079 + -7.13865 * tfactors.T9i + 6.51754 * tfactors.T913
                             + -0.295846 * tfactors.T9 + -0.00271436 * tfactors.T953;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  7.13865 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 6.51754 * tfactors.T923i
                                      + -0.295846 + (5.0/3.0) * -0.00271436 * tfactors.T923;

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

        rate_eval.screened_rates(k_S35_p_to_n_Cl35_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S35_p_to_n_Cl35_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S35_He4_to_n_Ar38_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // S35 + He4 --> n + Ar38

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  51.1473 + -65.58 * tfactors.T913i + 0.762185 * tfactors.T913
                             + -0.938779 * tfactors.T9 + 0.0668864 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -65.58 * tfactors.T943i + (1.0/3.0) * 0.762185 * tfactors.T923i
                                      + -0.938779 + (5.0/3.0) * 0.0668864 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_S35_He4_to_n_Ar38_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S35_He4_to_n_Ar38_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl33_n_to_p_S33_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cl33 + n --> p + S33

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  18.5772 + 1.39633 * tfactors.T913
                             + -0.136457 * tfactors.T9 + 0.00585594 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 1.39633 * tfactors.T923i
                                      + -0.136457 + (5.0/3.0) * 0.00585594 * tfactors.T923;

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

        rate_eval.screened_rates(k_Cl33_n_to_p_S33_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl33_n_to_p_S33_reaclib) = drate_dT;
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
    void rate_Cl35_n_to_p_S35_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cl35 + n --> p + S35

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  7.51079 + 6.51754 * tfactors.T913
                             + -0.295846 * tfactors.T9 + -0.00271436 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 6.51754 * tfactors.T923i
                                      + -0.295846 + (5.0/3.0) * -0.00271436 * tfactors.T923;

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

        rate_eval.screened_rates(k_Cl35_n_to_p_S35_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl35_n_to_p_S35_reaclib) = drate_dT;
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
    void rate_Cl37_p_to_n_Ar37_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cl37 + p --> n + Ar37

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  16.2094 + -18.5233 * tfactors.T9i + 1.64495 * tfactors.T913
                             + -0.154677 * tfactors.T9 + 0.0083347 * tfactors.T953;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  18.5233 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 1.64495 * tfactors.T923i
                                      + -0.154677 + (5.0/3.0) * 0.0083347 * tfactors.T923;

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

        rate_eval.screened_rates(k_Cl37_p_to_n_Ar37_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl37_p_to_n_Ar37_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl37_p_to_He4_S34_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cl37 + p --> He4 + S34

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // SM93 
        ln_set_rate =  80.6961 + -47.5185 * tfactors.T913i + -32.046 * tfactors.T913
                             + 4.43077 * tfactors.T9 + -0.368628 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -47.5185 * tfactors.T943i + (1.0/3.0) * -32.046 * tfactors.T923i
                                      + 4.43077 + (5.0/3.0) * -0.368628 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cl37_p_to_He4_S34_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl37_p_to_He4_S34_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar36_n_to_He4_S33_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ar36 + n --> He4 + S33

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  3.04033 + 5.44923 * tfactors.T913
                             + 1.17359 * tfactors.T9 + -0.175811 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 5.44923 * tfactors.T923i
                                      + 1.17359 + (5.0/3.0) * -0.175811 * tfactors.T923;

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

        rate_eval.screened_rates(k_Ar36_n_to_He4_S33_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar36_n_to_He4_S33_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar36_p_to_He4_Cl33_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ar36 + p --> He4 + Cl33

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  54.8009 + -50.6426 * tfactors.T9i + -68.1442 * tfactors.T913i + 0.291238 * tfactors.T913
                             + -0.791384 * tfactors.T9 + 0.0524823 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  50.6426 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -68.1442 * tfactors.T943i + (1.0/3.0) * 0.291238 * tfactors.T923i
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

        rate_eval.screened_rates(k_Ar36_p_to_He4_Cl33_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar36_p_to_He4_Cl33_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar36_He4_to_p_K39_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ar36 + He4 --> p + K39

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  20.6367 + -14.9533 * tfactors.T9i + -30.0732 * tfactors.T913i + 7.03263 * tfactors.T913
                             + -1.10085 * tfactors.T9 + 0.133768 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  14.9533 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -30.0732 * tfactors.T943i + (1.0/3.0) * 7.03263 * tfactors.T923i
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

        rate_eval.screened_rates(k_Ar36_He4_to_p_K39_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar36_He4_to_p_K39_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar37_n_to_p_Cl37_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ar37 + n --> p + Cl37

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  16.2094 + 1.64495 * tfactors.T913
                             + -0.154677 * tfactors.T9 + 0.0083347 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 1.64495 * tfactors.T923i
                                      + -0.154677 + (5.0/3.0) * 0.0083347 * tfactors.T923;

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

        rate_eval.screened_rates(k_Ar37_n_to_p_Cl37_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar37_n_to_p_Cl37_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar37_n_to_He4_S34_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ar37 + n --> He4 + S34

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // SM93 
        ln_set_rate =  18.409 + -1.69894 * tfactors.T913
                             + 0.371187 * tfactors.T9 + -0.0457242 * tfactors.T953 + 1.0 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -1.69894 * tfactors.T923i
                                      + 0.371187 + (5.0/3.0) * -0.0457242 * tfactors.T923 + 1.0 * tfactors.T9i;

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

        // SM93 
        ln_set_rate =  18.409 + 0.522368 * tfactors.T913
                             + -1.03944 * tfactors.T9 + 0.142735 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 0.522368 * tfactors.T923i
                                      + -1.03944 + (5.0/3.0) * 0.142735 * tfactors.T923;

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

        rate_eval.screened_rates(k_Ar37_n_to_He4_S34_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar37_n_to_He4_S34_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar37_p_to_n_K37_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ar37 + p --> n + K37

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  19.3703 + -80.417 * tfactors.T9i + 0.399131 * tfactors.T913
                             + 0.0819317 * tfactors.T9 + -0.0133489 * tfactors.T953;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  80.417 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 0.399131 * tfactors.T923i
                                      + 0.0819317 + (5.0/3.0) * -0.0133489 * tfactors.T923;

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

        rate_eval.screened_rates(k_Ar37_p_to_n_K37_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar37_p_to_n_K37_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar37_He4_to_n_Ca40_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ar37 + He4 --> n + Ca40

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  -12.0914 + -20.2822 * tfactors.T9i + 13.8882 * tfactors.T913
                             + 0.260223 * tfactors.T9 + -0.108063 * tfactors.T953;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  20.2822 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 13.8882 * tfactors.T923i
                                      + 0.260223 + (5.0/3.0) * -0.108063 * tfactors.T923;

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

        rate_eval.screened_rates(k_Ar37_He4_to_n_Ca40_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar37_He4_to_n_Ca40_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar38_n_to_He4_S35_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ar38 + n --> He4 + S35

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  53.7965 + -2.57381 * tfactors.T9i + -65.58 * tfactors.T913i + 0.762185 * tfactors.T913
                             + -0.938779 * tfactors.T9 + 0.0668864 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  2.57381 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -65.58 * tfactors.T943i + (1.0/3.0) * 0.762185 * tfactors.T923i
                                      + -0.938779 + (5.0/3.0) * 0.0668864 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ar38_n_to_He4_S35_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar38_n_to_He4_S35_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar38_p_to_He4_Cl35_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ar38 + p --> He4 + Cl35

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  53.7764 + -9.71246 * tfactors.T9i + -68.2848 * tfactors.T913i + 2.5993 * tfactors.T913
                             + -1.59144 * tfactors.T9 + 0.137745 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  9.71246 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -68.2848 * tfactors.T943i + (1.0/3.0) * 2.5993 * tfactors.T923i
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

        rate_eval.screened_rates(k_Ar38_p_to_He4_Cl35_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar38_p_to_He4_Cl35_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K37_n_to_p_Ar37_reaclib(T& rate_eval, const tf_t& tfactors) {

        // K37 + n --> p + Ar37

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  19.3703 + 0.399131 * tfactors.T913
                             + 0.0819317 * tfactors.T9 + -0.0133489 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 0.399131 * tfactors.T923i
                                      + 0.0819317 + (5.0/3.0) * -0.0133489 * tfactors.T923;

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

        rate_eval.screened_rates(k_K37_n_to_p_Ar37_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K37_n_to_p_Ar37_reaclib) = drate_dT;
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
    void rate_K39_He4_to_p_Ca42_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // K39 + He4 --> p + Ca42

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  -80.8734 + -1.43889 * tfactors.T9i + -31.1381 * tfactors.T913i + 117.868 * tfactors.T913
                             + -27.0455 * tfactors.T9 + 2.75661 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  1.43889 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -31.1381 * tfactors.T943i + (1.0/3.0) * 117.868 * tfactors.T923i
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

        rate_eval.screened_rates(k_K39_He4_to_p_Ca42_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K39_He4_to_p_Ca42_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca40_n_to_He4_Ar37_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ca40 + n --> He4 + Ar37

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  -9.43576 + 13.8882 * tfactors.T913
                             + 0.260223 * tfactors.T9 + -0.108063 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * 13.8882 * tfactors.T923i
                                      + 0.260223 + (5.0/3.0) * -0.108063 * tfactors.T923;

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

        rate_eval.screened_rates(k_Ca40_n_to_He4_Ar37_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca40_n_to_He4_Ar37_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca40_p_to_He4_K37_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ca40 + p --> He4 + K37

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  56.9522 + -60.1348 * tfactors.T9i + -73.6776 * tfactors.T913i + 0.493294 * tfactors.T913
                             + -0.902307 * tfactors.T9 + 0.0641391 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  60.1348 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -73.6776 * tfactors.T943i + (1.0/3.0) * 0.493294 * tfactors.T923i
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

        rate_eval.screened_rates(k_Ca40_p_to_He4_K37_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca40_p_to_He4_K37_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca40_He4_to_p_Sc43_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ca40 + He4 --> p + Sc43

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  35.6575 + -40.8757 * tfactors.T9i + -32.1734 * tfactors.T913i + 0.0296879 * tfactors.T913
                             + -0.95232 * tfactors.T9 + 0.129022 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  40.8757 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -32.1734 * tfactors.T943i + (1.0/3.0) * 0.0296879 * tfactors.T923i
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

        rate_eval.screened_rates(k_Ca40_He4_to_p_Sc43_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca40_He4_to_p_Sc43_reaclib) = drate_dT;
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
    void rate_Ca42_He4_to_p_Sc45_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ca42 + He4 --> p + Sc45

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  18.625 + -27.1933 * tfactors.T9i + -32.1843 * tfactors.T913i + 14.0733 * tfactors.T913
                             + -3.38786 * tfactors.T9 + 0.349494 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  27.1933 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -32.1843 * tfactors.T943i + (1.0/3.0) * 14.0733 * tfactors.T923i
                                      + -3.38786 + (5.0/3.0) * 0.349494 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ca42_He4_to_p_Sc45_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca42_He4_to_p_Sc45_reaclib) = drate_dT;
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
    void rate_Sc45_p_to_He4_Ca42_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Sc45 + p --> He4 + Ca42

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  17.8284 + -32.1843 * tfactors.T913i + 14.0733 * tfactors.T913
                             + -3.38786 * tfactors.T9 + 0.349494 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -32.1843 * tfactors.T943i + (1.0/3.0) * 14.0733 * tfactors.T923i
                                      + -3.38786 + (5.0/3.0) * 0.349494 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Sc45_p_to_He4_Ca42_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc45_p_to_He4_Ca42_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc45_He4_to_p_Ti48_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Sc45 + He4 --> p + Ti48

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  57.207 + -79.222 * tfactors.T913i + 0.0710184 * tfactors.T913
                             + -1.47173 * tfactors.T9 + 0.15944 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -79.222 * tfactors.T943i + (1.0/3.0) * 0.0710184 * tfactors.T923i
                                      + -1.47173 + (5.0/3.0) * 0.15944 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Sc45_He4_to_p_Ti48_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc45_He4_to_p_Ti48_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti44_p_to_He4_Sc41_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ti44 + p --> He4 + Sc41

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  60.8451 + -46.9055 * tfactors.T9i + -79.013 * tfactors.T913i + -0.989013 * tfactors.T913
                             + -0.463943 * tfactors.T9 + 0.00902542 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  46.9055 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -79.013 * tfactors.T943i + (1.0/3.0) * -0.989013 * tfactors.T923i
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

        rate_eval.screened_rates(k_Ti44_p_to_He4_Sc41_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti44_p_to_He4_Sc41_reaclib) = drate_dT;
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
    void rate_Ti46_p_to_He4_Sc43_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ti46 + p --> He4 + Sc43

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  60.7369 + -35.6425 * tfactors.T9i + -79.122 * tfactors.T913i + -0.495539 * tfactors.T913
                             + -0.618414 * tfactors.T9 + 0.0175873 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  35.6425 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -79.122 * tfactors.T943i + (1.0/3.0) * -0.495539 * tfactors.T923i
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

        rate_eval.screened_rates(k_Ti46_p_to_He4_Sc43_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti46_p_to_He4_Sc43_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti46_He4_to_p_V49_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ti46 + He4 --> p + V49

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  -19.3117 + -11.9588 * tfactors.T9i + -34.217 * tfactors.T913i + 47.2584 * tfactors.T913
                             + -9.18871 * tfactors.T9 + 0.87823 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  11.9588 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -34.217 * tfactors.T943i + (1.0/3.0) * 47.2584 * tfactors.T923i
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

        rate_eval.screened_rates(k_Ti46_He4_to_p_V49_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti46_He4_to_p_V49_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti48_p_to_He4_Sc45_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ti48 + p --> He4 + Sc45

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  60.5759 + -29.6597 * tfactors.T9i + -79.222 * tfactors.T913i + 0.0710184 * tfactors.T913
                             + -1.47173 * tfactors.T9 + 0.15944 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  29.6597 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -79.222 * tfactors.T943i + (1.0/3.0) * 0.0710184 * tfactors.T923i
                                      + -1.47173 + (5.0/3.0) * 0.15944 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ti48_p_to_He4_Sc45_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti48_p_to_He4_Sc45_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti48_He4_to_p_V51_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ti48 + He4 --> p + V51

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  -12.9707 + -13.3498 * tfactors.T9i + -34.2261 * tfactors.T913i + 38.6574 * tfactors.T913
                             + -7.18743 * tfactors.T9 + 0.692255 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  13.3498 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -34.2261 * tfactors.T943i + (1.0/3.0) * 38.6574 * tfactors.T923i
                                      + -7.18743 + (5.0/3.0) * 0.692255 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ti48_He4_to_p_V51_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti48_He4_to_p_V51_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V47_p_to_He4_Ti44_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // V47 + p --> He4 + Ti44

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // chw0r
        ln_set_rate =  -76.6143 + -6.02945 * tfactors.T9i + 70.2835 * tfactors.T913
                             + -7.99061 * tfactors.T9 + 0.486213 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  6.02945 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 70.2835 * tfactors.T923i
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

        rate_eval.screened_rates(k_V47_p_to_He4_Ti44_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V47_p_to_He4_Ti44_reaclib) = drate_dT;
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
    void rate_V51_p_to_He4_Ti48_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // V51 + p --> He4 + Ti48

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  -13.7548 + -34.2261 * tfactors.T913i + 38.6574 * tfactors.T913
                             + -7.18743 * tfactors.T9 + 0.692255 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -34.2261 * tfactors.T943i + (1.0/3.0) * 38.6574 * tfactors.T923i
                                      + -7.18743 + (5.0/3.0) * 0.692255 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_V51_p_to_He4_Ti48_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V51_p_to_He4_Ti48_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V51_He4_to_p_Cr54_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // V51 + He4 --> p + Cr54

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  -98.078 + -1.54216 * tfactors.T9i + -35.2237 * tfactors.T913i + 131.726 * tfactors.T913
                             + -29.7395 * tfactors.T9 + 3.03809 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  1.54216 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -35.2237 * tfactors.T943i + (1.0/3.0) * 131.726 * tfactors.T923i
                                      + -29.7395 + (5.0/3.0) * 3.03809 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_V51_He4_to_p_Cr54_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V51_He4_to_p_Cr54_reaclib) = drate_dT;
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
    void rate_Cr50_p_to_He4_V47_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cr50 + p --> He4 + V47

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  62.5904 + -39.3786 * tfactors.T9i + -84.2732 * tfactors.T913i + -1.10027 * tfactors.T913
                             + -0.487084 * tfactors.T9 + 0.00453269 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  39.3786 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -84.2732 * tfactors.T943i + (1.0/3.0) * -1.10027 * tfactors.T923i
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

        rate_eval.screened_rates(k_Cr50_p_to_He4_V47_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr50_p_to_He4_V47_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr50_He4_to_p_Mn53_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cr50 + He4 --> p + Mn53

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  -67.5619 + -5.0555 * tfactors.T9i + -36.1913 * tfactors.T913i + 99.0857 * tfactors.T913
                             + -20.8943 * tfactors.T9 + 2.06209 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  5.0555 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -36.1913 * tfactors.T943i + (1.0/3.0) * 99.0857 * tfactors.T923i
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

        rate_eval.screened_rates(k_Cr50_He4_to_p_Mn53_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr50_He4_to_p_Mn53_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr52_p_to_He4_V49_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cr52 + p --> He4 + V49

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  60.8046 + -30.1245 * tfactors.T9i + -84.3633 * tfactors.T913i + 2.73419 * tfactors.T913
                             + -2.05655 * tfactors.T9 + 0.203834 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  30.1245 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -84.3633 * tfactors.T943i + (1.0/3.0) * 2.73419 * tfactors.T923i
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

        rate_eval.screened_rates(k_Cr52_p_to_He4_V49_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr52_p_to_He4_V49_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr52_He4_to_p_Mn55_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cr52 + He4 --> p + Mn55

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  18.2793 + -29.8283 * tfactors.T9i + -36.1995 * tfactors.T913i + 11.023 * tfactors.T913
                             + -2.435 * tfactors.T9 + 0.306933 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  29.8283 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -36.1995 * tfactors.T943i + (1.0/3.0) * 11.023 * tfactors.T923i
                                      + -2.435 + (5.0/3.0) * 0.306933 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cr52_He4_to_p_Mn55_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr52_He4_to_p_Mn55_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr54_p_to_He4_V51_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Cr54 + p --> He4 + V51

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  -94.698 + -35.2237 * tfactors.T913i + 131.726 * tfactors.T913
                             + -29.7395 * tfactors.T9 + 3.03809 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -35.2237 * tfactors.T943i + (1.0/3.0) * 131.726 * tfactors.T923i
                                      + -29.7395 + (5.0/3.0) * 3.03809 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cr54_p_to_He4_V51_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr54_p_to_He4_V51_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn51_p_to_He4_Cr48_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Mn51 + p --> He4 + Cr48

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  58.7312 + -6.47654 * tfactors.T9i + -86.7459 * tfactors.T913i + 1.05653 * tfactors.T913
                             + -1.15757 * tfactors.T9 + 0.0877546 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  6.47654 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -86.7459 * tfactors.T943i + (1.0/3.0) * 1.05653 * tfactors.T923i
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

        rate_eval.screened_rates(k_Mn51_p_to_He4_Cr48_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn51_p_to_He4_Cr48_reaclib) = drate_dT;
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
    void rate_Mn53_He4_to_n_Co56_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Mn53 + He4 --> n + Co56

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  11.3813 + -49.8407 * tfactors.T9i + -2.81854 * tfactors.T913
                             + 1.84188 * tfactors.T9 + -0.164509 * tfactors.T953;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  49.8407 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -2.81854 * tfactors.T923i
                                      + 1.84188 + (5.0/3.0) * -0.164509 * tfactors.T923;

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

        rate_eval.screened_rates(k_Mn53_He4_to_n_Co56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn53_He4_to_n_Co56_reaclib) = drate_dT;
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
    void rate_Mn55_p_to_He4_Cr52_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Mn55 + p --> He4 + Cr52

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  17.7897 + -36.1995 * tfactors.T913i + 11.023 * tfactors.T913
                             + -2.435 * tfactors.T9 + 0.306933 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -36.1995 * tfactors.T943i + (1.0/3.0) * 11.023 * tfactors.T923i
                                      + -2.435 + (5.0/3.0) * 0.306933 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Mn55_p_to_He4_Cr52_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn55_p_to_He4_Cr52_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn55_He4_to_p_Fe58_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Mn55 + He4 --> p + Fe58

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  -77.6133 + -4.88801 * tfactors.T9i + -37.17 * tfactors.T913i + 110.878 * tfactors.T913
                             + -24.5555 * tfactors.T9 + 2.47003 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  4.88801 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -37.17 * tfactors.T943i + (1.0/3.0) * 110.878 * tfactors.T923i
                                      + -24.5555 + (5.0/3.0) * 2.47003 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Mn55_He4_to_p_Fe58_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn55_He4_to_p_Fe58_reaclib) = drate_dT;
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
    void rate_Fe54_p_to_He4_Mn51_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Fe54 + p --> He4 + Mn51

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  65.37 + -36.5201 * tfactors.T9i + -89.274 * tfactors.T913i + -0.862452 * tfactors.T913
                             + -0.635672 * tfactors.T9 + 0.0196464 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  36.5201 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -89.274 * tfactors.T943i + (1.0/3.0) * -0.862452 * tfactors.T923i
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

        rate_eval.screened_rates(k_Fe54_p_to_He4_Mn51_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe54_p_to_He4_Mn51_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe54_He4_to_p_Co57_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Fe54 + He4 --> p + Co57

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  -1.42195 + -20.5668 * tfactors.T9i + -38.1133 * tfactors.T913i + 29.3541 * tfactors.T913
                             + -4.75966 * tfactors.T9 + 0.40418 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  20.5668 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -38.1133 * tfactors.T943i + (1.0/3.0) * 29.3541 * tfactors.T923i
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

        rate_eval.screened_rates(k_Fe54_He4_to_p_Co57_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe54_He4_to_p_Co57_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe56_p_to_n_Co56_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Fe56 + p --> n + Co56

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  22.6511 + -62.0649 * tfactors.T9i + -1.13331 * tfactors.T913
                             + 0.347185 * tfactors.T9 + -0.0328879 * tfactors.T953;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  62.0649 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.13331 * tfactors.T923i
                                      + 0.347185 + (5.0/3.0) * -0.0328879 * tfactors.T923;

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

        rate_eval.screened_rates(k_Fe56_p_to_n_Co56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe56_p_to_n_Co56_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe56_p_to_He4_Mn53_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Fe56 + p --> He4 + Mn53

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  57.4007 + -12.2241 * tfactors.T9i + -89.3558 * tfactors.T913i + 9.25287 * tfactors.T913
                             + -3.70352 * tfactors.T9 + 0.372202 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  12.2241 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -89.3558 * tfactors.T943i + (1.0/3.0) * 9.25287 * tfactors.T923i
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

        rate_eval.screened_rates(k_Fe56_p_to_He4_Mn53_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe56_p_to_He4_Mn53_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe56_He4_to_p_Co59_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Fe56 + He4 --> p + Co59

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  20.2414 + -37.6102 * tfactors.T9i + -38.1208 * tfactors.T913i + 13.9086 * tfactors.T913
                             + -3.45838 * tfactors.T9 + 0.391302 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  37.6102 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -38.1208 * tfactors.T943i + (1.0/3.0) * 13.9086 * tfactors.T923i
                                      + -3.45838 + (5.0/3.0) * 0.391302 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Fe56_He4_to_p_Co59_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe56_He4_to_p_Co59_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe58_p_to_He4_Mn55_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Fe58 + p --> He4 + Mn55

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  -74.5149 + -37.17 * tfactors.T913i + 110.878 * tfactors.T913
                             + -24.5555 * tfactors.T9 + 2.47003 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -37.17 * tfactors.T943i + (1.0/3.0) * 110.878 * tfactors.T923i
                                      + -24.5555 + (5.0/3.0) * 2.47003 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Fe58_p_to_He4_Mn55_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe58_p_to_He4_Mn55_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co55_p_to_He4_Fe52_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Co55 + p --> He4 + Fe52

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  61.4434 + -9.65363 * tfactors.T9i + -91.6819 * tfactors.T913i + -0.329235 * tfactors.T913
                             + -0.780924 * tfactors.T9 + 0.0425179 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  9.65363 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -91.6819 * tfactors.T943i + (1.0/3.0) * -0.329235 * tfactors.T923i
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

        rate_eval.screened_rates(k_Co55_p_to_He4_Fe52_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co55_p_to_He4_Fe52_reaclib) = drate_dT;
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
    void rate_Co56_p_to_n_Ni56_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Co56 + p --> n + Ni56

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  13.3721 + -33.8622 * tfactors.T9i + 1.76846 * tfactors.T913
                             + 0.197992 * tfactors.T9 + -0.017494 * tfactors.T953;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  33.8622 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 1.76846 * tfactors.T923i
                                      + 0.197992 + (5.0/3.0) * -0.017494 * tfactors.T923;

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

        rate_eval.screened_rates(k_Co56_p_to_n_Ni56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co56_p_to_n_Ni56_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co56_He4_to_n_Cu59_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Co56 + He4 --> n + Cu59

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  8.89185 + -61.8349 * tfactors.T9i + -1.00882 * tfactors.T913
                             + 1.88993 * tfactors.T9 + -0.202319 * tfactors.T953;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  61.8349 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.00882 * tfactors.T923i
                                      + 1.88993 + (5.0/3.0) * -0.202319 * tfactors.T923;

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

        rate_eval.screened_rates(k_Co56_He4_to_n_Cu59_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co56_He4_to_n_Cu59_reaclib) = drate_dT;
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
    void rate_Co57_He4_to_p_Ni60_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Co57 + He4 --> p + Ni60

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  49.6113 + -94.2152 * tfactors.T913i + 16.1097 * tfactors.T913
                             + -4.84624 * tfactors.T9 + 0.441487 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -94.2152 * tfactors.T943i + (1.0/3.0) * 16.1097 * tfactors.T923i
                                      + -4.84624 + (5.0/3.0) * 0.441487 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Co57_He4_to_p_Ni60_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co57_He4_to_p_Ni60_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co59_p_to_He4_Fe56_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Co59 + p --> He4 + Fe56

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  19.47 + -38.1208 * tfactors.T913i + 13.9086 * tfactors.T913
                             + -3.45838 * tfactors.T9 + 0.391302 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + -(1.0/3.0) * -38.1208 * tfactors.T943i + (1.0/3.0) * 13.9086 * tfactors.T923i
                                      + -3.45838 + (5.0/3.0) * 0.391302 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Co59_p_to_He4_Fe56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co59_p_to_He4_Fe56_reaclib) = drate_dT;
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
    void rate_Ni56_He4_to_p_Cu59_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ni56 + He4 --> p + Cu59

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  20.8471 + -27.9728 * tfactors.T9i + -39.9808 * tfactors.T913i + 13.6773 * tfactors.T913
                             + -3.76429 * tfactors.T9 + 0.438096 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  27.9728 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -39.9808 * tfactors.T943i + (1.0/3.0) * 13.6773 * tfactors.T923i
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

        rate_eval.screened_rates(k_Ni56_He4_to_p_Cu59_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni56_He4_to_p_Cu59_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni58_p_to_He4_Co55_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ni58 + p --> He4 + Co55

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  63.6142 + -15.5045 * tfactors.T9i + -94.1404 * tfactors.T913i + 3.39179 * tfactors.T913
                             + -1.71062 * tfactors.T9 + 0.133003 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  15.5045 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -94.1404 * tfactors.T943i + (1.0/3.0) * 3.39179 * tfactors.T923i
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

        rate_eval.screened_rates(k_Ni58_p_to_He4_Co55_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni58_p_to_He4_Co55_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni60_p_to_He4_Co57_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ni60 + p --> He4 + Co57

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  53.0001 + -3.0615 * tfactors.T9i + -94.2152 * tfactors.T913i + 16.1097 * tfactors.T913
                             + -4.84624 * tfactors.T9 + 0.441487 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  3.0615 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -94.2152 * tfactors.T943i + (1.0/3.0) * 16.1097 * tfactors.T923i
                                      + -4.84624 + (5.0/3.0) * 0.441487 * tfactors.T923 + -0.666667 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ni60_p_to_He4_Co57_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni60_p_to_He4_Co57_reaclib) = drate_dT;
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
    void rate_Ni57_to_n_Ni56_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ni57 --> n + Ni56

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  38.3694 + -118.941 * tfactors.T9i + -1.19665 * tfactors.T913
                             + 0.507179 * tfactors.T9 + -0.074604 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  118.941 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.19665 * tfactors.T923i
                                      + 0.507179 + (5.0/3.0) * -0.074604 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ni57_to_n_Ni56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni57_to_n_Ni56_reaclib) = drate_dT;
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
    void rate_Ni58_to_n_Ni57_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ni58 --> n + Ni57

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  42.3391 + -141.772 * tfactors.T9i + -1.90814 * tfactors.T913
                             + 0.493188 * tfactors.T9 + -0.0684633 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  141.772 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.90814 * tfactors.T923i
                                      + 0.493188 + (5.0/3.0) * -0.0684633 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ni58_to_n_Ni57_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni58_to_n_Ni57_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni58_n_to_Ni59_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ni58 + n --> Ni59

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  8.63197 + 0.13279 * tfactors.T9i + -11.785 * tfactors.T913i + 19.5347 * tfactors.T913
                             + -0.857179 * tfactors.T9 + 0.00111653 * tfactors.T953 + -9.35642 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  -0.13279 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -11.785 * tfactors.T943i + (1.0/3.0) * 19.5347 * tfactors.T923i
                                      + -0.857179 + (5.0/3.0) * 0.00111653 * tfactors.T923 + -9.35642 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ni58_n_to_Ni59_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni58_n_to_Ni59_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni59_to_n_Ni58_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ni59 --> n + Ni58

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  30.9258 + -104.3 * tfactors.T9i + -11.785 * tfactors.T913i + 19.5347 * tfactors.T913
                             + -0.857179 * tfactors.T9 + 0.00111653 * tfactors.T953 + -7.85642 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  104.3 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -11.785 * tfactors.T943i + (1.0/3.0) * 19.5347 * tfactors.T923i
                                      + -0.857179 + (5.0/3.0) * 0.00111653 * tfactors.T923 + -7.85642 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ni59_to_n_Ni58_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni59_to_n_Ni58_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni59_n_to_Ni60_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ni59 + n --> Ni60

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  25.5762 + -0.110079 * tfactors.T9i + 9.92411 * tfactors.T913i + -20.503 * tfactors.T913
                             + 1.39981 * tfactors.T9 + -0.101641 * tfactors.T953 + 8.58204 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  0.110079 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 9.92411 * tfactors.T943i + (1.0/3.0) * -20.503 * tfactors.T923i
                                      + 1.39981 + (5.0/3.0) * -0.101641 * tfactors.T923 + 8.58204 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ni59_n_to_Ni60_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni59_n_to_Ni60_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni60_to_n_Ni59_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ni60 --> n + Ni59

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  50.6431 + -132.266 * tfactors.T9i + 9.92411 * tfactors.T913i + -20.503 * tfactors.T913
                             + 1.39981 * tfactors.T9 + -0.101641 * tfactors.T953 + 10.082 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  132.266 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 9.92411 * tfactors.T943i + (1.0/3.0) * -20.503 * tfactors.T923i
                                      + 1.39981 + (5.0/3.0) * -0.101641 * tfactors.T923 + 10.082 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ni60_to_n_Ni59_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni60_to_n_Ni59_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co57_n_to_Co58_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Co57 + n --> Co58

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  16.9961 + -1.03654 * tfactors.T913
                             + 0.342218 * tfactors.T9 + -0.0638367 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -1.03654 * tfactors.T923i
                                      + 0.342218 + (5.0/3.0) * -0.0638367 * tfactors.T923;

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

        rate_eval.screened_rates(k_Co57_n_to_Co58_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co57_n_to_Co58_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co58_to_n_Co57_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Co58 --> n + Co57

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  41.1458 + -99.4851 * tfactors.T9i + -1.03654 * tfactors.T913
                             + 0.342218 * tfactors.T9 + -0.0638367 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  99.4851 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.03654 * tfactors.T923i
                                      + 0.342218 + (5.0/3.0) * -0.0638367 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Co58_to_n_Co57_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co58_to_n_Co57_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co58_n_to_Co59_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Co58 + n --> Co59

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  18.1906 + -2.18418 * tfactors.T913
                             + 0.328629 * tfactors.T9 + -0.0407234 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -2.18418 * tfactors.T923i
                                      + 0.328629 + (5.0/3.0) * -0.0407234 * tfactors.T923;

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

        rate_eval.screened_rates(k_Co58_n_to_Co59_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co58_n_to_Co59_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co59_to_n_Co58_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Co59 --> n + Co58

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  41.4007 + -121.312 * tfactors.T9i + -2.18418 * tfactors.T913
                             + 0.328629 * tfactors.T9 + -0.0407234 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  121.312 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -2.18418 * tfactors.T923i
                                      + 0.328629 + (5.0/3.0) * -0.0407234 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Co59_to_n_Co58_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co59_to_n_Co58_reaclib) = drate_dT;
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
    void rate_Fe53_to_n_Fe52_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Fe53 --> n + Fe52

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  37.4863 + -123.994 * tfactors.T9i + -0.344319 * tfactors.T913
                             + 0.178277 * tfactors.T9 + -0.0334326 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  123.994 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.344319 * tfactors.T923i
                                      + 0.178277 + (5.0/3.0) * -0.0334326 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Fe53_to_n_Fe52_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe53_to_n_Fe52_reaclib) = drate_dT;
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
    void rate_Fe54_to_n_Fe53_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Fe54 --> n + Fe53

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  42.2106 + -155.251 * tfactors.T9i + -1.10421 * tfactors.T913
                             + 0.379905 * tfactors.T9 + -0.0581878 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  155.251 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.10421 * tfactors.T923i
                                      + 0.379905 + (5.0/3.0) * -0.0581878 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Fe54_to_n_Fe53_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe54_to_n_Fe53_reaclib) = drate_dT;
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
    void rate_Fe55_to_n_Fe54_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Fe55 --> n + Fe54

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  21.4833 + -107.843 * tfactors.T9i + -8.66617 * tfactors.T913i + 26.4472 * tfactors.T913
                             + -1.9222 * tfactors.T9 + 0.0986404 * tfactors.T953 + -8.28317 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  107.843 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -8.66617 * tfactors.T943i + (1.0/3.0) * 26.4472 * tfactors.T923i
                                      + -1.9222 + (5.0/3.0) * 0.0986404 * tfactors.T923 + -8.28317 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Fe55_to_n_Fe54_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe55_to_n_Fe54_reaclib) = drate_dT;
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
    void rate_Fe56_to_n_Fe55_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Fe56 --> n + Fe55

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  46.7852 + -130.023 * tfactors.T9i + 8.06062 * tfactors.T913i + -14.4809 * tfactors.T913
                             + 0.94252 * tfactors.T9 + -0.0776007 * tfactors.T953 + 7.97093 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  130.023 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 8.06062 * tfactors.T943i + (1.0/3.0) * -14.4809 * tfactors.T923i
                                      + 0.94252 + (5.0/3.0) * -0.0776007 * tfactors.T923 + 7.97093 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Fe56_to_n_Fe55_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe56_to_n_Fe55_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe56_n_to_Fe57_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Fe56 + n --> Fe57

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  35.834 + -0.0600873 * tfactors.T9i + 9.79414 * tfactors.T913i + -33.7366 * tfactors.T913
                             + 2.97173 * tfactors.T9 + -0.236498 * tfactors.T953 + 11.8711 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  0.0600873 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 9.79414 * tfactors.T943i + (1.0/3.0) * -33.7366 * tfactors.T923i
                                      + 2.97173 + (5.0/3.0) * -0.236498 * tfactors.T923 + 11.8711 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Fe56_n_to_Fe57_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe56_n_to_Fe57_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe57_to_n_Fe56_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Fe57 --> n + Fe56

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  58.8201 + -88.7916 * tfactors.T9i + 9.79414 * tfactors.T913i + -33.7366 * tfactors.T913
                             + 2.97173 * tfactors.T9 + -0.236498 * tfactors.T953 + 13.3711 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  88.7916 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 9.79414 * tfactors.T943i + (1.0/3.0) * -33.7366 * tfactors.T923i
                                      + 2.97173 + (5.0/3.0) * -0.236498 * tfactors.T923 + 13.3711 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Fe57_to_n_Fe56_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe57_to_n_Fe56_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe57_n_to_Fe58_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Fe57 + n --> Fe58

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  20.2493 + -0.0214497 * tfactors.T9i + 2.1643 * tfactors.T913i + -7.49584 * tfactors.T913
                             + 0.656354 * tfactors.T9 + -0.0549855 * tfactors.T953 + 2.38149 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  0.0214497 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 2.1643 * tfactors.T943i + (1.0/3.0) * -7.49584 * tfactors.T923i
                                      + 0.656354 + (5.0/3.0) * -0.0549855 * tfactors.T923 + 2.38149 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Fe57_n_to_Fe58_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe57_n_to_Fe58_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe58_to_n_Fe57_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Fe58 --> n + Fe57

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  44.6221 + -116.581 * tfactors.T9i + 2.1643 * tfactors.T913i + -7.49584 * tfactors.T913
                             + 0.656354 * tfactors.T9 + -0.0549855 * tfactors.T953 + 3.88149 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  116.581 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 2.1643 * tfactors.T943i + (1.0/3.0) * -7.49584 * tfactors.T923i
                                      + 0.656354 + (5.0/3.0) * -0.0549855 * tfactors.T923 + 3.88149 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Fe58_to_n_Fe57_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe58_to_n_Fe57_reaclib) = drate_dT;
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
    void rate_Mn52_to_n_Mn51_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Mn52 --> n + Mn51

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  40.109 + -122.258 * tfactors.T9i + -1.83611 * tfactors.T913
                             + 0.460384 * tfactors.T9 + -0.0584947 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  122.258 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.83611 * tfactors.T923i
                                      + 0.460384 + (5.0/3.0) * -0.0584947 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Mn52_to_n_Mn51_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn52_to_n_Mn51_reaclib) = drate_dT;
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
    void rate_Mn53_to_n_Mn52_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Mn53 --> n + Mn52

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  41.2568 + -139.878 * tfactors.T9i + -0.91084 * tfactors.T913
                             + 0.162511 * tfactors.T9 + -0.0279619 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  139.878 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.91084 * tfactors.T923i
                                      + 0.162511 + (5.0/3.0) * -0.0279619 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Mn53_to_n_Mn52_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn53_to_n_Mn52_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn53_n_to_Mn54_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Mn53 + n --> Mn54

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  17.3501 + -1.70171 * tfactors.T913
                             + 0.386422 * tfactors.T9 + -0.0536858 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -1.70171 * tfactors.T923i
                                      + 0.386422 + (5.0/3.0) * -0.0536858 * tfactors.T923;

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

        rate_eval.screened_rates(k_Mn53_n_to_Mn54_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn53_n_to_Mn54_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn54_to_n_Mn53_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Mn54 --> n + Mn53

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  41.1614 + -103.73 * tfactors.T9i + -1.70171 * tfactors.T913
                             + 0.386422 * tfactors.T9 + -0.0536858 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  103.73 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.70171 * tfactors.T923i
                                      + 0.386422 + (5.0/3.0) * -0.0536858 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Mn54_to_n_Mn53_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn54_to_n_Mn53_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn54_n_to_Mn55_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Mn54 + n --> Mn55

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  18.4705 + -2.56572 * tfactors.T913
                             + 0.352447 * tfactors.T9 + -0.0375878 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -2.56572 * tfactors.T923i
                                      + 0.352447 + (5.0/3.0) * -0.0375878 * tfactors.T923;

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

        rate_eval.screened_rates(k_Mn54_n_to_Mn55_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn54_n_to_Mn55_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn55_to_n_Mn54_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Mn55 --> n + Mn54

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  42.3029 + -118.673 * tfactors.T9i + -2.56572 * tfactors.T913
                             + 0.352447 * tfactors.T9 + -0.0375878 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  118.673 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -2.56572 * tfactors.T923i
                                      + 0.352447 + (5.0/3.0) * -0.0375878 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Mn55_to_n_Mn54_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn55_to_n_Mn54_reaclib) = drate_dT;
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
    void rate_Cr49_to_n_Cr48_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cr49 --> n + Cr48

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  37.7375 + -122.806 * tfactors.T9i + -0.000452357 * tfactors.T913
                             + 0.0412838 * tfactors.T9 + -0.0160865 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  122.806 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.000452357 * tfactors.T923i
                                      + 0.0412838 + (5.0/3.0) * -0.0160865 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cr49_to_n_Cr48_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr49_to_n_Cr48_reaclib) = drate_dT;
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
    void rate_Cr50_to_n_Cr49_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cr50 --> n + Cr49

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  43.3027 + -150.862 * tfactors.T9i + -1.98927 * tfactors.T913
                             + 0.465124 * tfactors.T9 + -0.0582249 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  150.862 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.98927 * tfactors.T923i
                                      + 0.465124 + (5.0/3.0) * -0.0582249 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cr50_to_n_Cr49_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr50_to_n_Cr49_reaclib) = drate_dT;
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
    void rate_Cr51_to_n_Cr50_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cr51 --> n + Cr50

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  36.7032 + -107.484 * tfactors.T9i + 0.815037 * tfactors.T913i + -0.708732 * tfactors.T913
                             + 0.386295 * tfactors.T9 + -0.0615661 * tfactors.T953 + 1.54546 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  107.484 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 0.815037 * tfactors.T943i + (1.0/3.0) * -0.708732 * tfactors.T923i
                                      + 0.386295 + (5.0/3.0) * -0.0615661 * tfactors.T923 + 1.54546 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cr51_to_n_Cr50_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr51_to_n_Cr50_reaclib) = drate_dT;
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
    void rate_Cr52_to_n_Cr51_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cr52 --> n + Cr51

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  45.8814 + -139.784 * tfactors.T9i + 6.95563 * tfactors.T913i + -11.4485 * tfactors.T913
                             + 0.736956 * tfactors.T9 + -0.0695323 * tfactors.T953 + 6.81244 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  139.784 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 6.95563 * tfactors.T943i + (1.0/3.0) * -11.4485 * tfactors.T923i
                                      + 0.736956 + (5.0/3.0) * -0.0695323 * tfactors.T923 + 6.81244 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cr52_to_n_Cr51_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr52_to_n_Cr51_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr52_n_to_Cr53_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cr52 + n --> Cr53

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  27.3183 + -0.0314027 * tfactors.T9i + 5.71728 * tfactors.T913i + -20.4149 * tfactors.T913
                             + 1.95304 * tfactors.T9 + -0.175538 * tfactors.T953 + 7.18293 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  0.0314027 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 5.71728 * tfactors.T943i + (1.0/3.0) * -20.4149 * tfactors.T923i
                                      + 1.95304 + (5.0/3.0) * -0.175538 * tfactors.T923 + 7.18293 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cr52_n_to_Cr53_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr52_n_to_Cr53_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr53_to_n_Cr52_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cr53 --> n + Cr52

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  49.6092 + -92.163 * tfactors.T9i + 5.71728 * tfactors.T913i + -20.4149 * tfactors.T913
                             + 1.95304 * tfactors.T9 + -0.175538 * tfactors.T953 + 8.68293 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  92.163 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 5.71728 * tfactors.T943i + (1.0/3.0) * -20.4149 * tfactors.T923i
                                      + 1.95304 + (5.0/3.0) * -0.175538 * tfactors.T923 + 8.68293 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cr53_to_n_Cr52_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr53_to_n_Cr52_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr53_n_to_Cr54_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cr53 + n --> Cr54

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  10.8793 + -0.00984544 * tfactors.T9i + -1.26805 * tfactors.T913i + 6.34126 * tfactors.T913
                             + -0.39287 * tfactors.T9 + 0.000487041 * tfactors.T953 + -2.46456 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  0.00984544 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -1.26805 * tfactors.T943i + (1.0/3.0) * 6.34126 * tfactors.T923i
                                      + -0.39287 + (5.0/3.0) * 0.000487041 * tfactors.T923 + -2.46456 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cr53_n_to_Cr54_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr53_n_to_Cr54_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr54_to_n_Cr53_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Cr54 --> n + Cr53

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  35.9433 + -112.797 * tfactors.T9i + -1.26805 * tfactors.T913i + 6.34126 * tfactors.T913
                             + -0.39287 * tfactors.T9 + 0.000487041 * tfactors.T953 + -0.96456 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  112.797 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -1.26805 * tfactors.T943i + (1.0/3.0) * 6.34126 * tfactors.T923i
                                      + -0.39287 + (5.0/3.0) * 0.000487041 * tfactors.T923 + -0.96456 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Cr54_to_n_Cr53_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr54_to_n_Cr53_reaclib) = drate_dT;
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
    void rate_V48_to_n_V47_reaclib(T& rate_eval, const tf_t& tfactors) {

        // V48 --> n + V47

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  40.4146 + -122.366 * tfactors.T9i + -2.40159 * tfactors.T913
                             + 0.594573 * tfactors.T9 + -0.0682896 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  122.366 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -2.40159 * tfactors.T923i
                                      + 0.594573 + (5.0/3.0) * -0.0682896 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_V48_to_n_V47_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V48_to_n_V47_reaclib) = drate_dT;
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
    void rate_V49_to_n_V48_reaclib(T& rate_eval, const tf_t& tfactors) {

        // V49 --> n + V48

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  41.6761 + -134.066 * tfactors.T9i + -0.917026 * tfactors.T913
                             + -0.109162 * tfactors.T9 + 0.0127488 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  134.066 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.917026 * tfactors.T923i
                                      + -0.109162 + (5.0/3.0) * 0.0127488 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_V49_to_n_V48_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V49_to_n_V48_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V49_n_to_V50_reaclib(T& rate_eval, const tf_t& tfactors) {

        // V49 + n --> V50

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  17.5445 + -1.8725 * tfactors.T913
                             + 0.291609 * tfactors.T9 + -0.0370216 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -1.8725 * tfactors.T923i
                                      + 0.291609 + (5.0/3.0) * -0.0370216 * tfactors.T923;

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

        rate_eval.screened_rates(k_V49_n_to_V50_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V49_n_to_V50_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V50_to_n_V49_reaclib(T& rate_eval, const tf_t& tfactors) {

        // V50 --> n + V49

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  40.7345 + -108.339 * tfactors.T9i + -1.8725 * tfactors.T913
                             + 0.291609 * tfactors.T9 + -0.0370216 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  108.339 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.8725 * tfactors.T923i
                                      + 0.291609 + (5.0/3.0) * -0.0370216 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_V50_to_n_V49_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V50_to_n_V49_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V50_n_to_V51_reaclib(T& rate_eval, const tf_t& tfactors) {

        // V50 + n --> V51

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  26.4771 + -0.118554 * tfactors.T9i + 10.9535 * tfactors.T913i + -23.1254 * tfactors.T913
                             + 1.53091 * tfactors.T9 + -0.0992995 * tfactors.T953 + 9.62203 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  0.118554 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 10.9535 * tfactors.T943i + (1.0/3.0) * -23.1254 * tfactors.T923i
                                      + 1.53091 + (5.0/3.0) * -0.0992995 * tfactors.T923 + 9.62203 * tfactors.T9i;

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

        rate_eval.screened_rates(k_V50_n_to_V51_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V50_n_to_V51_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V51_to_n_V50_reaclib(T& rate_eval, const tf_t& tfactors) {

        // V51 --> n + V50

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  50.6387 + -128.363 * tfactors.T9i + 10.9535 * tfactors.T913i + -23.1254 * tfactors.T913
                             + 1.53091 * tfactors.T9 + -0.0992995 * tfactors.T953 + 11.122 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  128.363 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 10.9535 * tfactors.T943i + (1.0/3.0) * -23.1254 * tfactors.T923i
                                      + 1.53091 + (5.0/3.0) * -0.0992995 * tfactors.T923 + 11.122 * tfactors.T9i;

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

        rate_eval.screened_rates(k_V51_to_n_V50_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V51_to_n_V50_reaclib) = drate_dT;
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
    void rate_Ti45_to_n_Ti44_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ti45 --> n + Ti44

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  36.8432 + -110.575 * tfactors.T9i + 1.01203 * tfactors.T913
                             + -0.201174 * tfactors.T9 + 0.00360954 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  110.575 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 1.01203 * tfactors.T923i
                                      + -0.201174 + (5.0/3.0) * 0.00360954 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ti45_to_n_Ti44_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti45_to_n_Ti44_reaclib) = drate_dT;
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
    void rate_Ti46_to_n_Ti45_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ti46 --> n + Ti45

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  42.8852 + -153.052 * tfactors.T9i + -1.55633 * tfactors.T913
                             + 0.300783 * tfactors.T9 + -0.0373696 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  153.052 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.55633 * tfactors.T923i
                                      + 0.300783 + (5.0/3.0) * -0.0373696 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ti46_to_n_Ti45_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti46_to_n_Ti45_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti46_n_to_Ti47_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ti46 + n --> Ti47

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  12.2761 + 0.0882347 * tfactors.T9i + -7.30817 * tfactors.T913i + 10.713 * tfactors.T913
                             + -0.537251 * tfactors.T9 + 0.0187304 * tfactors.T953 + -5.38893 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  -0.0882347 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -7.30817 * tfactors.T943i + (1.0/3.0) * 10.713 * tfactors.T923i
                                      + -0.537251 + (5.0/3.0) * 0.0187304 * tfactors.T923 + -5.38893 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ti46_n_to_Ti47_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti46_n_to_Ti47_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti47_to_n_Ti46_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ti47 --> n + Ti46

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  34.1579 + -102.963 * tfactors.T9i + -7.30817 * tfactors.T913i + 10.713 * tfactors.T913
                             + -0.537251 * tfactors.T9 + 0.0187304 * tfactors.T953 + -3.88893 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  102.963 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -7.30817 * tfactors.T943i + (1.0/3.0) * 10.713 * tfactors.T923i
                                      + -0.537251 + (5.0/3.0) * 0.0187304 * tfactors.T923 + -3.88893 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ti47_to_n_Ti46_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti47_to_n_Ti46_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti47_n_to_Ti48_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ti47 + n --> Ti48

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  14.8734 + 0.0181751 * tfactors.T9i + -2.34945 * tfactors.T913i + 3.2238 * tfactors.T913
                             + -0.0455487 * tfactors.T9 + -0.0153771 * tfactors.T953 + -2.10785 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  -0.0181751 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -2.34945 * tfactors.T943i + (1.0/3.0) * 3.2238 * tfactors.T923i
                                      + -0.0455487 + (5.0/3.0) * -0.0153771 * tfactors.T923 + -2.10785 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ti47_n_to_Ti48_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti47_n_to_Ti48_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti48_to_n_Ti47_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ti48 --> n + Ti47

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  40.3394 + -134.911 * tfactors.T9i + -2.34945 * tfactors.T913i + 3.2238 * tfactors.T913
                             + -0.0455487 * tfactors.T9 + -0.0153771 * tfactors.T953 + -0.60785 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  134.911 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -2.34945 * tfactors.T943i + (1.0/3.0) * 3.2238 * tfactors.T923i
                                      + -0.0455487 + (5.0/3.0) * -0.0153771 * tfactors.T923 + -0.60785 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ti48_to_n_Ti47_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti48_to_n_Ti47_reaclib) = drate_dT;
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
    void rate_Sc42_to_n_Sc41_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Sc42 --> n + Sc41

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  36.3505 + -134.034 * tfactors.T9i + -1.01509 * tfactors.T913
                             + 0.362316 * tfactors.T9 + -0.0412181 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  134.034 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.01509 * tfactors.T923i
                                      + 0.362316 + (5.0/3.0) * -0.0412181 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Sc42_to_n_Sc41_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc42_to_n_Sc41_reaclib) = drate_dT;
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
    void rate_Sc43_to_n_Sc42_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Sc43 --> n + Sc42

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  36.5353 + -140.855 * tfactors.T9i + -0.921604 * tfactors.T913
                             + 0.333782 * tfactors.T9 + -0.0404241 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  140.855 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.921604 * tfactors.T923i
                                      + 0.333782 + (5.0/3.0) * -0.0404241 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Sc43_to_n_Sc42_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc43_to_n_Sc42_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc43_n_to_Sc44_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Sc43 + n --> Sc44

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  15.7633 + -0.575865 * tfactors.T913
                             + 0.0565199 * tfactors.T9 + -0.0129886 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -0.575865 * tfactors.T923i
                                      + 0.0565199 + (5.0/3.0) * -0.0129886 * tfactors.T923;

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

        rate_eval.screened_rates(k_Sc43_n_to_Sc44_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc43_n_to_Sc44_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc44_to_n_Sc43_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Sc44 --> n + Sc43

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  39.9046 + -112.558 * tfactors.T9i + -0.575865 * tfactors.T913
                             + 0.0565199 * tfactors.T9 + -0.0129886 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  112.558 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.575865 * tfactors.T923i
                                      + 0.0565199 + (5.0/3.0) * -0.0129886 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Sc44_to_n_Sc43_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc44_to_n_Sc43_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc44_n_to_Sc45_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Sc44 + n --> Sc45

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  17.0705 + -1.31922 * tfactors.T913
                             + 0.167096 * tfactors.T9 + -0.0191676 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  + (1.0/3.0) * -1.31922 * tfactors.T923i
                                      + 0.167096 + (5.0/3.0) * -0.0191676 * tfactors.T923;

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

        rate_eval.screened_rates(k_Sc44_n_to_Sc45_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc44_n_to_Sc45_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc45_to_n_Sc44_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Sc45 --> n + Sc44

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ths8r
        ln_set_rate =  40.2726 + -131.398 * tfactors.T9i + -1.31922 * tfactors.T913
                             + 0.167096 * tfactors.T9 + -0.0191676 * tfactors.T953 + 1.5 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  131.398 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -1.31922 * tfactors.T923i
                                      + 0.167096 + (5.0/3.0) * -0.0191676 * tfactors.T923 + 1.5 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Sc45_to_n_Sc44_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc45_to_n_Sc44_reaclib) = drate_dT;
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
    void rate_Ca41_to_n_Ca40_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ca41 --> n + Ca40

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  24.6353 + -96.9664 * tfactors.T9i + -9.18636 * tfactors.T913i + 21.0095 * tfactors.T913
                             + -1.33326 * tfactors.T9 + 0.0769347 * tfactors.T953 + -7.16245 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  96.9664 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -9.18636 * tfactors.T943i + (1.0/3.0) * 21.0095 * tfactors.T923i
                                      + -1.33326 + (5.0/3.0) * 0.0769347 * tfactors.T923 + -7.16245 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ca41_to_n_Ca40_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca41_to_n_Ca40_reaclib) = drate_dT;
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
    void rate_Ca42_to_n_Ca41_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Ca42 --> n + Ca41

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // ks03 
        ln_set_rate =  46.1952 + -133.264 * tfactors.T9i + 3.7037 * tfactors.T913i + -9.99246 * tfactors.T913
                             + 1.05894 * tfactors.T9 + -0.0978564 * tfactors.T953 + 4.93604 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  133.264 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * 3.7037 * tfactors.T943i + (1.0/3.0) * -9.99246 * tfactors.T923i
                                      + 1.05894 + (5.0/3.0) * -0.0978564 * tfactors.T923 + 4.93604 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Ca42_to_n_Ca41_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca42_to_n_Ca41_reaclib) = drate_dT;
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
    void rate_Na23_p_to_C12_C12_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Na23 + p --> C12 + C12

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // cf88r
        ln_set_rate =  60.9438 + -26.0184 * tfactors.T9i + -84.165 * tfactors.T913i + -1.4191 * tfactors.T913
                             + -0.114619 * tfactors.T9 + -0.070307 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  26.0184 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -84.165 * tfactors.T943i + (1.0/3.0) * -1.4191 * tfactors.T923i
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

        rate_eval.screened_rates(k_Na23_p_to_C12_C12_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Na23_p_to_C12_C12_reaclib) = drate_dT;
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
    void rate_Mg24_to_p_Na23_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Mg24 --> p + Na23

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10r
        ln_set_rate =  34.0876 + -138.968 * tfactors.T9i + -0.360588 * tfactors.T913
                             + 1.4187 * tfactors.T9 + -0.184061 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  138.968 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * -0.360588 * tfactors.T923i
                                      + 1.4187 + (5.0/3.0) * -0.184061 * tfactors.T923;

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
        ln_set_rate =  20.0024 + -137.3 * tfactors.T9i;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  137.3 * tfactors.T9i * tfactors.T9i;

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
        ln_set_rate =  43.9357 + -135.688 * tfactors.T9i + -20.6428 * tfactors.T913i + 1.52954 * tfactors.T913
                             + 2.7487 * tfactors.T9 + -1.0 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  135.688 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -20.6428 * tfactors.T943i + (1.0/3.0) * 1.52954 * tfactors.T923i
                                      + 2.7487 + (5.0/3.0) * -1.0 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Mg24_to_p_Na23_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg24_to_p_Na23_reaclib) = drate_dT;
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
    void rate_Ne20_He4_to_C12_C12_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ne20 + He4 --> C12 + C12

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // cf88r
        ln_set_rate =  61.4748 + -53.6267 * tfactors.T9i + -84.165 * tfactors.T913i + -1.56627 * tfactors.T913
                             + -0.0736084 * tfactors.T9 + -0.072797 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  53.6267 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -84.165 * tfactors.T943i + (1.0/3.0) * -1.56627 * tfactors.T923i
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

        rate_eval.screened_rates(k_Ne20_He4_to_C12_C12_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne20_He4_to_C12_C12_reaclib) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ne20_He4_to_p_Na23_reaclib(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ne20 + He4 --> p + Na23

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10r
        ln_set_rate =  -6.37772 + -29.8896 * tfactors.T9i + 19.7297 * tfactors.T913
                             + -2.20987 * tfactors.T9 + 0.153374 * tfactors.T953 + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  29.8896 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 19.7297 * tfactors.T923i
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
        ln_set_rate =  0.227472 + -29.4348 * tfactors.T9i
                             + -1.5 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  29.4348 * tfactors.T9i * tfactors.T9i
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
        ln_set_rate =  19.1852 + -27.5738 * tfactors.T9i + -20.0024 * tfactors.T913i + 11.5988 * tfactors.T913
                             + -1.37398 * tfactors.T9 + -1.0 * tfactors.T953 + -0.666667 * tfactors.lnT9;

        ln_set_rate += log_scor;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  27.5738 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -20.0024 * tfactors.T943i + (1.0/3.0) * 11.5988 * tfactors.T923i
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

        rate_eval.screened_rates(k_Ne20_He4_to_p_Na23_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne20_He4_to_p_Na23_reaclib) = drate_dT;
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
    void rate_Mg24_to_He4_Ne20_reaclib(T& rate_eval, const tf_t& tfactors) {

        // Mg24 --> He4 + Ne20

        amrex::Real rate{}, drate_dT{};
        amrex::Real ln_set_rate{0.0};
        amrex::Real dln_set_rate_dT9{0.0};
        amrex::Real set_rate{0.0};

        // il10r
        ln_set_rate =  16.0203 + -120.895 * tfactors.T9i + 16.9229 * tfactors.T913
                             + -2.57325 * tfactors.T9 + 0.208997 * tfactors.T953;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  120.895 * tfactors.T9i * tfactors.T9i + (1.0/3.0) * 16.9229 * tfactors.T923i
                                      + -2.57325 + (5.0/3.0) * 0.208997 * tfactors.T923;

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
        ln_set_rate =  26.8017 + -117.334 * tfactors.T9i;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  117.334 * tfactors.T9i * tfactors.T9i;

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
        ln_set_rate =  -13.8869 + -110.62 * tfactors.T9i;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  110.62 * tfactors.T9i * tfactors.T9i;

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
        ln_set_rate =  49.3244 + -108.114 * tfactors.T9i + -46.2525 * tfactors.T913i + 5.58901 * tfactors.T913
                             + 7.61843 * tfactors.T9 + -3.683 * tfactors.T953 + 0.833333 * tfactors.lnT9;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            dln_set_rate_dT9 =  108.114 * tfactors.T9i * tfactors.T9i + -(1.0/3.0) * -46.2525 * tfactors.T943i + (1.0/3.0) * 5.58901 * tfactors.T923i
                                      + 7.61843 + (5.0/3.0) * -3.683 * tfactors.T923 + 0.833333 * tfactors.T9i;

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

        rate_eval.screened_rates(k_Mg24_to_He4_Ne20_reaclib) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg24_to_He4_Ne20_reaclib) = drate_dT;
        }
    }


}
#endif
