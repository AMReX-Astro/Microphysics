#ifndef MODIFIED_RATES_H
#define MODIFIED_RATES_H

#include <AMReX.H>
#include <AMReX_Print.H>

#include <tfactors.H>
#include <actual_network.H>
#include <rate_type.H>
#include <temperature_table_rates.H>
#include <reaclib_rates.H>

namespace modified_rates {

    using namespace Rates;

    using namespace reaclib_rates;
    using namespace temp_tabular;

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_C12_C12_to_Mg24_modified(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // C12 + C12 --> Mg24 (calls the underlying rate)
        // represents the sequence: c12(c12,n)mg23(n,γ)mg24

        rate_C12_C12_to_n_Mg23_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        rate_eval.screened_rates(k_C12_C12_to_Mg24_modified) = rate_eval.screened_rates(k_C12_C12_to_n_Mg23_reaclib);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_C12_C12_to_Mg24_modified) = rate_eval.dscreened_rates_dT(k_C12_C12_to_n_Mg23_reaclib);
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_O16_O16_to_S32_modified(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // O16 + O16 --> S32 (calls the underlying rate)
        // represents the sequence: o16(o16,n)s31(n,γ)s32

        rate_O16_O16_to_n_S31_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        rate_eval.screened_rates(k_O16_O16_to_S32_modified) = rate_eval.screened_rates(k_O16_O16_to_n_S31_reaclib);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_O16_O16_to_S32_modified) = rate_eval.dscreened_rates_dT(k_O16_O16_to_n_S31_reaclib);
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_O16_C12_to_Si28_modified(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // O16 + C12 --> Si28 (calls the underlying rate)
        // represents the sequence: o16(c12,n)si27(n,γ)si28

        rate_O16_C12_to_n_Si27_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        rate_eval.screened_rates(k_O16_C12_to_Si28_modified) = rate_eval.screened_rates(k_O16_C12_to_n_Si27_reaclib);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_O16_C12_to_Si28_modified) = rate_eval.dscreened_rates_dT(k_O16_C12_to_n_Si27_reaclib);
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Na22_1p5He4_to_Si28_modified(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Na22 + 1.5 He4 --> Si28 (calls the underlying rate)
        rate_Na22_He4_to_Al26_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        rate_eval.screened_rates(k_Na22_1p5He4_to_Si28_modified) = rate_eval.screened_rates(k_Na22_He4_to_Al26_reaclib);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Na22_1p5He4_to_Si28_modified) = rate_eval.dscreened_rates_dT(k_Na22_He4_to_Al26_reaclib);
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mg24_p_p_to_Na22_He4_beta_pos_modified(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Mg24 + p + p --> Na22 + He4 (calls the underlying rate)
        rate_Mg24_p_to_Al25_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        rate_eval.screened_rates(k_Mg24_p_p_to_Na22_He4_beta_pos_modified) = rate_eval.screened_rates(k_Mg24_p_to_Al25_reaclib);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg24_p_p_to_Na22_He4_beta_pos_modified) = rate_eval.dscreened_rates_dT(k_Mg24_p_to_Al25_reaclib);
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
            // C12_C12_to_Mg24_modified
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_C12_C12);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_C12);
            }
#endif
            rate_C12_C12_to_Mg24_modified(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // O16_O16_to_S32_modified
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_O16_O16);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_O16_O16);
            }
#endif
            rate_O16_O16_to_S32_modified(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // O16_C12_to_Si28_modified
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_C12_O16);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_O16);
            }
#endif
            rate_O16_C12_to_Si28_modified(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Na22_1p5He4_to_Si28_modified
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_Na22);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_Na22);
            }
#endif
            rate_Na22_1p5He4_to_Si28_modified(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }

        {
            // Mg24_p_p_to_Na22_He4_beta_pos_modified
            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_Mg24);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_Mg24);
            }
#endif
            rate_Mg24_p_p_to_Na22_He4_beta_pos_modified(rate_eval, tfactors, log_scor, dlog_scor_dT);
        }


    }
}
#endif
