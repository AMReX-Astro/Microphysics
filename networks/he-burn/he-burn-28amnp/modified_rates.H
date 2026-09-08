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

    using namespace temp_tabular;

    template <int do_T_derivatives>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_He4_N14_to_Ne20_modified(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

        // N14 + 1.5 He4 --> Ne20 (calls the underlying rate)

        rate_He4_N14_to_F18_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
    }

    template <int do_T_derivatives>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_p_O16_to_N14_He4_modified(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

        // O16 + p + p --> N14 + He4 (calls the underlying rate)

        rate_p_O16_to_F17_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
    }

    template <int do_T_derivatives>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_C12_C12_to_Mg24_modified(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

        // C12 + C12 --> Mg24 (calls the underlying rate)

        rate_C12_C12_to_n_Mg23_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
    }

    template <int do_T_derivatives>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_O16_O16_to_S32_modified(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

        // O16 + O16 --> S32 (calls the underlying rate)

        rate_O16_O16_to_n_S31_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
    }

    template <int do_T_derivatives>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_C12_O16_to_Si28_modified(const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT, amrex::Real& rate, amrex::Real& drate_dT) {

        // O16 + C12 --> Si28 (calls the underlying rate)

        rate_C12_O16_to_n_Si27_reaclib<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
    }


    template <int do_T_derivatives, typename T>
    AMREX_GPU_HOST_DEVICE AMREX_INLINE
    void
    fill_rates([[maybe_unused]] const tf_t& tfactors,
               [[maybe_unused]] T& rate_eval)
    {

        [[maybe_unused]] amrex::Real rate;
        [[maybe_unused]] amrex::Real drate_dT;

        {
            // He4_N14_to_Ne20_modified

            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_He4_N14);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_He4_N14);
            }
#endif
            rate_He4_N14_to_Ne20_modified<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
            rate_eval.screened_rates(k_He4_N14_to_Ne20_modified) = rate;
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                rate_eval.dscreened_rates_dT(k_He4_N14_to_Ne20_modified) = drate_dT;
            }
        }

        {
            // p_O16_to_N14_He4_modified

            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_p_O16);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_p_O16);
            }
#endif
            rate_p_O16_to_N14_He4_modified<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
            rate_eval.screened_rates(k_p_O16_to_N14_He4_modified) = rate;
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                rate_eval.dscreened_rates_dT(k_p_O16_to_N14_He4_modified) = drate_dT;
            }
        }

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
            rate_C12_C12_to_Mg24_modified<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
            rate_eval.screened_rates(k_C12_C12_to_Mg24_modified) = rate;
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                rate_eval.dscreened_rates_dT(k_C12_C12_to_Mg24_modified) = drate_dT;
            }
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
            rate_O16_O16_to_S32_modified<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
            rate_eval.screened_rates(k_O16_O16_to_S32_modified) = rate;
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                rate_eval.dscreened_rates_dT(k_O16_O16_to_S32_modified) = drate_dT;
            }
        }

        {
            // C12_O16_to_Si28_modified

            amrex::Real log_scor {0.0_rt};
            amrex::Real dlog_scor_dT {0.0_rt};
#ifdef SCREENING
            log_scor = rate_eval.log_screen(k_C12_O16);
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                dlog_scor_dT = rate_eval.dlog_screen_dT(k_C12_O16);
            }
#endif
            rate_C12_O16_to_Si28_modified<do_T_derivatives>(tfactors, log_scor, dlog_scor_dT, rate, drate_dT);
            rate_eval.screened_rates(k_C12_O16_to_Si28_modified) = rate;
            if constexpr (std::is_same_v<T, rate_derivs_t>) {
                rate_eval.dscreened_rates_dT(k_C12_O16_to_Si28_modified) = drate_dT;
            }
        }


    }
}
#endif
