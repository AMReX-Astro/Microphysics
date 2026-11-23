#ifndef MODIFIED_RATES_H
#define MODIFIED_RATES_H

#include <AMReX.H>
#include <AMReX_Print.H>

#include <tfactors.H>
#include <actual_network.H>
#include <rate_type.H>
#include <temperature_table_rates.H>
#include <reaclib_rates.H>
#include <beta_limited_rates.H>

namespace modified_rates {

    using namespace Rates;

    using namespace reaclib_rates;
    using namespace temp_tabular;
    using namespace beta_limited_rates;

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_C12_C12_to_Mg24_modified(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // C12 + C12 --> Mg24 (calls the underlying rate)
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
        rate_O16_C12_to_n_Si27_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        rate_eval.screened_rates(k_O16_C12_to_Si28_modified) = rate_eval.screened_rates(k_O16_C12_to_n_Si27_reaclib);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_O16_C12_to_Si28_modified) = rate_eval.dscreened_rates_dT(k_O16_C12_to_n_Si27_reaclib);
        }
    }


}
#endif
