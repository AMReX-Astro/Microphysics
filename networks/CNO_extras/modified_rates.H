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
    void rate_Ne19_3p_to_Mg22_beta_pos_modified(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ne19 + 3 p --> Mg22 (calls the underlying rate)
        // represents the sequence: Ne19(p,γ)Na20(p,γ)Mg21(β+)Na21(p,γ)Mg22

        rate_Ne19_p_to_Na20_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        rate_eval.screened_rates(k_Ne19_3p_to_Mg22_beta_pos_modified) = rate_eval.screened_rates(k_Ne19_p_to_Na20_reaclib);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne19_3p_to_Mg22_beta_pos_modified) = rate_eval.dscreened_rates_dT(k_Ne19_p_to_Na20_reaclib);
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ne20_p_p_to_Mg22_modified(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // Ne20 + p + p --> Mg22 (calls the underlying rate)
        // represents the sequence: Ne20(p,γ)Na21(p,γ)Mg22

        rate_Ne20_p_to_Na21_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        rate_eval.screened_rates(k_Ne20_p_p_to_Mg22_modified) = rate_eval.screened_rates(k_Ne20_p_to_Na21_reaclib);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne20_p_p_to_Mg22_modified) = rate_eval.dscreened_rates_dT(k_Ne20_p_to_Na21_reaclib);
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_p_p_p_to_He3_beta_pos_modified(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // 3 p --> He3 (calls the underlying rate)
        // represents the sequence: p(p,e⁺ν)d(p,γ)He3

        rate_p_p_to_d_beta_pos_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        rate_eval.screened_rates(k_p_p_p_to_He3_beta_pos_modified) = rate_eval.screened_rates(k_p_p_to_d_beta_pos_reaclib);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_p_p_to_He3_beta_pos_modified) = rate_eval.dscreened_rates_dT(k_p_p_to_d_beta_pos_reaclib);
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_p_p_p_to_He3_electron_capture_modified(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // 3 p --> He3 (calls the underlying rate)
        // represents the sequence: p(pe⁻,ν)d(p,γ)He3

        rate_p_p_to_d_electron_capture_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        rate_eval.screened_rates(k_p_p_p_to_He3_electron_capture_modified) = rate_eval.screened_rates(k_p_p_to_d_electron_capture_reaclib);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_p_p_to_He3_electron_capture_modified) = rate_eval.dscreened_rates_dT(k_p_p_to_d_electron_capture_reaclib);
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_He3_He4_p_to_He4_He4_beta_pos_modified(T& rate_eval, const tf_t& tfactors, const amrex::Real log_scor, const amrex::Real dlog_scor_dT) {

        // He3 + He4 + p --> He4 + He4 (calls the underlying rate)
        // represents the sequence: He4(He3,γ)Be7(e⁻,ν)Li7(p,α)He4

        rate_He3_He4_to_Be7_reaclib(rate_eval, tfactors, log_scor, dlog_scor_dT);
        rate_eval.screened_rates(k_He3_He4_p_to_He4_He4_beta_pos_modified) = rate_eval.screened_rates(k_He3_He4_to_Be7_reaclib);
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_He3_He4_p_to_He4_He4_beta_pos_modified) = rate_eval.dscreened_rates_dT(k_He3_He4_to_Be7_reaclib);
        }
    }


}
#endif
