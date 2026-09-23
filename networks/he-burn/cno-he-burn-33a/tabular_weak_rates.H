#ifndef TABULAR_WEAK_RATES_H
#define TABULAR_WEAK_RATES_H

#include <AMReX.H>

#include <table_rates.H>

namespace tabular_weak_rates {

    using namespace Rates;
    using namespace Species;

    using namespace rate_tables;

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void
    fill_rates([[maybe_unused]] const amrex::Real temp,
               [[maybe_unused]] const amrex::Real rhoy,
               [[maybe_unused]] const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
               [[maybe_unused]] T& rate_eval)
    {

        // Calculate tabular weak rates

        [[maybe_unused]] amrex::Real rate, drate_dt, edot_nu, edot_gamma;

        rate_eval.enuc_weak = 0.0_rt;

        const amrex::Real log_temp = std::log10(temp);
        const amrex::Real log_rhoy = std::log10(rhoy);

        tabular_evaluate(j_F17_O17_meta, j_F17_O17_rhoy, j_F17_O17_temp, j_F17_O17_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_F17_to_O17_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_F17_to_O17_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(F17) * (edot_nu + edot_gamma);

        tabular_evaluate(j_O17_F17_meta, j_O17_F17_rhoy, j_O17_F17_temp, j_O17_F17_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_O17_to_F17_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_O17_to_F17_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(O17) * (edot_nu + edot_gamma);

        tabular_evaluate(j_F18_Ne18_meta, j_F18_Ne18_rhoy, j_F18_Ne18_temp, j_F18_Ne18_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_F18_to_Ne18_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_F18_to_Ne18_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(F18) * (edot_nu + edot_gamma);

        tabular_evaluate(j_F18_O18_meta, j_F18_O18_rhoy, j_F18_O18_temp, j_F18_O18_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_F18_to_O18_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_F18_to_O18_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(F18) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ne18_F18_meta, j_Ne18_F18_rhoy, j_Ne18_F18_temp, j_Ne18_F18_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ne18_to_F18_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne18_to_F18_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ne18) * (edot_nu + edot_gamma);

        tabular_evaluate(j_O18_F18_meta, j_O18_F18_rhoy, j_O18_F18_temp, j_O18_F18_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_O18_to_F18_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_O18_to_F18_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(O18) * (edot_nu + edot_gamma);

        tabular_evaluate(j_F19_Ne19_meta, j_F19_Ne19_rhoy, j_F19_Ne19_temp, j_F19_Ne19_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_F19_to_Ne19_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_F19_to_Ne19_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(F19) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ne19_F19_meta, j_Ne19_F19_rhoy, j_Ne19_F19_temp, j_Ne19_F19_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ne19_to_F19_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne19_to_F19_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ne19) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Mg22_Na22_meta, j_Mg22_Na22_rhoy, j_Mg22_Na22_temp, j_Mg22_Na22_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Mg22_to_Na22_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg22_to_Na22_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Mg22) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Na22_Mg22_meta, j_Na22_Mg22_rhoy, j_Na22_Mg22_temp, j_Na22_Mg22_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Na22_to_Mg22_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Na22_to_Mg22_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Na22) * (edot_nu + edot_gamma);

    }

}
#endif
