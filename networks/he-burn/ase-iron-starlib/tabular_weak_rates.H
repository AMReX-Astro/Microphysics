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

        tabular_evaluate(j_Co56_Fe56_meta, j_Co56_Fe56_rhoy, j_Co56_Fe56_temp, j_Co56_Fe56_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Co56_to_Fe56_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co56_to_Fe56_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Co56) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Co56_Ni56_meta, j_Co56_Ni56_rhoy, j_Co56_Ni56_temp, j_Co56_Ni56_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Co56_to_Ni56_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co56_to_Ni56_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Co56) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Fe56_Co56_meta, j_Fe56_Co56_rhoy, j_Fe56_Co56_temp, j_Fe56_Co56_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Fe56_to_Co56_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe56_to_Co56_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Fe56) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ni56_Co56_meta, j_Ni56_Co56_rhoy, j_Ni56_Co56_temp, j_Ni56_Co56_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ni56_to_Co56_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni56_to_Co56_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ni56) * (edot_nu + edot_gamma);

        tabular_evaluate(j_n_p_meta, j_n_p_rhoy, j_n_p_temp, j_n_p_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_n_to_p_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_n_to_p_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(N) * (edot_nu + edot_gamma);

        tabular_evaluate(j_p_n_meta, j_p_n_rhoy, j_p_n_temp, j_p_n_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_p_to_n_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_to_n_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(H1) * (edot_nu + edot_gamma);

    }

}
#endif
