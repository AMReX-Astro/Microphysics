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

        tabular_evaluate(j_Na21_Ne21_meta, j_Na21_Ne21_rhoy, j_Na21_Ne21_temp, j_Na21_Ne21_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Na21_to_Ne21_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Na21_to_Ne21_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Na21) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ne21_Na21_meta, j_Ne21_Na21_rhoy, j_Ne21_Na21_temp, j_Ne21_Na21_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ne21_to_Na21_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne21_to_Na21_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ne21) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Na22_Ne22_meta, j_Na22_Ne22_rhoy, j_Na22_Ne22_temp, j_Na22_Ne22_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Na22_to_Ne22_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Na22_to_Ne22_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Na22) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ne22_Na22_meta, j_Ne22_Na22_rhoy, j_Ne22_Na22_temp, j_Ne22_Na22_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ne22_to_Na22_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne22_to_Na22_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ne22) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Mg23_Na23_meta, j_Mg23_Na23_rhoy, j_Mg23_Na23_temp, j_Mg23_Na23_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Mg23_to_Na23_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg23_to_Na23_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Mg23) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Na23_Mg23_meta, j_Na23_Mg23_rhoy, j_Na23_Mg23_temp, j_Na23_Mg23_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Na23_to_Mg23_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Na23_to_Mg23_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Na23) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Al25_Mg25_meta, j_Al25_Mg25_rhoy, j_Al25_Mg25_temp, j_Al25_Mg25_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Al25_to_Mg25_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Al25_to_Mg25_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Al25) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Mg25_Al25_meta, j_Mg25_Al25_rhoy, j_Mg25_Al25_temp, j_Mg25_Al25_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Mg25_to_Al25_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg25_to_Al25_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Mg25) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Al26_Mg26_meta, j_Al26_Mg26_rhoy, j_Al26_Mg26_temp, j_Al26_Mg26_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Al26_to_Mg26_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Al26_to_Mg26_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Al26) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Mg26_Al26_meta, j_Mg26_Al26_rhoy, j_Mg26_Al26_temp, j_Mg26_Al26_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Mg26_to_Al26_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg26_to_Al26_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Mg26) * (edot_nu + edot_gamma);

        tabular_evaluate(j_P29_Si29_meta, j_P29_Si29_rhoy, j_P29_Si29_temp, j_P29_Si29_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_P29_to_Si29_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P29_to_Si29_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(P29) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Si29_P29_meta, j_Si29_P29_rhoy, j_Si29_P29_temp, j_Si29_P29_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Si29_to_P29_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si29_to_P29_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Si29) * (edot_nu + edot_gamma);

        tabular_evaluate(j_P30_Si30_meta, j_P30_Si30_rhoy, j_P30_Si30_temp, j_P30_Si30_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_P30_to_Si30_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P30_to_Si30_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(P30) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Si30_P30_meta, j_Si30_P30_rhoy, j_Si30_P30_temp, j_Si30_P30_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Si30_to_P30_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si30_to_P30_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Si30) * (edot_nu + edot_gamma);

        tabular_evaluate(j_P31_Si31_meta, j_P31_Si31_rhoy, j_P31_Si31_temp, j_P31_Si31_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_P31_to_Si31_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P31_to_Si31_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(P31) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Si31_P31_meta, j_Si31_P31_rhoy, j_Si31_P31_temp, j_Si31_P31_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Si31_to_P31_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si31_to_P31_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Si31) * (edot_nu + edot_gamma);

        tabular_evaluate(j_P32_S32_meta, j_P32_S32_rhoy, j_P32_S32_temp, j_P32_S32_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_P32_to_S32_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P32_to_S32_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(P32) * (edot_nu + edot_gamma);

        tabular_evaluate(j_P32_Si32_meta, j_P32_Si32_rhoy, j_P32_Si32_temp, j_P32_Si32_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_P32_to_Si32_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P32_to_Si32_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(P32) * (edot_nu + edot_gamma);

        tabular_evaluate(j_S32_P32_meta, j_S32_P32_rhoy, j_S32_P32_temp, j_S32_P32_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_S32_to_P32_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S32_to_P32_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(S32) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Si32_P32_meta, j_Si32_P32_rhoy, j_Si32_P32_temp, j_Si32_P32_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Si32_to_P32_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si32_to_P32_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Si32) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cl33_S33_meta, j_Cl33_S33_rhoy, j_Cl33_S33_temp, j_Cl33_S33_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cl33_to_S33_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl33_to_S33_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cl33) * (edot_nu + edot_gamma);

        tabular_evaluate(j_P33_S33_meta, j_P33_S33_rhoy, j_P33_S33_temp, j_P33_S33_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_P33_to_S33_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P33_to_S33_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(P33) * (edot_nu + edot_gamma);

        tabular_evaluate(j_S33_Cl33_meta, j_S33_Cl33_rhoy, j_S33_Cl33_temp, j_S33_Cl33_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_S33_to_Cl33_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S33_to_Cl33_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(S33) * (edot_nu + edot_gamma);

        tabular_evaluate(j_S33_P33_meta, j_S33_P33_rhoy, j_S33_P33_temp, j_S33_P33_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_S33_to_P33_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S33_to_P33_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(S33) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cl34_S34_meta, j_Cl34_S34_rhoy, j_Cl34_S34_temp, j_Cl34_S34_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cl34_to_S34_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl34_to_S34_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cl34) * (edot_nu + edot_gamma);

        tabular_evaluate(j_S34_Cl34_meta, j_S34_Cl34_rhoy, j_S34_Cl34_temp, j_S34_Cl34_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_S34_to_Cl34_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S34_to_Cl34_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(S34) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cl35_S35_meta, j_Cl35_S35_rhoy, j_Cl35_S35_temp, j_Cl35_S35_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cl35_to_S35_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl35_to_S35_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cl35) * (edot_nu + edot_gamma);

        tabular_evaluate(j_S35_Cl35_meta, j_S35_Cl35_rhoy, j_S35_Cl35_temp, j_S35_Cl35_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_S35_to_Cl35_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S35_to_Cl35_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(S35) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ar36_Cl36_meta, j_Ar36_Cl36_rhoy, j_Ar36_Cl36_temp, j_Ar36_Cl36_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ar36_to_Cl36_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar36_to_Cl36_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ar36) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cl36_Ar36_meta, j_Cl36_Ar36_rhoy, j_Cl36_Ar36_temp, j_Cl36_Ar36_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cl36_to_Ar36_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl36_to_Ar36_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cl36) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cl36_S36_meta, j_Cl36_S36_rhoy, j_Cl36_S36_temp, j_Cl36_S36_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cl36_to_S36_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl36_to_S36_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cl36) * (edot_nu + edot_gamma);

        tabular_evaluate(j_S36_Cl36_meta, j_S36_Cl36_rhoy, j_S36_Cl36_temp, j_S36_Cl36_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_S36_to_Cl36_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S36_to_Cl36_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(S36) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ar37_Cl37_meta, j_Ar37_Cl37_rhoy, j_Ar37_Cl37_temp, j_Ar37_Cl37_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ar37_to_Cl37_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar37_to_Cl37_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ar37) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ar37_K37_meta, j_Ar37_K37_rhoy, j_Ar37_K37_temp, j_Ar37_K37_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ar37_to_K37_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar37_to_K37_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ar37) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cl37_Ar37_meta, j_Cl37_Ar37_rhoy, j_Cl37_Ar37_temp, j_Cl37_Ar37_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cl37_to_Ar37_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl37_to_Ar37_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cl37) * (edot_nu + edot_gamma);

        tabular_evaluate(j_K37_Ar37_meta, j_K37_Ar37_rhoy, j_K37_Ar37_temp, j_K37_Ar37_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_K37_to_Ar37_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K37_to_Ar37_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(K37) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ar38_K38_meta, j_Ar38_K38_rhoy, j_Ar38_K38_temp, j_Ar38_K38_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ar38_to_K38_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar38_to_K38_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ar38) * (edot_nu + edot_gamma);

        tabular_evaluate(j_K38_Ar38_meta, j_K38_Ar38_rhoy, j_K38_Ar38_temp, j_K38_Ar38_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_K38_to_Ar38_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K38_to_Ar38_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(K38) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ar39_K39_meta, j_Ar39_K39_rhoy, j_Ar39_K39_temp, j_Ar39_K39_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ar39_to_K39_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar39_to_K39_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ar39) * (edot_nu + edot_gamma);

        tabular_evaluate(j_K39_Ar39_meta, j_K39_Ar39_rhoy, j_K39_Ar39_temp, j_K39_Ar39_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_K39_to_Ar39_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K39_to_Ar39_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(K39) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ar40_K40_meta, j_Ar40_K40_rhoy, j_Ar40_K40_temp, j_Ar40_K40_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ar40_to_K40_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar40_to_K40_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ar40) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ca40_K40_meta, j_Ca40_K40_rhoy, j_Ca40_K40_temp, j_Ca40_K40_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ca40_to_K40_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca40_to_K40_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ca40) * (edot_nu + edot_gamma);

        tabular_evaluate(j_K40_Ar40_meta, j_K40_Ar40_rhoy, j_K40_Ar40_temp, j_K40_Ar40_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_K40_to_Ar40_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K40_to_Ar40_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(K40) * (edot_nu + edot_gamma);

        tabular_evaluate(j_K40_Ca40_meta, j_K40_Ca40_rhoy, j_K40_Ca40_temp, j_K40_Ca40_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_K40_to_Ca40_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K40_to_Ca40_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(K40) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ca41_K41_meta, j_Ca41_K41_rhoy, j_Ca41_K41_temp, j_Ca41_K41_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ca41_to_K41_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca41_to_K41_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ca41) * (edot_nu + edot_gamma);

        tabular_evaluate(j_K41_Ca41_meta, j_K41_Ca41_rhoy, j_K41_Ca41_temp, j_K41_Ca41_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_K41_to_Ca41_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K41_to_Ca41_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(K41) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ca43_Sc43_meta, j_Ca43_Sc43_rhoy, j_Ca43_Sc43_temp, j_Ca43_Sc43_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ca43_to_Sc43_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca43_to_Sc43_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ca43) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Sc43_Ca43_meta, j_Sc43_Ca43_rhoy, j_Sc43_Ca43_temp, j_Sc43_Ca43_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Sc43_to_Ca43_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc43_to_Ca43_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Sc43) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ca44_Sc44_meta, j_Ca44_Sc44_rhoy, j_Ca44_Sc44_temp, j_Ca44_Sc44_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ca44_to_Sc44_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca44_to_Sc44_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ca44) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Sc44_Ca44_meta, j_Sc44_Ca44_rhoy, j_Sc44_Ca44_temp, j_Sc44_Ca44_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Sc44_to_Ca44_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc44_to_Ca44_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Sc44) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Sc44_Ti44_meta, j_Sc44_Ti44_rhoy, j_Sc44_Ti44_temp, j_Sc44_Ti44_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Sc44_to_Ti44_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc44_to_Ti44_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Sc44) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ti44_Sc44_meta, j_Ti44_Sc44_rhoy, j_Ti44_Sc44_temp, j_Ti44_Sc44_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ti44_to_Sc44_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti44_to_Sc44_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ti44) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Co53_Fe53_meta, j_Co53_Fe53_rhoy, j_Co53_Fe53_temp, j_Co53_Fe53_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Co53_to_Fe53_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co53_to_Fe53_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Co53) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Fe53_Co53_meta, j_Fe53_Co53_rhoy, j_Fe53_Co53_temp, j_Fe53_Co53_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Fe53_to_Co53_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe53_to_Co53_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Fe53) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cu57_Ni57_meta, j_Cu57_Ni57_rhoy, j_Cu57_Ni57_temp, j_Cu57_Ni57_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cu57_to_Ni57_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu57_to_Ni57_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cu57) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ni57_Cu57_meta, j_Ni57_Cu57_rhoy, j_Ni57_Cu57_temp, j_Ni57_Cu57_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ni57_to_Cu57_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni57_to_Cu57_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ni57) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ca45_Sc45_meta, j_Ca45_Sc45_rhoy, j_Ca45_Sc45_temp, j_Ca45_Sc45_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ca45_to_Sc45_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca45_to_Sc45_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ca45) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Sc45_Ca45_meta, j_Sc45_Ca45_rhoy, j_Sc45_Ca45_temp, j_Sc45_Ca45_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Sc45_to_Ca45_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc45_to_Ca45_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Sc45) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Sc45_Ti45_meta, j_Sc45_Ti45_rhoy, j_Sc45_Ti45_temp, j_Sc45_Ti45_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Sc45_to_Ti45_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc45_to_Ti45_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Sc45) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ti45_Sc45_meta, j_Ti45_Sc45_rhoy, j_Ti45_Sc45_temp, j_Ti45_Sc45_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ti45_to_Sc45_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti45_to_Sc45_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ti45) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ca46_Sc46_meta, j_Ca46_Sc46_rhoy, j_Ca46_Sc46_temp, j_Ca46_Sc46_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ca46_to_Sc46_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca46_to_Sc46_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ca46) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Sc46_Ca46_meta, j_Sc46_Ca46_rhoy, j_Sc46_Ca46_temp, j_Sc46_Ca46_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Sc46_to_Ca46_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc46_to_Ca46_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Sc46) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Sc46_Ti46_meta, j_Sc46_Ti46_rhoy, j_Sc46_Ti46_temp, j_Sc46_Ti46_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Sc46_to_Ti46_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc46_to_Ti46_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Sc46) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ti46_Sc46_meta, j_Ti46_Sc46_rhoy, j_Ti46_Sc46_temp, j_Ti46_Sc46_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ti46_to_Sc46_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti46_to_Sc46_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ti46) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ti46_V46_meta, j_Ti46_V46_rhoy, j_Ti46_V46_temp, j_Ti46_V46_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ti46_to_V46_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti46_to_V46_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ti46) * (edot_nu + edot_gamma);

        tabular_evaluate(j_V46_Ti46_meta, j_V46_Ti46_rhoy, j_V46_Ti46_temp, j_V46_Ti46_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_V46_to_Ti46_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V46_to_Ti46_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(V46) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ca47_Sc47_meta, j_Ca47_Sc47_rhoy, j_Ca47_Sc47_temp, j_Ca47_Sc47_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ca47_to_Sc47_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca47_to_Sc47_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ca47) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Sc47_Ca47_meta, j_Sc47_Ca47_rhoy, j_Sc47_Ca47_temp, j_Sc47_Ca47_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Sc47_to_Ca47_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc47_to_Ca47_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Sc47) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Sc47_Ti47_meta, j_Sc47_Ti47_rhoy, j_Sc47_Ti47_temp, j_Sc47_Ti47_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Sc47_to_Ti47_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc47_to_Ti47_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Sc47) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ti47_Sc47_meta, j_Ti47_Sc47_rhoy, j_Ti47_Sc47_temp, j_Ti47_Sc47_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ti47_to_Sc47_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti47_to_Sc47_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ti47) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ti47_V47_meta, j_Ti47_V47_rhoy, j_Ti47_V47_temp, j_Ti47_V47_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ti47_to_V47_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti47_to_V47_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ti47) * (edot_nu + edot_gamma);

        tabular_evaluate(j_V47_Ti47_meta, j_V47_Ti47_rhoy, j_V47_Ti47_temp, j_V47_Ti47_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_V47_to_Ti47_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V47_to_Ti47_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(V47) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ca48_Sc48_meta, j_Ca48_Sc48_rhoy, j_Ca48_Sc48_temp, j_Ca48_Sc48_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ca48_to_Sc48_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca48_to_Sc48_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ca48) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cr48_V48_meta, j_Cr48_V48_rhoy, j_Cr48_V48_temp, j_Cr48_V48_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cr48_to_V48_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr48_to_V48_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cr48) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Sc48_Ca48_meta, j_Sc48_Ca48_rhoy, j_Sc48_Ca48_temp, j_Sc48_Ca48_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Sc48_to_Ca48_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc48_to_Ca48_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Sc48) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Sc48_Ti48_meta, j_Sc48_Ti48_rhoy, j_Sc48_Ti48_temp, j_Sc48_Ti48_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Sc48_to_Ti48_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc48_to_Ti48_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Sc48) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ti48_Sc48_meta, j_Ti48_Sc48_rhoy, j_Ti48_Sc48_temp, j_Ti48_Sc48_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ti48_to_Sc48_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti48_to_Sc48_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ti48) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ti48_V48_meta, j_Ti48_V48_rhoy, j_Ti48_V48_temp, j_Ti48_V48_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ti48_to_V48_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti48_to_V48_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ti48) * (edot_nu + edot_gamma);

        tabular_evaluate(j_V48_Cr48_meta, j_V48_Cr48_rhoy, j_V48_Cr48_temp, j_V48_Cr48_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_V48_to_Cr48_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V48_to_Cr48_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(V48) * (edot_nu + edot_gamma);

        tabular_evaluate(j_V48_Ti48_meta, j_V48_Ti48_rhoy, j_V48_Ti48_temp, j_V48_Ti48_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_V48_to_Ti48_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V48_to_Ti48_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(V48) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cr49_V49_meta, j_Cr49_V49_rhoy, j_Cr49_V49_temp, j_Cr49_V49_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cr49_to_V49_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr49_to_V49_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cr49) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Sc49_Ti49_meta, j_Sc49_Ti49_rhoy, j_Sc49_Ti49_temp, j_Sc49_Ti49_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Sc49_to_Ti49_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc49_to_Ti49_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Sc49) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ti49_Sc49_meta, j_Ti49_Sc49_rhoy, j_Ti49_Sc49_temp, j_Ti49_Sc49_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ti49_to_Sc49_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti49_to_Sc49_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ti49) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ti49_V49_meta, j_Ti49_V49_rhoy, j_Ti49_V49_temp, j_Ti49_V49_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ti49_to_V49_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti49_to_V49_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ti49) * (edot_nu + edot_gamma);

        tabular_evaluate(j_V49_Cr49_meta, j_V49_Cr49_rhoy, j_V49_Cr49_temp, j_V49_Cr49_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_V49_to_Cr49_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V49_to_Cr49_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(V49) * (edot_nu + edot_gamma);

        tabular_evaluate(j_V49_Ti49_meta, j_V49_Ti49_rhoy, j_V49_Ti49_temp, j_V49_Ti49_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_V49_to_Ti49_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V49_to_Ti49_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(V49) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cr50_Mn50_meta, j_Cr50_Mn50_rhoy, j_Cr50_Mn50_temp, j_Cr50_Mn50_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cr50_to_Mn50_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr50_to_Mn50_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cr50) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cr50_V50_meta, j_Cr50_V50_rhoy, j_Cr50_V50_temp, j_Cr50_V50_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cr50_to_V50_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr50_to_V50_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cr50) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Mn50_Cr50_meta, j_Mn50_Cr50_rhoy, j_Mn50_Cr50_temp, j_Mn50_Cr50_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Mn50_to_Cr50_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn50_to_Cr50_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Mn50) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ti50_V50_meta, j_Ti50_V50_rhoy, j_Ti50_V50_temp, j_Ti50_V50_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ti50_to_V50_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti50_to_V50_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ti50) * (edot_nu + edot_gamma);

        tabular_evaluate(j_V50_Cr50_meta, j_V50_Cr50_rhoy, j_V50_Cr50_temp, j_V50_Cr50_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_V50_to_Cr50_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V50_to_Cr50_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(V50) * (edot_nu + edot_gamma);

        tabular_evaluate(j_V50_Ti50_meta, j_V50_Ti50_rhoy, j_V50_Ti50_temp, j_V50_Ti50_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_V50_to_Ti50_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V50_to_Ti50_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(V50) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cr51_Mn51_meta, j_Cr51_Mn51_rhoy, j_Cr51_Mn51_temp, j_Cr51_Mn51_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cr51_to_Mn51_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr51_to_Mn51_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cr51) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cr51_V51_meta, j_Cr51_V51_rhoy, j_Cr51_V51_temp, j_Cr51_V51_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cr51_to_V51_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr51_to_V51_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cr51) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Mn51_Cr51_meta, j_Mn51_Cr51_rhoy, j_Mn51_Cr51_temp, j_Mn51_Cr51_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Mn51_to_Cr51_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn51_to_Cr51_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Mn51) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ti51_V51_meta, j_Ti51_V51_rhoy, j_Ti51_V51_temp, j_Ti51_V51_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ti51_to_V51_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti51_to_V51_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ti51) * (edot_nu + edot_gamma);

        tabular_evaluate(j_V51_Cr51_meta, j_V51_Cr51_rhoy, j_V51_Cr51_temp, j_V51_Cr51_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_V51_to_Cr51_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V51_to_Cr51_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(V51) * (edot_nu + edot_gamma);

        tabular_evaluate(j_V51_Ti51_meta, j_V51_Ti51_rhoy, j_V51_Ti51_temp, j_V51_Ti51_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_V51_to_Ti51_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V51_to_Ti51_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(V51) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cr52_Mn52_meta, j_Cr52_Mn52_rhoy, j_Cr52_Mn52_temp, j_Cr52_Mn52_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cr52_to_Mn52_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr52_to_Mn52_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cr52) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cr52_V52_meta, j_Cr52_V52_rhoy, j_Cr52_V52_temp, j_Cr52_V52_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cr52_to_V52_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr52_to_V52_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cr52) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Fe52_Mn52_meta, j_Fe52_Mn52_rhoy, j_Fe52_Mn52_temp, j_Fe52_Mn52_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Fe52_to_Mn52_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe52_to_Mn52_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Fe52) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Mn52_Cr52_meta, j_Mn52_Cr52_rhoy, j_Mn52_Cr52_temp, j_Mn52_Cr52_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Mn52_to_Cr52_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn52_to_Cr52_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Mn52) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Mn52_Fe52_meta, j_Mn52_Fe52_rhoy, j_Mn52_Fe52_temp, j_Mn52_Fe52_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Mn52_to_Fe52_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn52_to_Fe52_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Mn52) * (edot_nu + edot_gamma);

        tabular_evaluate(j_V52_Cr52_meta, j_V52_Cr52_rhoy, j_V52_Cr52_temp, j_V52_Cr52_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_V52_to_Cr52_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V52_to_Cr52_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(V52) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cr53_Mn53_meta, j_Cr53_Mn53_rhoy, j_Cr53_Mn53_temp, j_Cr53_Mn53_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cr53_to_Mn53_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr53_to_Mn53_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cr53) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Fe53_Mn53_meta, j_Fe53_Mn53_rhoy, j_Fe53_Mn53_temp, j_Fe53_Mn53_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Fe53_to_Mn53_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe53_to_Mn53_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Fe53) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Mn53_Cr53_meta, j_Mn53_Cr53_rhoy, j_Mn53_Cr53_temp, j_Mn53_Cr53_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Mn53_to_Cr53_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn53_to_Cr53_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Mn53) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Mn53_Fe53_meta, j_Mn53_Fe53_rhoy, j_Mn53_Fe53_temp, j_Mn53_Fe53_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Mn53_to_Fe53_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn53_to_Fe53_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Mn53) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Co54_Fe54_meta, j_Co54_Fe54_rhoy, j_Co54_Fe54_temp, j_Co54_Fe54_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Co54_to_Fe54_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co54_to_Fe54_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Co54) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cr54_Mn54_meta, j_Cr54_Mn54_rhoy, j_Cr54_Mn54_temp, j_Cr54_Mn54_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cr54_to_Mn54_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr54_to_Mn54_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cr54) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Fe54_Co54_meta, j_Fe54_Co54_rhoy, j_Fe54_Co54_temp, j_Fe54_Co54_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Fe54_to_Co54_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe54_to_Co54_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Fe54) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Fe54_Mn54_meta, j_Fe54_Mn54_rhoy, j_Fe54_Mn54_temp, j_Fe54_Mn54_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Fe54_to_Mn54_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe54_to_Mn54_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Fe54) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Mn54_Cr54_meta, j_Mn54_Cr54_rhoy, j_Mn54_Cr54_temp, j_Mn54_Cr54_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Mn54_to_Cr54_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn54_to_Cr54_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Mn54) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Mn54_Fe54_meta, j_Mn54_Fe54_rhoy, j_Mn54_Fe54_temp, j_Mn54_Fe54_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Mn54_to_Fe54_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn54_to_Fe54_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Mn54) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Co55_Fe55_meta, j_Co55_Fe55_rhoy, j_Co55_Fe55_temp, j_Co55_Fe55_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Co55_to_Fe55_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co55_to_Fe55_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Co55) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Fe55_Co55_meta, j_Fe55_Co55_rhoy, j_Fe55_Co55_temp, j_Fe55_Co55_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Fe55_to_Co55_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe55_to_Co55_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Fe55) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Fe55_Mn55_meta, j_Fe55_Mn55_rhoy, j_Fe55_Mn55_temp, j_Fe55_Mn55_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Fe55_to_Mn55_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe55_to_Mn55_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Fe55) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Mn55_Fe55_meta, j_Mn55_Fe55_rhoy, j_Mn55_Fe55_temp, j_Mn55_Fe55_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Mn55_to_Fe55_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn55_to_Fe55_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Mn55) * (edot_nu + edot_gamma);

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

        tabular_evaluate(j_Co57_Fe57_meta, j_Co57_Fe57_rhoy, j_Co57_Fe57_temp, j_Co57_Fe57_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Co57_to_Fe57_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co57_to_Fe57_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Co57) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Co57_Ni57_meta, j_Co57_Ni57_rhoy, j_Co57_Ni57_temp, j_Co57_Ni57_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Co57_to_Ni57_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co57_to_Ni57_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Co57) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Fe57_Co57_meta, j_Fe57_Co57_rhoy, j_Fe57_Co57_temp, j_Fe57_Co57_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Fe57_to_Co57_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe57_to_Co57_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Fe57) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ni57_Co57_meta, j_Ni57_Co57_rhoy, j_Ni57_Co57_temp, j_Ni57_Co57_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ni57_to_Co57_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni57_to_Co57_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ni57) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Co58_Fe58_meta, j_Co58_Fe58_rhoy, j_Co58_Fe58_temp, j_Co58_Fe58_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Co58_to_Fe58_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co58_to_Fe58_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Co58) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Co58_Ni58_meta, j_Co58_Ni58_rhoy, j_Co58_Ni58_temp, j_Co58_Ni58_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Co58_to_Ni58_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co58_to_Ni58_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Co58) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cu58_Ni58_meta, j_Cu58_Ni58_rhoy, j_Cu58_Ni58_temp, j_Cu58_Ni58_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cu58_to_Ni58_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu58_to_Ni58_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cu58) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Fe58_Co58_meta, j_Fe58_Co58_rhoy, j_Fe58_Co58_temp, j_Fe58_Co58_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Fe58_to_Co58_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe58_to_Co58_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Fe58) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ni58_Co58_meta, j_Ni58_Co58_rhoy, j_Ni58_Co58_temp, j_Ni58_Co58_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ni58_to_Co58_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni58_to_Co58_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ni58) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ni58_Cu58_meta, j_Ni58_Cu58_rhoy, j_Ni58_Cu58_temp, j_Ni58_Cu58_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ni58_to_Cu58_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni58_to_Cu58_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ni58) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Co59_Ni59_meta, j_Co59_Ni59_rhoy, j_Co59_Ni59_temp, j_Co59_Ni59_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Co59_to_Ni59_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co59_to_Ni59_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Co59) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cu59_Ni59_meta, j_Cu59_Ni59_rhoy, j_Cu59_Ni59_temp, j_Cu59_Ni59_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cu59_to_Ni59_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu59_to_Ni59_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cu59) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ni59_Co59_meta, j_Ni59_Co59_rhoy, j_Ni59_Co59_temp, j_Ni59_Co59_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ni59_to_Co59_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni59_to_Co59_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ni59) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ni59_Cu59_meta, j_Ni59_Cu59_rhoy, j_Ni59_Cu59_temp, j_Ni59_Cu59_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ni59_to_Cu59_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni59_to_Cu59_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ni59) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cu60_Ni60_meta, j_Cu60_Ni60_rhoy, j_Cu60_Ni60_temp, j_Cu60_Ni60_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cu60_to_Ni60_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu60_to_Ni60_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cu60) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cu60_Zn60_meta, j_Cu60_Zn60_rhoy, j_Cu60_Zn60_temp, j_Cu60_Zn60_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cu60_to_Zn60_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu60_to_Zn60_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cu60) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ni60_Cu60_meta, j_Ni60_Cu60_rhoy, j_Ni60_Cu60_temp, j_Ni60_Cu60_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ni60_to_Cu60_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni60_to_Cu60_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ni60) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Zn60_Cu60_meta, j_Zn60_Cu60_rhoy, j_Zn60_Cu60_temp, j_Zn60_Cu60_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Zn60_to_Cu60_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn60_to_Cu60_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Zn60) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cu61_Ni61_meta, j_Cu61_Ni61_rhoy, j_Cu61_Ni61_temp, j_Cu61_Ni61_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cu61_to_Ni61_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu61_to_Ni61_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cu61) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cu61_Zn61_meta, j_Cu61_Zn61_rhoy, j_Cu61_Zn61_temp, j_Cu61_Zn61_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cu61_to_Zn61_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu61_to_Zn61_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cu61) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ni61_Cu61_meta, j_Ni61_Cu61_rhoy, j_Ni61_Cu61_temp, j_Ni61_Cu61_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ni61_to_Cu61_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni61_to_Cu61_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ni61) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Zn61_Cu61_meta, j_Zn61_Cu61_rhoy, j_Zn61_Cu61_temp, j_Zn61_Cu61_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Zn61_to_Cu61_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn61_to_Cu61_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Zn61) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cu62_Ni62_meta, j_Cu62_Ni62_rhoy, j_Cu62_Ni62_temp, j_Cu62_Ni62_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cu62_to_Ni62_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu62_to_Ni62_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cu62) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cu62_Zn62_meta, j_Cu62_Zn62_rhoy, j_Cu62_Zn62_temp, j_Cu62_Zn62_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cu62_to_Zn62_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu62_to_Zn62_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cu62) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ga62_Zn62_meta, j_Ga62_Zn62_rhoy, j_Ga62_Zn62_temp, j_Ga62_Zn62_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ga62_to_Zn62_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ga62_to_Zn62_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ga62) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ni62_Cu62_meta, j_Ni62_Cu62_rhoy, j_Ni62_Cu62_temp, j_Ni62_Cu62_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ni62_to_Cu62_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni62_to_Cu62_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ni62) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Zn62_Cu62_meta, j_Zn62_Cu62_rhoy, j_Zn62_Cu62_temp, j_Zn62_Cu62_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Zn62_to_Cu62_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn62_to_Cu62_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Zn62) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Zn62_Ga62_meta, j_Zn62_Ga62_rhoy, j_Zn62_Ga62_temp, j_Zn62_Ga62_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Zn62_to_Ga62_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn62_to_Ga62_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Zn62) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cu63_Ni63_meta, j_Cu63_Ni63_rhoy, j_Cu63_Ni63_temp, j_Cu63_Ni63_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cu63_to_Ni63_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu63_to_Ni63_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cu63) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cu63_Zn63_meta, j_Cu63_Zn63_rhoy, j_Cu63_Zn63_temp, j_Cu63_Zn63_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cu63_to_Zn63_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu63_to_Zn63_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cu63) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ga63_Zn63_meta, j_Ga63_Zn63_rhoy, j_Ga63_Zn63_temp, j_Ga63_Zn63_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ga63_to_Zn63_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ga63_to_Zn63_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ga63) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ni63_Cu63_meta, j_Ni63_Cu63_rhoy, j_Ni63_Cu63_temp, j_Ni63_Cu63_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ni63_to_Cu63_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni63_to_Cu63_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ni63) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Zn63_Cu63_meta, j_Zn63_Cu63_rhoy, j_Zn63_Cu63_temp, j_Zn63_Cu63_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Zn63_to_Cu63_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn63_to_Cu63_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Zn63) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Zn63_Ga63_meta, j_Zn63_Ga63_rhoy, j_Zn63_Ga63_temp, j_Zn63_Ga63_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Zn63_to_Ga63_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn63_to_Ga63_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Zn63) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cu64_Ni64_meta, j_Cu64_Ni64_rhoy, j_Cu64_Ni64_temp, j_Cu64_Ni64_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cu64_to_Ni64_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu64_to_Ni64_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cu64) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cu64_Zn64_meta, j_Cu64_Zn64_rhoy, j_Cu64_Zn64_temp, j_Cu64_Zn64_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cu64_to_Zn64_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu64_to_Zn64_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cu64) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ga64_Ge64_meta, j_Ga64_Ge64_rhoy, j_Ga64_Ge64_temp, j_Ga64_Ge64_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ga64_to_Ge64_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ga64_to_Ge64_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ga64) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ga64_Zn64_meta, j_Ga64_Zn64_rhoy, j_Ga64_Zn64_temp, j_Ga64_Zn64_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ga64_to_Zn64_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ga64_to_Zn64_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ga64) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ge64_Ga64_meta, j_Ge64_Ga64_rhoy, j_Ge64_Ga64_temp, j_Ge64_Ga64_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ge64_to_Ga64_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ge64_to_Ga64_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ge64) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Ni64_Cu64_meta, j_Ni64_Cu64_rhoy, j_Ni64_Cu64_temp, j_Ni64_Cu64_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Ni64_to_Cu64_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni64_to_Cu64_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Ni64) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Zn64_Cu64_meta, j_Zn64_Cu64_rhoy, j_Zn64_Cu64_temp, j_Zn64_Cu64_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Zn64_to_Cu64_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn64_to_Cu64_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Zn64) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Zn64_Ga64_meta, j_Zn64_Ga64_rhoy, j_Zn64_Ga64_temp, j_Zn64_Ga64_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Zn64_to_Ga64_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn64_to_Ga64_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Zn64) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Cu65_Zn65_meta, j_Cu65_Zn65_rhoy, j_Cu65_Zn65_temp, j_Cu65_Zn65_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Cu65_to_Zn65_beta_neg_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu65_to_Zn65_beta_neg_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Cu65) * (edot_nu + edot_gamma);

        tabular_evaluate(j_Zn65_Cu65_meta, j_Zn65_Cu65_rhoy, j_Zn65_Cu65_temp, j_Zn65_Cu65_data,
                         log_rhoy, log_temp, temp, rate, drate_dt, edot_nu, edot_gamma);
        rate_eval.screened_rates(k_Zn65_to_Cu65_electron_capture_weaktab) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn65_to_Cu65_electron_capture_weaktab) = drate_dt;
        }
        rate_eval.enuc_weak += C::n_A * Y(Zn65) * (edot_nu + edot_gamma);

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
