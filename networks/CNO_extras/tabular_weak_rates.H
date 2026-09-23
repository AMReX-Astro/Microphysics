#ifndef TABULAR_WEAK_RATES_H
#define TABULAR_WEAK_RATES_H

#include <cmath>

#include <AMReX.H>

#include <actual_network.H>
#include <fundamental_constants.H>
#include <network_properties.H>
#include <rate_type.H>
#include <table_rates.H>

namespace tabular_weak_rates {

    using namespace Rates;
    using namespace Species;

    using namespace rate_tables;

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_F17_to_O17_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // F17 --> O17

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_F17_O17_meta, j_F17_O17_rhoy, j_F17_O17_temp, j_F17_O17_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_F17_to_O17_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_F17_to_O17_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(F17) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_O17_to_F17_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // O17 --> F17

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_O17_F17_meta, j_O17_F17_rhoy, j_O17_F17_temp, j_O17_F17_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_O17_to_F17_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_O17_to_F17_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(O17) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_F18_to_Ne18_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // F18 --> Ne18

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_F18_Ne18_meta, j_F18_Ne18_rhoy, j_F18_Ne18_temp, j_F18_Ne18_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_F18_to_Ne18_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_F18_to_Ne18_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(F18) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_F18_to_O18_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // F18 --> O18

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_F18_O18_meta, j_F18_O18_rhoy, j_F18_O18_temp, j_F18_O18_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_F18_to_O18_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_F18_to_O18_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(F18) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ne18_to_F18_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ne18 --> F18

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ne18_F18_meta, j_Ne18_F18_rhoy, j_Ne18_F18_temp, j_Ne18_F18_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ne18_to_F18_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne18_to_F18_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ne18) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_O18_to_F18_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // O18 --> F18

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_O18_F18_meta, j_O18_F18_rhoy, j_O18_F18_temp, j_O18_F18_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_O18_to_F18_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_O18_to_F18_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(O18) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_F19_to_Ne19_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // F19 --> Ne19

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_F19_Ne19_meta, j_F19_Ne19_rhoy, j_F19_Ne19_temp, j_F19_Ne19_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_F19_to_Ne19_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_F19_to_Ne19_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(F19) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ne19_to_F19_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ne19 --> F19

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ne19_F19_meta, j_Ne19_F19_rhoy, j_Ne19_F19_temp, j_Ne19_F19_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ne19_to_F19_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne19_to_F19_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ne19) * (table_values.enu + table_values.gamma);
    }


    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void
    fill_rates([[maybe_unused]] const amrex::Real temp,
               [[maybe_unused]] const amrex::Real rhoy,
               [[maybe_unused]] const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
               [[maybe_unused]] T& rate_eval)
    {

        // Calculate tabular weak rates

        const amrex::Real log_temp = std::log10(temp);
        const amrex::Real log_rhoy = std::log10(rhoy);

        rate_eval.enuc_weak = 0.0_rt;

        {
            // F17_to_O17_electron_capture_weaktab
            rate_F17_to_O17_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // O17_to_F17_beta_neg_weaktab
            rate_O17_to_F17_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // F18_to_Ne18_beta_neg_weaktab
            rate_F18_to_Ne18_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // F18_to_O18_electron_capture_weaktab
            rate_F18_to_O18_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ne18_to_F18_electron_capture_weaktab
            rate_Ne18_to_F18_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // O18_to_F18_beta_neg_weaktab
            rate_O18_to_F18_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // F19_to_Ne19_beta_neg_weaktab
            rate_F19_to_Ne19_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ne19_to_F19_electron_capture_weaktab
            rate_Ne19_to_F19_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

    }

}
#endif
