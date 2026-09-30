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
    void rate_F20_to_Ne20_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // F20 --> Ne20

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_F20_Ne20_meta, j_F20_Ne20_rhoy, j_F20_Ne20_temp, j_F20_Ne20_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_F20_to_Ne20_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_F20_to_Ne20_beta_neg_weaktab) = table_values.drate_dT;
            rate_eval.denuc_weak_dY(F20) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(F20) * table_values.denu_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(F20) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_F20_to_O20_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // F20 --> O20

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_F20_O20_meta, j_F20_O20_rhoy, j_F20_O20_temp, j_F20_O20_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_F20_to_O20_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_F20_to_O20_electron_capture_weaktab) = table_values.drate_dT;
            rate_eval.denuc_weak_dY(F20) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(F20) * table_values.denu_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(F20) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ne20_to_F20_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ne20 --> F20

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ne20_F20_meta, j_Ne20_F20_rhoy, j_Ne20_F20_temp, j_Ne20_F20_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ne20_to_F20_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne20_to_F20_electron_capture_weaktab) = table_values.drate_dT;
            rate_eval.denuc_weak_dY(Ne20) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ne20) * table_values.denu_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ne20) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_O20_to_F20_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // O20 --> F20

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_O20_F20_meta, j_O20_F20_rhoy, j_O20_F20_temp, j_O20_F20_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_O20_to_F20_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_O20_to_F20_beta_neg_weaktab) = table_values.drate_dT;
            rate_eval.denuc_weak_dY(O20) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(O20) * table_values.denu_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(O20) * (table_values.enu + table_values.gamma);
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

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.denuc_weak_dT = 0.0_rt;
            for (int n = 1; n <= NumSpec; ++n) {
                rate_eval.denuc_weak_dY(n) = 0.0_rt;
            }
        }

        {
            // F20_to_Ne20_beta_neg_weaktab
            rate_F20_to_Ne20_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // F20_to_O20_electron_capture_weaktab
            rate_F20_to_O20_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ne20_to_F20_electron_capture_weaktab
            rate_Ne20_to_F20_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // O20_to_F20_beta_neg_weaktab
            rate_O20_to_F20_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

    }

}
#endif
