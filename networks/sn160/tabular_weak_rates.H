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
    void rate_Na21_to_Ne21_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Na21 --> Ne21

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Na21_Ne21_meta, j_Na21_Ne21_rhoy, j_Na21_Ne21_temp, j_Na21_Ne21_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Na21_to_Ne21_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Na21_to_Ne21_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Na21) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ne21_to_Na21_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ne21 --> Na21

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ne21_Na21_meta, j_Ne21_Na21_rhoy, j_Ne21_Na21_temp, j_Ne21_Na21_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ne21_to_Na21_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne21_to_Na21_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ne21) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Na22_to_Ne22_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Na22 --> Ne22

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Na22_Ne22_meta, j_Na22_Ne22_rhoy, j_Na22_Ne22_temp, j_Na22_Ne22_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Na22_to_Ne22_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Na22_to_Ne22_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Na22) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ne22_to_Na22_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ne22 --> Na22

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ne22_Na22_meta, j_Ne22_Na22_rhoy, j_Ne22_Na22_temp, j_Ne22_Na22_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ne22_to_Na22_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne22_to_Na22_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ne22) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mg23_to_Na23_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Mg23 --> Na23

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Mg23_Na23_meta, j_Mg23_Na23_rhoy, j_Mg23_Na23_temp, j_Mg23_Na23_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Mg23_to_Na23_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg23_to_Na23_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Mg23) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Na23_to_Mg23_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Na23 --> Mg23

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Na23_Mg23_meta, j_Na23_Mg23_rhoy, j_Na23_Mg23_temp, j_Na23_Mg23_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Na23_to_Mg23_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Na23_to_Mg23_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Na23) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Al25_to_Mg25_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Al25 --> Mg25

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Al25_Mg25_meta, j_Al25_Mg25_rhoy, j_Al25_Mg25_temp, j_Al25_Mg25_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Al25_to_Mg25_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Al25_to_Mg25_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Al25) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mg25_to_Al25_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Mg25 --> Al25

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Mg25_Al25_meta, j_Mg25_Al25_rhoy, j_Mg25_Al25_temp, j_Mg25_Al25_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Mg25_to_Al25_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg25_to_Al25_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Mg25) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Al26_to_Mg26_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Al26 --> Mg26

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Al26_Mg26_meta, j_Al26_Mg26_rhoy, j_Al26_Mg26_temp, j_Al26_Mg26_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Al26_to_Mg26_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Al26_to_Mg26_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Al26) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mg26_to_Al26_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Mg26 --> Al26

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Mg26_Al26_meta, j_Mg26_Al26_rhoy, j_Mg26_Al26_temp, j_Mg26_Al26_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Mg26_to_Al26_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg26_to_Al26_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Mg26) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P29_to_Si29_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // P29 --> Si29

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_P29_Si29_meta, j_P29_Si29_rhoy, j_P29_Si29_temp, j_P29_Si29_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_P29_to_Si29_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P29_to_Si29_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(P29) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si29_to_P29_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Si29 --> P29

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Si29_P29_meta, j_Si29_P29_rhoy, j_Si29_P29_temp, j_Si29_P29_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Si29_to_P29_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si29_to_P29_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Si29) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P30_to_Si30_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // P30 --> Si30

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_P30_Si30_meta, j_P30_Si30_rhoy, j_P30_Si30_temp, j_P30_Si30_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_P30_to_Si30_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P30_to_Si30_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(P30) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si30_to_P30_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Si30 --> P30

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Si30_P30_meta, j_Si30_P30_rhoy, j_Si30_P30_temp, j_Si30_P30_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Si30_to_P30_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si30_to_P30_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Si30) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P31_to_Si31_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // P31 --> Si31

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_P31_Si31_meta, j_P31_Si31_rhoy, j_P31_Si31_temp, j_P31_Si31_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_P31_to_Si31_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P31_to_Si31_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(P31) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si31_to_P31_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Si31 --> P31

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Si31_P31_meta, j_Si31_P31_rhoy, j_Si31_P31_temp, j_Si31_P31_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Si31_to_P31_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si31_to_P31_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Si31) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P32_to_S32_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // P32 --> S32

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_P32_S32_meta, j_P32_S32_rhoy, j_P32_S32_temp, j_P32_S32_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_P32_to_S32_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P32_to_S32_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(P32) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P32_to_Si32_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // P32 --> Si32

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_P32_Si32_meta, j_P32_Si32_rhoy, j_P32_Si32_temp, j_P32_Si32_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_P32_to_Si32_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P32_to_Si32_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(P32) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S32_to_P32_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // S32 --> P32

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_S32_P32_meta, j_S32_P32_rhoy, j_S32_P32_temp, j_S32_P32_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_S32_to_P32_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S32_to_P32_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(S32) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si32_to_P32_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Si32 --> P32

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Si32_P32_meta, j_Si32_P32_rhoy, j_Si32_P32_temp, j_Si32_P32_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Si32_to_P32_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si32_to_P32_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Si32) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl33_to_S33_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cl33 --> S33

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cl33_S33_meta, j_Cl33_S33_rhoy, j_Cl33_S33_temp, j_Cl33_S33_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cl33_to_S33_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl33_to_S33_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cl33) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P33_to_S33_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // P33 --> S33

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_P33_S33_meta, j_P33_S33_rhoy, j_P33_S33_temp, j_P33_S33_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_P33_to_S33_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P33_to_S33_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(P33) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S33_to_Cl33_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // S33 --> Cl33

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_S33_Cl33_meta, j_S33_Cl33_rhoy, j_S33_Cl33_temp, j_S33_Cl33_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_S33_to_Cl33_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S33_to_Cl33_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(S33) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S33_to_P33_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // S33 --> P33

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_S33_P33_meta, j_S33_P33_rhoy, j_S33_P33_temp, j_S33_P33_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_S33_to_P33_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S33_to_P33_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(S33) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl34_to_S34_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cl34 --> S34

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cl34_S34_meta, j_Cl34_S34_rhoy, j_Cl34_S34_temp, j_Cl34_S34_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cl34_to_S34_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl34_to_S34_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cl34) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S34_to_Cl34_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // S34 --> Cl34

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_S34_Cl34_meta, j_S34_Cl34_rhoy, j_S34_Cl34_temp, j_S34_Cl34_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_S34_to_Cl34_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S34_to_Cl34_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(S34) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl35_to_S35_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cl35 --> S35

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cl35_S35_meta, j_Cl35_S35_rhoy, j_Cl35_S35_temp, j_Cl35_S35_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cl35_to_S35_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl35_to_S35_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cl35) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S35_to_Cl35_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // S35 --> Cl35

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_S35_Cl35_meta, j_S35_Cl35_rhoy, j_S35_Cl35_temp, j_S35_Cl35_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_S35_to_Cl35_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S35_to_Cl35_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(S35) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar36_to_Cl36_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ar36 --> Cl36

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ar36_Cl36_meta, j_Ar36_Cl36_rhoy, j_Ar36_Cl36_temp, j_Ar36_Cl36_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ar36_to_Cl36_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar36_to_Cl36_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ar36) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl36_to_Ar36_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cl36 --> Ar36

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cl36_Ar36_meta, j_Cl36_Ar36_rhoy, j_Cl36_Ar36_temp, j_Cl36_Ar36_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cl36_to_Ar36_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl36_to_Ar36_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cl36) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl36_to_S36_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cl36 --> S36

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cl36_S36_meta, j_Cl36_S36_rhoy, j_Cl36_S36_temp, j_Cl36_S36_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cl36_to_S36_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl36_to_S36_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cl36) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S36_to_Cl36_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // S36 --> Cl36

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_S36_Cl36_meta, j_S36_Cl36_rhoy, j_S36_Cl36_temp, j_S36_Cl36_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_S36_to_Cl36_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S36_to_Cl36_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(S36) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar37_to_Cl37_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ar37 --> Cl37

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ar37_Cl37_meta, j_Ar37_Cl37_rhoy, j_Ar37_Cl37_temp, j_Ar37_Cl37_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ar37_to_Cl37_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar37_to_Cl37_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ar37) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar37_to_K37_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ar37 --> K37

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ar37_K37_meta, j_Ar37_K37_rhoy, j_Ar37_K37_temp, j_Ar37_K37_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ar37_to_K37_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar37_to_K37_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ar37) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl37_to_Ar37_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cl37 --> Ar37

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cl37_Ar37_meta, j_Cl37_Ar37_rhoy, j_Cl37_Ar37_temp, j_Cl37_Ar37_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cl37_to_Ar37_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl37_to_Ar37_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cl37) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K37_to_Ar37_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // K37 --> Ar37

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_K37_Ar37_meta, j_K37_Ar37_rhoy, j_K37_Ar37_temp, j_K37_Ar37_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_K37_to_Ar37_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K37_to_Ar37_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(K37) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar38_to_K38_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ar38 --> K38

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ar38_K38_meta, j_Ar38_K38_rhoy, j_Ar38_K38_temp, j_Ar38_K38_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ar38_to_K38_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar38_to_K38_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ar38) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K38_to_Ar38_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // K38 --> Ar38

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_K38_Ar38_meta, j_K38_Ar38_rhoy, j_K38_Ar38_temp, j_K38_Ar38_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_K38_to_Ar38_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K38_to_Ar38_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(K38) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar39_to_K39_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ar39 --> K39

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ar39_K39_meta, j_Ar39_K39_rhoy, j_Ar39_K39_temp, j_Ar39_K39_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ar39_to_K39_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar39_to_K39_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ar39) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K39_to_Ar39_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // K39 --> Ar39

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_K39_Ar39_meta, j_K39_Ar39_rhoy, j_K39_Ar39_temp, j_K39_Ar39_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_K39_to_Ar39_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K39_to_Ar39_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(K39) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar40_to_K40_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ar40 --> K40

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ar40_K40_meta, j_Ar40_K40_rhoy, j_Ar40_K40_temp, j_Ar40_K40_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ar40_to_K40_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar40_to_K40_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ar40) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca40_to_K40_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ca40 --> K40

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ca40_K40_meta, j_Ca40_K40_rhoy, j_Ca40_K40_temp, j_Ca40_K40_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ca40_to_K40_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca40_to_K40_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ca40) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K40_to_Ar40_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // K40 --> Ar40

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_K40_Ar40_meta, j_K40_Ar40_rhoy, j_K40_Ar40_temp, j_K40_Ar40_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_K40_to_Ar40_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K40_to_Ar40_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(K40) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K40_to_Ca40_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // K40 --> Ca40

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_K40_Ca40_meta, j_K40_Ca40_rhoy, j_K40_Ca40_temp, j_K40_Ca40_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_K40_to_Ca40_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K40_to_Ca40_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(K40) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca41_to_K41_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ca41 --> K41

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ca41_K41_meta, j_Ca41_K41_rhoy, j_Ca41_K41_temp, j_Ca41_K41_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ca41_to_K41_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca41_to_K41_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ca41) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K41_to_Ca41_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // K41 --> Ca41

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_K41_Ca41_meta, j_K41_Ca41_rhoy, j_K41_Ca41_temp, j_K41_Ca41_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_K41_to_Ca41_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K41_to_Ca41_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(K41) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca43_to_Sc43_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ca43 --> Sc43

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ca43_Sc43_meta, j_Ca43_Sc43_rhoy, j_Ca43_Sc43_temp, j_Ca43_Sc43_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ca43_to_Sc43_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca43_to_Sc43_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ca43) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc43_to_Ca43_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Sc43 --> Ca43

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Sc43_Ca43_meta, j_Sc43_Ca43_rhoy, j_Sc43_Ca43_temp, j_Sc43_Ca43_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Sc43_to_Ca43_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc43_to_Ca43_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Sc43) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca44_to_Sc44_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ca44 --> Sc44

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ca44_Sc44_meta, j_Ca44_Sc44_rhoy, j_Ca44_Sc44_temp, j_Ca44_Sc44_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ca44_to_Sc44_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca44_to_Sc44_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ca44) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc44_to_Ca44_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Sc44 --> Ca44

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Sc44_Ca44_meta, j_Sc44_Ca44_rhoy, j_Sc44_Ca44_temp, j_Sc44_Ca44_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Sc44_to_Ca44_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc44_to_Ca44_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Sc44) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc44_to_Ti44_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Sc44 --> Ti44

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Sc44_Ti44_meta, j_Sc44_Ti44_rhoy, j_Sc44_Ti44_temp, j_Sc44_Ti44_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Sc44_to_Ti44_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc44_to_Ti44_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Sc44) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti44_to_Sc44_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ti44 --> Sc44

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ti44_Sc44_meta, j_Ti44_Sc44_rhoy, j_Ti44_Sc44_temp, j_Ti44_Sc44_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ti44_to_Sc44_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti44_to_Sc44_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ti44) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co53_to_Fe53_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Co53 --> Fe53

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Co53_Fe53_meta, j_Co53_Fe53_rhoy, j_Co53_Fe53_temp, j_Co53_Fe53_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Co53_to_Fe53_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co53_to_Fe53_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Co53) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe53_to_Co53_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Fe53 --> Co53

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Fe53_Co53_meta, j_Fe53_Co53_rhoy, j_Fe53_Co53_temp, j_Fe53_Co53_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Fe53_to_Co53_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe53_to_Co53_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Fe53) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu57_to_Ni57_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu57 --> Ni57

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu57_Ni57_meta, j_Cu57_Ni57_rhoy, j_Cu57_Ni57_temp, j_Cu57_Ni57_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu57_to_Ni57_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu57_to_Ni57_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu57) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni57_to_Cu57_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ni57 --> Cu57

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ni57_Cu57_meta, j_Ni57_Cu57_rhoy, j_Ni57_Cu57_temp, j_Ni57_Cu57_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ni57_to_Cu57_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni57_to_Cu57_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ni57) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca45_to_Sc45_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ca45 --> Sc45

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ca45_Sc45_meta, j_Ca45_Sc45_rhoy, j_Ca45_Sc45_temp, j_Ca45_Sc45_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ca45_to_Sc45_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca45_to_Sc45_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ca45) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc45_to_Ca45_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Sc45 --> Ca45

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Sc45_Ca45_meta, j_Sc45_Ca45_rhoy, j_Sc45_Ca45_temp, j_Sc45_Ca45_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Sc45_to_Ca45_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc45_to_Ca45_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Sc45) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc45_to_Ti45_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Sc45 --> Ti45

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Sc45_Ti45_meta, j_Sc45_Ti45_rhoy, j_Sc45_Ti45_temp, j_Sc45_Ti45_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Sc45_to_Ti45_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc45_to_Ti45_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Sc45) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti45_to_Sc45_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ti45 --> Sc45

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ti45_Sc45_meta, j_Ti45_Sc45_rhoy, j_Ti45_Sc45_temp, j_Ti45_Sc45_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ti45_to_Sc45_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti45_to_Sc45_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ti45) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca46_to_Sc46_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ca46 --> Sc46

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ca46_Sc46_meta, j_Ca46_Sc46_rhoy, j_Ca46_Sc46_temp, j_Ca46_Sc46_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ca46_to_Sc46_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca46_to_Sc46_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ca46) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc46_to_Ca46_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Sc46 --> Ca46

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Sc46_Ca46_meta, j_Sc46_Ca46_rhoy, j_Sc46_Ca46_temp, j_Sc46_Ca46_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Sc46_to_Ca46_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc46_to_Ca46_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Sc46) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc46_to_Ti46_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Sc46 --> Ti46

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Sc46_Ti46_meta, j_Sc46_Ti46_rhoy, j_Sc46_Ti46_temp, j_Sc46_Ti46_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Sc46_to_Ti46_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc46_to_Ti46_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Sc46) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti46_to_Sc46_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ti46 --> Sc46

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ti46_Sc46_meta, j_Ti46_Sc46_rhoy, j_Ti46_Sc46_temp, j_Ti46_Sc46_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ti46_to_Sc46_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti46_to_Sc46_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ti46) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti46_to_V46_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ti46 --> V46

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ti46_V46_meta, j_Ti46_V46_rhoy, j_Ti46_V46_temp, j_Ti46_V46_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ti46_to_V46_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti46_to_V46_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ti46) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V46_to_Ti46_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // V46 --> Ti46

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_V46_Ti46_meta, j_V46_Ti46_rhoy, j_V46_Ti46_temp, j_V46_Ti46_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_V46_to_Ti46_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V46_to_Ti46_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(V46) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca47_to_Sc47_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ca47 --> Sc47

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ca47_Sc47_meta, j_Ca47_Sc47_rhoy, j_Ca47_Sc47_temp, j_Ca47_Sc47_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ca47_to_Sc47_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca47_to_Sc47_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ca47) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc47_to_Ca47_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Sc47 --> Ca47

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Sc47_Ca47_meta, j_Sc47_Ca47_rhoy, j_Sc47_Ca47_temp, j_Sc47_Ca47_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Sc47_to_Ca47_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc47_to_Ca47_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Sc47) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc47_to_Ti47_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Sc47 --> Ti47

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Sc47_Ti47_meta, j_Sc47_Ti47_rhoy, j_Sc47_Ti47_temp, j_Sc47_Ti47_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Sc47_to_Ti47_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc47_to_Ti47_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Sc47) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti47_to_Sc47_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ti47 --> Sc47

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ti47_Sc47_meta, j_Ti47_Sc47_rhoy, j_Ti47_Sc47_temp, j_Ti47_Sc47_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ti47_to_Sc47_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti47_to_Sc47_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ti47) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti47_to_V47_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ti47 --> V47

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ti47_V47_meta, j_Ti47_V47_rhoy, j_Ti47_V47_temp, j_Ti47_V47_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ti47_to_V47_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti47_to_V47_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ti47) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V47_to_Ti47_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // V47 --> Ti47

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_V47_Ti47_meta, j_V47_Ti47_rhoy, j_V47_Ti47_temp, j_V47_Ti47_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_V47_to_Ti47_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V47_to_Ti47_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(V47) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca48_to_Sc48_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ca48 --> Sc48

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ca48_Sc48_meta, j_Ca48_Sc48_rhoy, j_Ca48_Sc48_temp, j_Ca48_Sc48_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ca48_to_Sc48_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca48_to_Sc48_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ca48) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr48_to_V48_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cr48 --> V48

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cr48_V48_meta, j_Cr48_V48_rhoy, j_Cr48_V48_temp, j_Cr48_V48_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cr48_to_V48_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr48_to_V48_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cr48) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc48_to_Ca48_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Sc48 --> Ca48

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Sc48_Ca48_meta, j_Sc48_Ca48_rhoy, j_Sc48_Ca48_temp, j_Sc48_Ca48_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Sc48_to_Ca48_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc48_to_Ca48_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Sc48) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc48_to_Ti48_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Sc48 --> Ti48

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Sc48_Ti48_meta, j_Sc48_Ti48_rhoy, j_Sc48_Ti48_temp, j_Sc48_Ti48_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Sc48_to_Ti48_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc48_to_Ti48_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Sc48) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti48_to_Sc48_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ti48 --> Sc48

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ti48_Sc48_meta, j_Ti48_Sc48_rhoy, j_Ti48_Sc48_temp, j_Ti48_Sc48_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ti48_to_Sc48_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti48_to_Sc48_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ti48) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti48_to_V48_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ti48 --> V48

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ti48_V48_meta, j_Ti48_V48_rhoy, j_Ti48_V48_temp, j_Ti48_V48_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ti48_to_V48_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti48_to_V48_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ti48) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V48_to_Cr48_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // V48 --> Cr48

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_V48_Cr48_meta, j_V48_Cr48_rhoy, j_V48_Cr48_temp, j_V48_Cr48_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_V48_to_Cr48_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V48_to_Cr48_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(V48) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V48_to_Ti48_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // V48 --> Ti48

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_V48_Ti48_meta, j_V48_Ti48_rhoy, j_V48_Ti48_temp, j_V48_Ti48_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_V48_to_Ti48_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V48_to_Ti48_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(V48) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr49_to_V49_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cr49 --> V49

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cr49_V49_meta, j_Cr49_V49_rhoy, j_Cr49_V49_temp, j_Cr49_V49_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cr49_to_V49_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr49_to_V49_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cr49) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc49_to_Ti49_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Sc49 --> Ti49

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Sc49_Ti49_meta, j_Sc49_Ti49_rhoy, j_Sc49_Ti49_temp, j_Sc49_Ti49_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Sc49_to_Ti49_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc49_to_Ti49_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Sc49) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti49_to_Sc49_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ti49 --> Sc49

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ti49_Sc49_meta, j_Ti49_Sc49_rhoy, j_Ti49_Sc49_temp, j_Ti49_Sc49_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ti49_to_Sc49_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti49_to_Sc49_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ti49) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti49_to_V49_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ti49 --> V49

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ti49_V49_meta, j_Ti49_V49_rhoy, j_Ti49_V49_temp, j_Ti49_V49_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ti49_to_V49_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti49_to_V49_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ti49) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V49_to_Cr49_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // V49 --> Cr49

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_V49_Cr49_meta, j_V49_Cr49_rhoy, j_V49_Cr49_temp, j_V49_Cr49_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_V49_to_Cr49_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V49_to_Cr49_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(V49) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V49_to_Ti49_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // V49 --> Ti49

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_V49_Ti49_meta, j_V49_Ti49_rhoy, j_V49_Ti49_temp, j_V49_Ti49_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_V49_to_Ti49_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V49_to_Ti49_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(V49) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr50_to_Mn50_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cr50 --> Mn50

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cr50_Mn50_meta, j_Cr50_Mn50_rhoy, j_Cr50_Mn50_temp, j_Cr50_Mn50_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cr50_to_Mn50_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr50_to_Mn50_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cr50) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr50_to_V50_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cr50 --> V50

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cr50_V50_meta, j_Cr50_V50_rhoy, j_Cr50_V50_temp, j_Cr50_V50_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cr50_to_V50_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr50_to_V50_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cr50) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn50_to_Cr50_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Mn50 --> Cr50

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Mn50_Cr50_meta, j_Mn50_Cr50_rhoy, j_Mn50_Cr50_temp, j_Mn50_Cr50_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Mn50_to_Cr50_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn50_to_Cr50_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Mn50) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti50_to_V50_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ti50 --> V50

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ti50_V50_meta, j_Ti50_V50_rhoy, j_Ti50_V50_temp, j_Ti50_V50_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ti50_to_V50_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti50_to_V50_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ti50) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V50_to_Cr50_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // V50 --> Cr50

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_V50_Cr50_meta, j_V50_Cr50_rhoy, j_V50_Cr50_temp, j_V50_Cr50_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_V50_to_Cr50_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V50_to_Cr50_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(V50) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V50_to_Ti50_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // V50 --> Ti50

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_V50_Ti50_meta, j_V50_Ti50_rhoy, j_V50_Ti50_temp, j_V50_Ti50_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_V50_to_Ti50_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V50_to_Ti50_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(V50) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr51_to_Mn51_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cr51 --> Mn51

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cr51_Mn51_meta, j_Cr51_Mn51_rhoy, j_Cr51_Mn51_temp, j_Cr51_Mn51_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cr51_to_Mn51_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr51_to_Mn51_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cr51) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr51_to_V51_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cr51 --> V51

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cr51_V51_meta, j_Cr51_V51_rhoy, j_Cr51_V51_temp, j_Cr51_V51_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cr51_to_V51_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr51_to_V51_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cr51) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn51_to_Cr51_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Mn51 --> Cr51

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Mn51_Cr51_meta, j_Mn51_Cr51_rhoy, j_Mn51_Cr51_temp, j_Mn51_Cr51_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Mn51_to_Cr51_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn51_to_Cr51_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Mn51) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti51_to_V51_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ti51 --> V51

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ti51_V51_meta, j_Ti51_V51_rhoy, j_Ti51_V51_temp, j_Ti51_V51_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ti51_to_V51_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti51_to_V51_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ti51) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V51_to_Cr51_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // V51 --> Cr51

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_V51_Cr51_meta, j_V51_Cr51_rhoy, j_V51_Cr51_temp, j_V51_Cr51_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_V51_to_Cr51_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V51_to_Cr51_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(V51) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V51_to_Ti51_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // V51 --> Ti51

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_V51_Ti51_meta, j_V51_Ti51_rhoy, j_V51_Ti51_temp, j_V51_Ti51_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_V51_to_Ti51_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V51_to_Ti51_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(V51) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr52_to_Mn52_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cr52 --> Mn52

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cr52_Mn52_meta, j_Cr52_Mn52_rhoy, j_Cr52_Mn52_temp, j_Cr52_Mn52_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cr52_to_Mn52_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr52_to_Mn52_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cr52) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr52_to_V52_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cr52 --> V52

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cr52_V52_meta, j_Cr52_V52_rhoy, j_Cr52_V52_temp, j_Cr52_V52_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cr52_to_V52_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr52_to_V52_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cr52) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe52_to_Mn52_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Fe52 --> Mn52

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Fe52_Mn52_meta, j_Fe52_Mn52_rhoy, j_Fe52_Mn52_temp, j_Fe52_Mn52_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Fe52_to_Mn52_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe52_to_Mn52_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Fe52) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn52_to_Cr52_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Mn52 --> Cr52

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Mn52_Cr52_meta, j_Mn52_Cr52_rhoy, j_Mn52_Cr52_temp, j_Mn52_Cr52_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Mn52_to_Cr52_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn52_to_Cr52_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Mn52) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn52_to_Fe52_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Mn52 --> Fe52

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Mn52_Fe52_meta, j_Mn52_Fe52_rhoy, j_Mn52_Fe52_temp, j_Mn52_Fe52_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Mn52_to_Fe52_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn52_to_Fe52_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Mn52) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V52_to_Cr52_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // V52 --> Cr52

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_V52_Cr52_meta, j_V52_Cr52_rhoy, j_V52_Cr52_temp, j_V52_Cr52_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_V52_to_Cr52_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V52_to_Cr52_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(V52) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr53_to_Mn53_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cr53 --> Mn53

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cr53_Mn53_meta, j_Cr53_Mn53_rhoy, j_Cr53_Mn53_temp, j_Cr53_Mn53_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cr53_to_Mn53_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr53_to_Mn53_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cr53) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe53_to_Mn53_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Fe53 --> Mn53

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Fe53_Mn53_meta, j_Fe53_Mn53_rhoy, j_Fe53_Mn53_temp, j_Fe53_Mn53_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Fe53_to_Mn53_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe53_to_Mn53_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Fe53) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn53_to_Cr53_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Mn53 --> Cr53

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Mn53_Cr53_meta, j_Mn53_Cr53_rhoy, j_Mn53_Cr53_temp, j_Mn53_Cr53_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Mn53_to_Cr53_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn53_to_Cr53_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Mn53) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn53_to_Fe53_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Mn53 --> Fe53

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Mn53_Fe53_meta, j_Mn53_Fe53_rhoy, j_Mn53_Fe53_temp, j_Mn53_Fe53_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Mn53_to_Fe53_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn53_to_Fe53_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Mn53) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co54_to_Fe54_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Co54 --> Fe54

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Co54_Fe54_meta, j_Co54_Fe54_rhoy, j_Co54_Fe54_temp, j_Co54_Fe54_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Co54_to_Fe54_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co54_to_Fe54_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Co54) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr54_to_Mn54_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cr54 --> Mn54

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cr54_Mn54_meta, j_Cr54_Mn54_rhoy, j_Cr54_Mn54_temp, j_Cr54_Mn54_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cr54_to_Mn54_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr54_to_Mn54_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cr54) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe54_to_Co54_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Fe54 --> Co54

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Fe54_Co54_meta, j_Fe54_Co54_rhoy, j_Fe54_Co54_temp, j_Fe54_Co54_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Fe54_to_Co54_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe54_to_Co54_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Fe54) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe54_to_Mn54_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Fe54 --> Mn54

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Fe54_Mn54_meta, j_Fe54_Mn54_rhoy, j_Fe54_Mn54_temp, j_Fe54_Mn54_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Fe54_to_Mn54_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe54_to_Mn54_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Fe54) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn54_to_Cr54_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Mn54 --> Cr54

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Mn54_Cr54_meta, j_Mn54_Cr54_rhoy, j_Mn54_Cr54_temp, j_Mn54_Cr54_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Mn54_to_Cr54_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn54_to_Cr54_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Mn54) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn54_to_Fe54_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Mn54 --> Fe54

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Mn54_Fe54_meta, j_Mn54_Fe54_rhoy, j_Mn54_Fe54_temp, j_Mn54_Fe54_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Mn54_to_Fe54_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn54_to_Fe54_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Mn54) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co55_to_Fe55_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Co55 --> Fe55

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Co55_Fe55_meta, j_Co55_Fe55_rhoy, j_Co55_Fe55_temp, j_Co55_Fe55_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Co55_to_Fe55_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co55_to_Fe55_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Co55) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe55_to_Co55_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Fe55 --> Co55

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Fe55_Co55_meta, j_Fe55_Co55_rhoy, j_Fe55_Co55_temp, j_Fe55_Co55_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Fe55_to_Co55_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe55_to_Co55_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Fe55) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe55_to_Mn55_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Fe55 --> Mn55

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Fe55_Mn55_meta, j_Fe55_Mn55_rhoy, j_Fe55_Mn55_temp, j_Fe55_Mn55_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Fe55_to_Mn55_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe55_to_Mn55_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Fe55) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn55_to_Fe55_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Mn55 --> Fe55

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Mn55_Fe55_meta, j_Mn55_Fe55_rhoy, j_Mn55_Fe55_temp, j_Mn55_Fe55_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Mn55_to_Fe55_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn55_to_Fe55_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Mn55) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co56_to_Fe56_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Co56 --> Fe56

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Co56_Fe56_meta, j_Co56_Fe56_rhoy, j_Co56_Fe56_temp, j_Co56_Fe56_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Co56_to_Fe56_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co56_to_Fe56_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Co56) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co56_to_Ni56_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Co56 --> Ni56

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Co56_Ni56_meta, j_Co56_Ni56_rhoy, j_Co56_Ni56_temp, j_Co56_Ni56_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Co56_to_Ni56_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co56_to_Ni56_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Co56) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe56_to_Co56_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Fe56 --> Co56

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Fe56_Co56_meta, j_Fe56_Co56_rhoy, j_Fe56_Co56_temp, j_Fe56_Co56_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Fe56_to_Co56_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe56_to_Co56_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Fe56) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni56_to_Co56_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ni56 --> Co56

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ni56_Co56_meta, j_Ni56_Co56_rhoy, j_Ni56_Co56_temp, j_Ni56_Co56_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ni56_to_Co56_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni56_to_Co56_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ni56) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co57_to_Fe57_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Co57 --> Fe57

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Co57_Fe57_meta, j_Co57_Fe57_rhoy, j_Co57_Fe57_temp, j_Co57_Fe57_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Co57_to_Fe57_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co57_to_Fe57_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Co57) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co57_to_Ni57_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Co57 --> Ni57

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Co57_Ni57_meta, j_Co57_Ni57_rhoy, j_Co57_Ni57_temp, j_Co57_Ni57_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Co57_to_Ni57_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co57_to_Ni57_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Co57) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe57_to_Co57_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Fe57 --> Co57

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Fe57_Co57_meta, j_Fe57_Co57_rhoy, j_Fe57_Co57_temp, j_Fe57_Co57_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Fe57_to_Co57_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe57_to_Co57_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Fe57) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni57_to_Co57_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ni57 --> Co57

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ni57_Co57_meta, j_Ni57_Co57_rhoy, j_Ni57_Co57_temp, j_Ni57_Co57_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ni57_to_Co57_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni57_to_Co57_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ni57) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co58_to_Fe58_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Co58 --> Fe58

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Co58_Fe58_meta, j_Co58_Fe58_rhoy, j_Co58_Fe58_temp, j_Co58_Fe58_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Co58_to_Fe58_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co58_to_Fe58_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Co58) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co58_to_Ni58_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Co58 --> Ni58

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Co58_Ni58_meta, j_Co58_Ni58_rhoy, j_Co58_Ni58_temp, j_Co58_Ni58_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Co58_to_Ni58_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co58_to_Ni58_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Co58) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu58_to_Ni58_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu58 --> Ni58

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu58_Ni58_meta, j_Cu58_Ni58_rhoy, j_Cu58_Ni58_temp, j_Cu58_Ni58_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu58_to_Ni58_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu58_to_Ni58_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu58) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe58_to_Co58_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Fe58 --> Co58

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Fe58_Co58_meta, j_Fe58_Co58_rhoy, j_Fe58_Co58_temp, j_Fe58_Co58_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Fe58_to_Co58_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe58_to_Co58_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Fe58) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni58_to_Co58_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ni58 --> Co58

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ni58_Co58_meta, j_Ni58_Co58_rhoy, j_Ni58_Co58_temp, j_Ni58_Co58_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ni58_to_Co58_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni58_to_Co58_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ni58) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni58_to_Cu58_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ni58 --> Cu58

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ni58_Cu58_meta, j_Ni58_Cu58_rhoy, j_Ni58_Cu58_temp, j_Ni58_Cu58_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ni58_to_Cu58_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni58_to_Cu58_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ni58) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co59_to_Ni59_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Co59 --> Ni59

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Co59_Ni59_meta, j_Co59_Ni59_rhoy, j_Co59_Ni59_temp, j_Co59_Ni59_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Co59_to_Ni59_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co59_to_Ni59_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Co59) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu59_to_Ni59_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu59 --> Ni59

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu59_Ni59_meta, j_Cu59_Ni59_rhoy, j_Cu59_Ni59_temp, j_Cu59_Ni59_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu59_to_Ni59_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu59_to_Ni59_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu59) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni59_to_Co59_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ni59 --> Co59

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ni59_Co59_meta, j_Ni59_Co59_rhoy, j_Ni59_Co59_temp, j_Ni59_Co59_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ni59_to_Co59_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni59_to_Co59_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ni59) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni59_to_Cu59_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ni59 --> Cu59

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ni59_Cu59_meta, j_Ni59_Cu59_rhoy, j_Ni59_Cu59_temp, j_Ni59_Cu59_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ni59_to_Cu59_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni59_to_Cu59_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ni59) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu60_to_Ni60_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu60 --> Ni60

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu60_Ni60_meta, j_Cu60_Ni60_rhoy, j_Cu60_Ni60_temp, j_Cu60_Ni60_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu60_to_Ni60_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu60_to_Ni60_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu60) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu60_to_Zn60_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu60 --> Zn60

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu60_Zn60_meta, j_Cu60_Zn60_rhoy, j_Cu60_Zn60_temp, j_Cu60_Zn60_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu60_to_Zn60_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu60_to_Zn60_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu60) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni60_to_Cu60_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ni60 --> Cu60

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ni60_Cu60_meta, j_Ni60_Cu60_rhoy, j_Ni60_Cu60_temp, j_Ni60_Cu60_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ni60_to_Cu60_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni60_to_Cu60_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ni60) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Zn60_to_Cu60_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Zn60 --> Cu60

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Zn60_Cu60_meta, j_Zn60_Cu60_rhoy, j_Zn60_Cu60_temp, j_Zn60_Cu60_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Zn60_to_Cu60_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn60_to_Cu60_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Zn60) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu61_to_Ni61_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu61 --> Ni61

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu61_Ni61_meta, j_Cu61_Ni61_rhoy, j_Cu61_Ni61_temp, j_Cu61_Ni61_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu61_to_Ni61_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu61_to_Ni61_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu61) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu61_to_Zn61_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu61 --> Zn61

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu61_Zn61_meta, j_Cu61_Zn61_rhoy, j_Cu61_Zn61_temp, j_Cu61_Zn61_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu61_to_Zn61_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu61_to_Zn61_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu61) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni61_to_Cu61_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ni61 --> Cu61

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ni61_Cu61_meta, j_Ni61_Cu61_rhoy, j_Ni61_Cu61_temp, j_Ni61_Cu61_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ni61_to_Cu61_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni61_to_Cu61_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ni61) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Zn61_to_Cu61_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Zn61 --> Cu61

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Zn61_Cu61_meta, j_Zn61_Cu61_rhoy, j_Zn61_Cu61_temp, j_Zn61_Cu61_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Zn61_to_Cu61_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn61_to_Cu61_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Zn61) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu62_to_Ni62_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu62 --> Ni62

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu62_Ni62_meta, j_Cu62_Ni62_rhoy, j_Cu62_Ni62_temp, j_Cu62_Ni62_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu62_to_Ni62_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu62_to_Ni62_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu62) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu62_to_Zn62_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu62 --> Zn62

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu62_Zn62_meta, j_Cu62_Zn62_rhoy, j_Cu62_Zn62_temp, j_Cu62_Zn62_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu62_to_Zn62_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu62_to_Zn62_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu62) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ga62_to_Zn62_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ga62 --> Zn62

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ga62_Zn62_meta, j_Ga62_Zn62_rhoy, j_Ga62_Zn62_temp, j_Ga62_Zn62_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ga62_to_Zn62_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ga62_to_Zn62_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ga62) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni62_to_Cu62_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ni62 --> Cu62

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ni62_Cu62_meta, j_Ni62_Cu62_rhoy, j_Ni62_Cu62_temp, j_Ni62_Cu62_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ni62_to_Cu62_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni62_to_Cu62_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ni62) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Zn62_to_Cu62_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Zn62 --> Cu62

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Zn62_Cu62_meta, j_Zn62_Cu62_rhoy, j_Zn62_Cu62_temp, j_Zn62_Cu62_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Zn62_to_Cu62_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn62_to_Cu62_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Zn62) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Zn62_to_Ga62_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Zn62 --> Ga62

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Zn62_Ga62_meta, j_Zn62_Ga62_rhoy, j_Zn62_Ga62_temp, j_Zn62_Ga62_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Zn62_to_Ga62_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn62_to_Ga62_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Zn62) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu63_to_Ni63_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu63 --> Ni63

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu63_Ni63_meta, j_Cu63_Ni63_rhoy, j_Cu63_Ni63_temp, j_Cu63_Ni63_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu63_to_Ni63_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu63_to_Ni63_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu63) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu63_to_Zn63_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu63 --> Zn63

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu63_Zn63_meta, j_Cu63_Zn63_rhoy, j_Cu63_Zn63_temp, j_Cu63_Zn63_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu63_to_Zn63_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu63_to_Zn63_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu63) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ga63_to_Zn63_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ga63 --> Zn63

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ga63_Zn63_meta, j_Ga63_Zn63_rhoy, j_Ga63_Zn63_temp, j_Ga63_Zn63_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ga63_to_Zn63_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ga63_to_Zn63_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ga63) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni63_to_Cu63_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ni63 --> Cu63

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ni63_Cu63_meta, j_Ni63_Cu63_rhoy, j_Ni63_Cu63_temp, j_Ni63_Cu63_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ni63_to_Cu63_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni63_to_Cu63_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ni63) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Zn63_to_Cu63_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Zn63 --> Cu63

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Zn63_Cu63_meta, j_Zn63_Cu63_rhoy, j_Zn63_Cu63_temp, j_Zn63_Cu63_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Zn63_to_Cu63_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn63_to_Cu63_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Zn63) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Zn63_to_Ga63_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Zn63 --> Ga63

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Zn63_Ga63_meta, j_Zn63_Ga63_rhoy, j_Zn63_Ga63_temp, j_Zn63_Ga63_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Zn63_to_Ga63_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn63_to_Ga63_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Zn63) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu64_to_Ni64_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu64 --> Ni64

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu64_Ni64_meta, j_Cu64_Ni64_rhoy, j_Cu64_Ni64_temp, j_Cu64_Ni64_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu64_to_Ni64_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu64_to_Ni64_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu64) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu64_to_Zn64_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu64 --> Zn64

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu64_Zn64_meta, j_Cu64_Zn64_rhoy, j_Cu64_Zn64_temp, j_Cu64_Zn64_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu64_to_Zn64_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu64_to_Zn64_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu64) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ga64_to_Ge64_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ga64 --> Ge64

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ga64_Ge64_meta, j_Ga64_Ge64_rhoy, j_Ga64_Ge64_temp, j_Ga64_Ge64_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ga64_to_Ge64_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ga64_to_Ge64_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ga64) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ga64_to_Zn64_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ga64 --> Zn64

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ga64_Zn64_meta, j_Ga64_Zn64_rhoy, j_Ga64_Zn64_temp, j_Ga64_Zn64_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ga64_to_Zn64_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ga64_to_Zn64_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ga64) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ge64_to_Ga64_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ge64 --> Ga64

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ge64_Ga64_meta, j_Ge64_Ga64_rhoy, j_Ge64_Ga64_temp, j_Ge64_Ga64_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ge64_to_Ga64_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ge64_to_Ga64_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ge64) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni64_to_Cu64_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ni64 --> Cu64

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ni64_Cu64_meta, j_Ni64_Cu64_rhoy, j_Ni64_Cu64_temp, j_Ni64_Cu64_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ni64_to_Cu64_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni64_to_Cu64_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ni64) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Zn64_to_Cu64_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Zn64 --> Cu64

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Zn64_Cu64_meta, j_Zn64_Cu64_rhoy, j_Zn64_Cu64_temp, j_Zn64_Cu64_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Zn64_to_Cu64_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn64_to_Cu64_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Zn64) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Zn64_to_Ga64_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Zn64 --> Ga64

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Zn64_Ga64_meta, j_Zn64_Ga64_rhoy, j_Zn64_Ga64_temp, j_Zn64_Ga64_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Zn64_to_Ga64_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn64_to_Ga64_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Zn64) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu65_to_Zn65_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu65 --> Zn65

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu65_Zn65_meta, j_Cu65_Zn65_rhoy, j_Cu65_Zn65_temp, j_Cu65_Zn65_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu65_to_Zn65_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu65_to_Zn65_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu65) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Zn65_to_Cu65_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Zn65 --> Cu65

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Zn65_Cu65_meta, j_Zn65_Cu65_rhoy, j_Zn65_Cu65_temp, j_Zn65_Cu65_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Zn65_to_Cu65_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn65_to_Cu65_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(Zn65) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_n_to_p_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // n --> p

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_n_p_meta, j_n_p_rhoy, j_n_p_temp, j_n_p_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_n_to_p_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_n_to_p_beta_neg_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(N) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_p_to_n_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // p --> n

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_p_n_meta, j_p_n_rhoy, j_p_n_temp, j_p_n_data,
                                            log_rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_p_to_n_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_to_n_electron_capture_weaktab) = table_values.drate_dT;
        }

        rate_eval.enuc_weak += C::n_A * Y(H1) * (table_values.enu + table_values.gamma);
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
            // Na21_to_Ne21_electron_capture_weaktab
            rate_Na21_to_Ne21_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ne21_to_Na21_beta_neg_weaktab
            rate_Ne21_to_Na21_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Na22_to_Ne22_electron_capture_weaktab
            rate_Na22_to_Ne22_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ne22_to_Na22_beta_neg_weaktab
            rate_Ne22_to_Na22_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Mg23_to_Na23_electron_capture_weaktab
            rate_Mg23_to_Na23_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Na23_to_Mg23_beta_neg_weaktab
            rate_Na23_to_Mg23_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Al25_to_Mg25_electron_capture_weaktab
            rate_Al25_to_Mg25_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Mg25_to_Al25_beta_neg_weaktab
            rate_Mg25_to_Al25_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Al26_to_Mg26_electron_capture_weaktab
            rate_Al26_to_Mg26_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Mg26_to_Al26_beta_neg_weaktab
            rate_Mg26_to_Al26_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // P29_to_Si29_electron_capture_weaktab
            rate_P29_to_Si29_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Si29_to_P29_beta_neg_weaktab
            rate_Si29_to_P29_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // P30_to_Si30_electron_capture_weaktab
            rate_P30_to_Si30_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Si30_to_P30_beta_neg_weaktab
            rate_Si30_to_P30_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // P31_to_Si31_electron_capture_weaktab
            rate_P31_to_Si31_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Si31_to_P31_beta_neg_weaktab
            rate_Si31_to_P31_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // P32_to_S32_beta_neg_weaktab
            rate_P32_to_S32_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // P32_to_Si32_electron_capture_weaktab
            rate_P32_to_Si32_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // S32_to_P32_electron_capture_weaktab
            rate_S32_to_P32_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Si32_to_P32_beta_neg_weaktab
            rate_Si32_to_P32_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cl33_to_S33_electron_capture_weaktab
            rate_Cl33_to_S33_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // P33_to_S33_beta_neg_weaktab
            rate_P33_to_S33_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // S33_to_Cl33_beta_neg_weaktab
            rate_S33_to_Cl33_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // S33_to_P33_electron_capture_weaktab
            rate_S33_to_P33_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cl34_to_S34_electron_capture_weaktab
            rate_Cl34_to_S34_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // S34_to_Cl34_beta_neg_weaktab
            rate_S34_to_Cl34_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cl35_to_S35_electron_capture_weaktab
            rate_Cl35_to_S35_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // S35_to_Cl35_beta_neg_weaktab
            rate_S35_to_Cl35_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ar36_to_Cl36_electron_capture_weaktab
            rate_Ar36_to_Cl36_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cl36_to_Ar36_beta_neg_weaktab
            rate_Cl36_to_Ar36_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cl36_to_S36_electron_capture_weaktab
            rate_Cl36_to_S36_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // S36_to_Cl36_beta_neg_weaktab
            rate_S36_to_Cl36_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ar37_to_Cl37_electron_capture_weaktab
            rate_Ar37_to_Cl37_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ar37_to_K37_beta_neg_weaktab
            rate_Ar37_to_K37_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cl37_to_Ar37_beta_neg_weaktab
            rate_Cl37_to_Ar37_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // K37_to_Ar37_electron_capture_weaktab
            rate_K37_to_Ar37_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ar38_to_K38_beta_neg_weaktab
            rate_Ar38_to_K38_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // K38_to_Ar38_electron_capture_weaktab
            rate_K38_to_Ar38_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ar39_to_K39_beta_neg_weaktab
            rate_Ar39_to_K39_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // K39_to_Ar39_electron_capture_weaktab
            rate_K39_to_Ar39_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ar40_to_K40_beta_neg_weaktab
            rate_Ar40_to_K40_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ca40_to_K40_electron_capture_weaktab
            rate_Ca40_to_K40_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // K40_to_Ar40_electron_capture_weaktab
            rate_K40_to_Ar40_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // K40_to_Ca40_beta_neg_weaktab
            rate_K40_to_Ca40_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ca41_to_K41_electron_capture_weaktab
            rate_Ca41_to_K41_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // K41_to_Ca41_beta_neg_weaktab
            rate_K41_to_Ca41_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ca43_to_Sc43_beta_neg_weaktab
            rate_Ca43_to_Sc43_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Sc43_to_Ca43_electron_capture_weaktab
            rate_Sc43_to_Ca43_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ca44_to_Sc44_beta_neg_weaktab
            rate_Ca44_to_Sc44_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Sc44_to_Ca44_electron_capture_weaktab
            rate_Sc44_to_Ca44_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Sc44_to_Ti44_beta_neg_weaktab
            rate_Sc44_to_Ti44_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ti44_to_Sc44_electron_capture_weaktab
            rate_Ti44_to_Sc44_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Co53_to_Fe53_electron_capture_weaktab
            rate_Co53_to_Fe53_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Fe53_to_Co53_beta_neg_weaktab
            rate_Fe53_to_Co53_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cu57_to_Ni57_electron_capture_weaktab
            rate_Cu57_to_Ni57_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ni57_to_Cu57_beta_neg_weaktab
            rate_Ni57_to_Cu57_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ca45_to_Sc45_beta_neg_weaktab
            rate_Ca45_to_Sc45_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Sc45_to_Ca45_electron_capture_weaktab
            rate_Sc45_to_Ca45_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Sc45_to_Ti45_beta_neg_weaktab
            rate_Sc45_to_Ti45_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ti45_to_Sc45_electron_capture_weaktab
            rate_Ti45_to_Sc45_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ca46_to_Sc46_beta_neg_weaktab
            rate_Ca46_to_Sc46_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Sc46_to_Ca46_electron_capture_weaktab
            rate_Sc46_to_Ca46_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Sc46_to_Ti46_beta_neg_weaktab
            rate_Sc46_to_Ti46_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ti46_to_Sc46_electron_capture_weaktab
            rate_Ti46_to_Sc46_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ti46_to_V46_beta_neg_weaktab
            rate_Ti46_to_V46_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // V46_to_Ti46_electron_capture_weaktab
            rate_V46_to_Ti46_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ca47_to_Sc47_beta_neg_weaktab
            rate_Ca47_to_Sc47_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Sc47_to_Ca47_electron_capture_weaktab
            rate_Sc47_to_Ca47_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Sc47_to_Ti47_beta_neg_weaktab
            rate_Sc47_to_Ti47_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ti47_to_Sc47_electron_capture_weaktab
            rate_Ti47_to_Sc47_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ti47_to_V47_beta_neg_weaktab
            rate_Ti47_to_V47_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // V47_to_Ti47_electron_capture_weaktab
            rate_V47_to_Ti47_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ca48_to_Sc48_beta_neg_weaktab
            rate_Ca48_to_Sc48_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cr48_to_V48_electron_capture_weaktab
            rate_Cr48_to_V48_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Sc48_to_Ca48_electron_capture_weaktab
            rate_Sc48_to_Ca48_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Sc48_to_Ti48_beta_neg_weaktab
            rate_Sc48_to_Ti48_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ti48_to_Sc48_electron_capture_weaktab
            rate_Ti48_to_Sc48_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ti48_to_V48_beta_neg_weaktab
            rate_Ti48_to_V48_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // V48_to_Cr48_beta_neg_weaktab
            rate_V48_to_Cr48_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // V48_to_Ti48_electron_capture_weaktab
            rate_V48_to_Ti48_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cr49_to_V49_electron_capture_weaktab
            rate_Cr49_to_V49_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Sc49_to_Ti49_beta_neg_weaktab
            rate_Sc49_to_Ti49_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ti49_to_Sc49_electron_capture_weaktab
            rate_Ti49_to_Sc49_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ti49_to_V49_beta_neg_weaktab
            rate_Ti49_to_V49_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // V49_to_Cr49_beta_neg_weaktab
            rate_V49_to_Cr49_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // V49_to_Ti49_electron_capture_weaktab
            rate_V49_to_Ti49_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cr50_to_Mn50_beta_neg_weaktab
            rate_Cr50_to_Mn50_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cr50_to_V50_electron_capture_weaktab
            rate_Cr50_to_V50_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Mn50_to_Cr50_electron_capture_weaktab
            rate_Mn50_to_Cr50_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ti50_to_V50_beta_neg_weaktab
            rate_Ti50_to_V50_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // V50_to_Cr50_beta_neg_weaktab
            rate_V50_to_Cr50_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // V50_to_Ti50_electron_capture_weaktab
            rate_V50_to_Ti50_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cr51_to_Mn51_beta_neg_weaktab
            rate_Cr51_to_Mn51_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cr51_to_V51_electron_capture_weaktab
            rate_Cr51_to_V51_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Mn51_to_Cr51_electron_capture_weaktab
            rate_Mn51_to_Cr51_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ti51_to_V51_beta_neg_weaktab
            rate_Ti51_to_V51_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // V51_to_Cr51_beta_neg_weaktab
            rate_V51_to_Cr51_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // V51_to_Ti51_electron_capture_weaktab
            rate_V51_to_Ti51_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cr52_to_Mn52_beta_neg_weaktab
            rate_Cr52_to_Mn52_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cr52_to_V52_electron_capture_weaktab
            rate_Cr52_to_V52_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Fe52_to_Mn52_electron_capture_weaktab
            rate_Fe52_to_Mn52_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Mn52_to_Cr52_electron_capture_weaktab
            rate_Mn52_to_Cr52_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Mn52_to_Fe52_beta_neg_weaktab
            rate_Mn52_to_Fe52_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // V52_to_Cr52_beta_neg_weaktab
            rate_V52_to_Cr52_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cr53_to_Mn53_beta_neg_weaktab
            rate_Cr53_to_Mn53_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Fe53_to_Mn53_electron_capture_weaktab
            rate_Fe53_to_Mn53_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Mn53_to_Cr53_electron_capture_weaktab
            rate_Mn53_to_Cr53_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Mn53_to_Fe53_beta_neg_weaktab
            rate_Mn53_to_Fe53_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Co54_to_Fe54_electron_capture_weaktab
            rate_Co54_to_Fe54_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cr54_to_Mn54_beta_neg_weaktab
            rate_Cr54_to_Mn54_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Fe54_to_Co54_beta_neg_weaktab
            rate_Fe54_to_Co54_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Fe54_to_Mn54_electron_capture_weaktab
            rate_Fe54_to_Mn54_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Mn54_to_Cr54_electron_capture_weaktab
            rate_Mn54_to_Cr54_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Mn54_to_Fe54_beta_neg_weaktab
            rate_Mn54_to_Fe54_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Co55_to_Fe55_electron_capture_weaktab
            rate_Co55_to_Fe55_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Fe55_to_Co55_beta_neg_weaktab
            rate_Fe55_to_Co55_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Fe55_to_Mn55_electron_capture_weaktab
            rate_Fe55_to_Mn55_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Mn55_to_Fe55_beta_neg_weaktab
            rate_Mn55_to_Fe55_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Co56_to_Fe56_electron_capture_weaktab
            rate_Co56_to_Fe56_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Co56_to_Ni56_beta_neg_weaktab
            rate_Co56_to_Ni56_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Fe56_to_Co56_beta_neg_weaktab
            rate_Fe56_to_Co56_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ni56_to_Co56_electron_capture_weaktab
            rate_Ni56_to_Co56_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Co57_to_Fe57_electron_capture_weaktab
            rate_Co57_to_Fe57_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Co57_to_Ni57_beta_neg_weaktab
            rate_Co57_to_Ni57_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Fe57_to_Co57_beta_neg_weaktab
            rate_Fe57_to_Co57_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ni57_to_Co57_electron_capture_weaktab
            rate_Ni57_to_Co57_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Co58_to_Fe58_electron_capture_weaktab
            rate_Co58_to_Fe58_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Co58_to_Ni58_beta_neg_weaktab
            rate_Co58_to_Ni58_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cu58_to_Ni58_electron_capture_weaktab
            rate_Cu58_to_Ni58_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Fe58_to_Co58_beta_neg_weaktab
            rate_Fe58_to_Co58_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ni58_to_Co58_electron_capture_weaktab
            rate_Ni58_to_Co58_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ni58_to_Cu58_beta_neg_weaktab
            rate_Ni58_to_Cu58_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Co59_to_Ni59_beta_neg_weaktab
            rate_Co59_to_Ni59_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cu59_to_Ni59_electron_capture_weaktab
            rate_Cu59_to_Ni59_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ni59_to_Co59_electron_capture_weaktab
            rate_Ni59_to_Co59_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ni59_to_Cu59_beta_neg_weaktab
            rate_Ni59_to_Cu59_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cu60_to_Ni60_electron_capture_weaktab
            rate_Cu60_to_Ni60_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cu60_to_Zn60_beta_neg_weaktab
            rate_Cu60_to_Zn60_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ni60_to_Cu60_beta_neg_weaktab
            rate_Ni60_to_Cu60_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Zn60_to_Cu60_electron_capture_weaktab
            rate_Zn60_to_Cu60_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cu61_to_Ni61_electron_capture_weaktab
            rate_Cu61_to_Ni61_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cu61_to_Zn61_beta_neg_weaktab
            rate_Cu61_to_Zn61_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ni61_to_Cu61_beta_neg_weaktab
            rate_Ni61_to_Cu61_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Zn61_to_Cu61_electron_capture_weaktab
            rate_Zn61_to_Cu61_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cu62_to_Ni62_electron_capture_weaktab
            rate_Cu62_to_Ni62_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cu62_to_Zn62_beta_neg_weaktab
            rate_Cu62_to_Zn62_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ga62_to_Zn62_electron_capture_weaktab
            rate_Ga62_to_Zn62_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ni62_to_Cu62_beta_neg_weaktab
            rate_Ni62_to_Cu62_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Zn62_to_Cu62_electron_capture_weaktab
            rate_Zn62_to_Cu62_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Zn62_to_Ga62_beta_neg_weaktab
            rate_Zn62_to_Ga62_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cu63_to_Ni63_electron_capture_weaktab
            rate_Cu63_to_Ni63_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cu63_to_Zn63_beta_neg_weaktab
            rate_Cu63_to_Zn63_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ga63_to_Zn63_electron_capture_weaktab
            rate_Ga63_to_Zn63_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ni63_to_Cu63_beta_neg_weaktab
            rate_Ni63_to_Cu63_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Zn63_to_Cu63_electron_capture_weaktab
            rate_Zn63_to_Cu63_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Zn63_to_Ga63_beta_neg_weaktab
            rate_Zn63_to_Ga63_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cu64_to_Ni64_electron_capture_weaktab
            rate_Cu64_to_Ni64_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cu64_to_Zn64_beta_neg_weaktab
            rate_Cu64_to_Zn64_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ga64_to_Ge64_beta_neg_weaktab
            rate_Ga64_to_Ge64_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ga64_to_Zn64_electron_capture_weaktab
            rate_Ga64_to_Zn64_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ge64_to_Ga64_electron_capture_weaktab
            rate_Ge64_to_Ga64_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Ni64_to_Cu64_beta_neg_weaktab
            rate_Ni64_to_Cu64_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Zn64_to_Cu64_electron_capture_weaktab
            rate_Zn64_to_Cu64_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Zn64_to_Ga64_beta_neg_weaktab
            rate_Zn64_to_Ga64_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Cu65_to_Zn65_beta_neg_weaktab
            rate_Cu65_to_Zn65_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // Zn65_to_Cu65_electron_capture_weaktab
            rate_Zn65_to_Cu65_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // n_to_p_beta_neg_weaktab
            rate_n_to_p_beta_neg_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

        {
            // p_to_n_electron_capture_weaktab
            rate_p_to_n_electron_capture_weaktab(rate_eval, temp, log_temp, log_rhoy, Y);
        }

    }

}
#endif
