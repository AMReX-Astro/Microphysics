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
    void rate_Na21_to_Ne21_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Na21 --> Ne21

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Na21_Ne21_meta, j_Na21_Ne21_rhoy, j_Na21_Ne21_temp, j_Na21_Ne21_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Na21_to_Ne21_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Na21_to_Ne21_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Na21 and Ne21
            // this is the derivative of ∂Y(Na21)/∂t = -Y(Na21) λ and
            //                           ∂Y(Ne21)/∂t = +Y(Na21) λ
            rate_eval.dweak_ydot_dYe(Na21) -= rho * Y(Na21) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ne21) += rho * Y(Na21) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Na21) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Na21) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Na21) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Na21) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ne21_to_Na21_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ne21 --> Na21

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ne21_Na21_meta, j_Ne21_Na21_rhoy, j_Ne21_Na21_temp, j_Ne21_Na21_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ne21_to_Na21_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne21_to_Na21_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ne21 and Na21
            // this is the derivative of ∂Y(Ne21)/∂t = -Y(Ne21) λ and
            //                           ∂Y(Na21)/∂t = +Y(Ne21) λ
            rate_eval.dweak_ydot_dYe(Ne21) -= rho * Y(Ne21) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Na21) += rho * Y(Ne21) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ne21) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ne21) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ne21) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ne21) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Na22_to_Ne22_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Na22 --> Ne22

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Na22_Ne22_meta, j_Na22_Ne22_rhoy, j_Na22_Ne22_temp, j_Na22_Ne22_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Na22_to_Ne22_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Na22_to_Ne22_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Na22 and Ne22
            // this is the derivative of ∂Y(Na22)/∂t = -Y(Na22) λ and
            //                           ∂Y(Ne22)/∂t = +Y(Na22) λ
            rate_eval.dweak_ydot_dYe(Na22) -= rho * Y(Na22) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ne22) += rho * Y(Na22) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Na22) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Na22) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Na22) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Na22) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ne22_to_Na22_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ne22 --> Na22

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ne22_Na22_meta, j_Ne22_Na22_rhoy, j_Ne22_Na22_temp, j_Ne22_Na22_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ne22_to_Na22_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne22_to_Na22_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ne22 and Na22
            // this is the derivative of ∂Y(Ne22)/∂t = -Y(Ne22) λ and
            //                           ∂Y(Na22)/∂t = +Y(Ne22) λ
            rate_eval.dweak_ydot_dYe(Ne22) -= rho * Y(Ne22) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Na22) += rho * Y(Ne22) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ne22) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ne22) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ne22) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ne22) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mg23_to_Na23_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Mg23 --> Na23

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Mg23_Na23_meta, j_Mg23_Na23_rhoy, j_Mg23_Na23_temp, j_Mg23_Na23_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Mg23_to_Na23_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg23_to_Na23_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Mg23 and Na23
            // this is the derivative of ∂Y(Mg23)/∂t = -Y(Mg23) λ and
            //                           ∂Y(Na23)/∂t = +Y(Mg23) λ
            rate_eval.dweak_ydot_dYe(Mg23) -= rho * Y(Mg23) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Na23) += rho * Y(Mg23) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Mg23) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Mg23) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Mg23) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Mg23) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Na23_to_Mg23_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Na23 --> Mg23

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Na23_Mg23_meta, j_Na23_Mg23_rhoy, j_Na23_Mg23_temp, j_Na23_Mg23_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Na23_to_Mg23_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Na23_to_Mg23_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Na23 and Mg23
            // this is the derivative of ∂Y(Na23)/∂t = -Y(Na23) λ and
            //                           ∂Y(Mg23)/∂t = +Y(Na23) λ
            rate_eval.dweak_ydot_dYe(Na23) -= rho * Y(Na23) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Mg23) += rho * Y(Na23) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Na23) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Na23) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Na23) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Na23) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Al25_to_Mg25_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Al25 --> Mg25

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Al25_Mg25_meta, j_Al25_Mg25_rhoy, j_Al25_Mg25_temp, j_Al25_Mg25_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Al25_to_Mg25_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Al25_to_Mg25_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Al25 and Mg25
            // this is the derivative of ∂Y(Al25)/∂t = -Y(Al25) λ and
            //                           ∂Y(Mg25)/∂t = +Y(Al25) λ
            rate_eval.dweak_ydot_dYe(Al25) -= rho * Y(Al25) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Mg25) += rho * Y(Al25) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Al25) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Al25) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Al25) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Al25) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mg25_to_Al25_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Mg25 --> Al25

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Mg25_Al25_meta, j_Mg25_Al25_rhoy, j_Mg25_Al25_temp, j_Mg25_Al25_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Mg25_to_Al25_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg25_to_Al25_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Mg25 and Al25
            // this is the derivative of ∂Y(Mg25)/∂t = -Y(Mg25) λ and
            //                           ∂Y(Al25)/∂t = +Y(Mg25) λ
            rate_eval.dweak_ydot_dYe(Mg25) -= rho * Y(Mg25) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Al25) += rho * Y(Mg25) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Mg25) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Mg25) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Mg25) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Mg25) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Al26_to_Mg26_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Al26 --> Mg26

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Al26_Mg26_meta, j_Al26_Mg26_rhoy, j_Al26_Mg26_temp, j_Al26_Mg26_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Al26_to_Mg26_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Al26_to_Mg26_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Al26 and Mg26
            // this is the derivative of ∂Y(Al26)/∂t = -Y(Al26) λ and
            //                           ∂Y(Mg26)/∂t = +Y(Al26) λ
            rate_eval.dweak_ydot_dYe(Al26) -= rho * Y(Al26) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Mg26) += rho * Y(Al26) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Al26) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Al26) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Al26) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Al26) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mg26_to_Al26_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Mg26 --> Al26

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Mg26_Al26_meta, j_Mg26_Al26_rhoy, j_Mg26_Al26_temp, j_Mg26_Al26_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Mg26_to_Al26_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg26_to_Al26_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Mg26 and Al26
            // this is the derivative of ∂Y(Mg26)/∂t = -Y(Mg26) λ and
            //                           ∂Y(Al26)/∂t = +Y(Mg26) λ
            rate_eval.dweak_ydot_dYe(Mg26) -= rho * Y(Mg26) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Al26) += rho * Y(Mg26) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Mg26) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Mg26) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Mg26) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Mg26) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P29_to_Si29_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // P29 --> Si29

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_P29_Si29_meta, j_P29_Si29_rhoy, j_P29_Si29_temp, j_P29_Si29_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_P29_to_Si29_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P29_to_Si29_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to P29 and Si29
            // this is the derivative of ∂Y(P29)/∂t = -Y(P29) λ and
            //                           ∂Y(Si29)/∂t = +Y(P29) λ
            rate_eval.dweak_ydot_dYe(P29) -= rho * Y(P29) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Si29) += rho * Y(P29) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(P29) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(P29) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(P29) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(P29) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si29_to_P29_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Si29 --> P29

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Si29_P29_meta, j_Si29_P29_rhoy, j_Si29_P29_temp, j_Si29_P29_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Si29_to_P29_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si29_to_P29_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Si29 and P29
            // this is the derivative of ∂Y(Si29)/∂t = -Y(Si29) λ and
            //                           ∂Y(P29)/∂t = +Y(Si29) λ
            rate_eval.dweak_ydot_dYe(Si29) -= rho * Y(Si29) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(P29) += rho * Y(Si29) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Si29) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Si29) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Si29) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Si29) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P30_to_Si30_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // P30 --> Si30

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_P30_Si30_meta, j_P30_Si30_rhoy, j_P30_Si30_temp, j_P30_Si30_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_P30_to_Si30_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P30_to_Si30_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to P30 and Si30
            // this is the derivative of ∂Y(P30)/∂t = -Y(P30) λ and
            //                           ∂Y(Si30)/∂t = +Y(P30) λ
            rate_eval.dweak_ydot_dYe(P30) -= rho * Y(P30) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Si30) += rho * Y(P30) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(P30) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(P30) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(P30) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(P30) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si30_to_P30_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Si30 --> P30

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Si30_P30_meta, j_Si30_P30_rhoy, j_Si30_P30_temp, j_Si30_P30_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Si30_to_P30_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si30_to_P30_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Si30 and P30
            // this is the derivative of ∂Y(Si30)/∂t = -Y(Si30) λ and
            //                           ∂Y(P30)/∂t = +Y(Si30) λ
            rate_eval.dweak_ydot_dYe(Si30) -= rho * Y(Si30) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(P30) += rho * Y(Si30) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Si30) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Si30) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Si30) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Si30) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P31_to_Si31_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // P31 --> Si31

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_P31_Si31_meta, j_P31_Si31_rhoy, j_P31_Si31_temp, j_P31_Si31_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_P31_to_Si31_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P31_to_Si31_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to P31 and Si31
            // this is the derivative of ∂Y(P31)/∂t = -Y(P31) λ and
            //                           ∂Y(Si31)/∂t = +Y(P31) λ
            rate_eval.dweak_ydot_dYe(P31) -= rho * Y(P31) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Si31) += rho * Y(P31) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(P31) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(P31) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(P31) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(P31) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si31_to_P31_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Si31 --> P31

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Si31_P31_meta, j_Si31_P31_rhoy, j_Si31_P31_temp, j_Si31_P31_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Si31_to_P31_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si31_to_P31_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Si31 and P31
            // this is the derivative of ∂Y(Si31)/∂t = -Y(Si31) λ and
            //                           ∂Y(P31)/∂t = +Y(Si31) λ
            rate_eval.dweak_ydot_dYe(Si31) -= rho * Y(Si31) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(P31) += rho * Y(Si31) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Si31) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Si31) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Si31) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Si31) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P32_to_S32_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // P32 --> S32

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_P32_S32_meta, j_P32_S32_rhoy, j_P32_S32_temp, j_P32_S32_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_P32_to_S32_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P32_to_S32_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to P32 and S32
            // this is the derivative of ∂Y(P32)/∂t = -Y(P32) λ and
            //                           ∂Y(S32)/∂t = +Y(P32) λ
            rate_eval.dweak_ydot_dYe(P32) -= rho * Y(P32) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(S32) += rho * Y(P32) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(P32) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(P32) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(P32) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(P32) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P32_to_Si32_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // P32 --> Si32

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_P32_Si32_meta, j_P32_Si32_rhoy, j_P32_Si32_temp, j_P32_Si32_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_P32_to_Si32_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P32_to_Si32_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to P32 and Si32
            // this is the derivative of ∂Y(P32)/∂t = -Y(P32) λ and
            //                           ∂Y(Si32)/∂t = +Y(P32) λ
            rate_eval.dweak_ydot_dYe(P32) -= rho * Y(P32) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Si32) += rho * Y(P32) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(P32) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(P32) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(P32) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(P32) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S32_to_P32_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // S32 --> P32

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_S32_P32_meta, j_S32_P32_rhoy, j_S32_P32_temp, j_S32_P32_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_S32_to_P32_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S32_to_P32_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to S32 and P32
            // this is the derivative of ∂Y(S32)/∂t = -Y(S32) λ and
            //                           ∂Y(P32)/∂t = +Y(S32) λ
            rate_eval.dweak_ydot_dYe(S32) -= rho * Y(S32) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(P32) += rho * Y(S32) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(S32) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(S32) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(S32) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(S32) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si32_to_P32_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Si32 --> P32

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Si32_P32_meta, j_Si32_P32_rhoy, j_Si32_P32_temp, j_Si32_P32_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Si32_to_P32_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si32_to_P32_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Si32 and P32
            // this is the derivative of ∂Y(Si32)/∂t = -Y(Si32) λ and
            //                           ∂Y(P32)/∂t = +Y(Si32) λ
            rate_eval.dweak_ydot_dYe(Si32) -= rho * Y(Si32) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(P32) += rho * Y(Si32) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Si32) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Si32) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Si32) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Si32) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl33_to_S33_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cl33 --> S33

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cl33_S33_meta, j_Cl33_S33_rhoy, j_Cl33_S33_temp, j_Cl33_S33_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cl33_to_S33_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl33_to_S33_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cl33 and S33
            // this is the derivative of ∂Y(Cl33)/∂t = -Y(Cl33) λ and
            //                           ∂Y(S33)/∂t = +Y(Cl33) λ
            rate_eval.dweak_ydot_dYe(Cl33) -= rho * Y(Cl33) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(S33) += rho * Y(Cl33) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cl33) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cl33) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cl33) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cl33) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_P33_to_S33_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // P33 --> S33

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_P33_S33_meta, j_P33_S33_rhoy, j_P33_S33_temp, j_P33_S33_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_P33_to_S33_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_P33_to_S33_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to P33 and S33
            // this is the derivative of ∂Y(P33)/∂t = -Y(P33) λ and
            //                           ∂Y(S33)/∂t = +Y(P33) λ
            rate_eval.dweak_ydot_dYe(P33) -= rho * Y(P33) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(S33) += rho * Y(P33) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(P33) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(P33) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(P33) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(P33) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S33_to_Cl33_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // S33 --> Cl33

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_S33_Cl33_meta, j_S33_Cl33_rhoy, j_S33_Cl33_temp, j_S33_Cl33_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_S33_to_Cl33_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S33_to_Cl33_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to S33 and Cl33
            // this is the derivative of ∂Y(S33)/∂t = -Y(S33) λ and
            //                           ∂Y(Cl33)/∂t = +Y(S33) λ
            rate_eval.dweak_ydot_dYe(S33) -= rho * Y(S33) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cl33) += rho * Y(S33) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(S33) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(S33) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(S33) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(S33) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S33_to_P33_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // S33 --> P33

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_S33_P33_meta, j_S33_P33_rhoy, j_S33_P33_temp, j_S33_P33_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_S33_to_P33_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S33_to_P33_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to S33 and P33
            // this is the derivative of ∂Y(S33)/∂t = -Y(S33) λ and
            //                           ∂Y(P33)/∂t = +Y(S33) λ
            rate_eval.dweak_ydot_dYe(S33) -= rho * Y(S33) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(P33) += rho * Y(S33) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(S33) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(S33) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(S33) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(S33) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl34_to_S34_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cl34 --> S34

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cl34_S34_meta, j_Cl34_S34_rhoy, j_Cl34_S34_temp, j_Cl34_S34_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cl34_to_S34_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl34_to_S34_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cl34 and S34
            // this is the derivative of ∂Y(Cl34)/∂t = -Y(Cl34) λ and
            //                           ∂Y(S34)/∂t = +Y(Cl34) λ
            rate_eval.dweak_ydot_dYe(Cl34) -= rho * Y(Cl34) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(S34) += rho * Y(Cl34) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cl34) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cl34) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cl34) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cl34) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S34_to_Cl34_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // S34 --> Cl34

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_S34_Cl34_meta, j_S34_Cl34_rhoy, j_S34_Cl34_temp, j_S34_Cl34_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_S34_to_Cl34_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S34_to_Cl34_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to S34 and Cl34
            // this is the derivative of ∂Y(S34)/∂t = -Y(S34) λ and
            //                           ∂Y(Cl34)/∂t = +Y(S34) λ
            rate_eval.dweak_ydot_dYe(S34) -= rho * Y(S34) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cl34) += rho * Y(S34) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(S34) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(S34) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(S34) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(S34) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl35_to_S35_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cl35 --> S35

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cl35_S35_meta, j_Cl35_S35_rhoy, j_Cl35_S35_temp, j_Cl35_S35_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cl35_to_S35_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl35_to_S35_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cl35 and S35
            // this is the derivative of ∂Y(Cl35)/∂t = -Y(Cl35) λ and
            //                           ∂Y(S35)/∂t = +Y(Cl35) λ
            rate_eval.dweak_ydot_dYe(Cl35) -= rho * Y(Cl35) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(S35) += rho * Y(Cl35) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cl35) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cl35) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cl35) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cl35) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S35_to_Cl35_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // S35 --> Cl35

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_S35_Cl35_meta, j_S35_Cl35_rhoy, j_S35_Cl35_temp, j_S35_Cl35_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_S35_to_Cl35_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S35_to_Cl35_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to S35 and Cl35
            // this is the derivative of ∂Y(S35)/∂t = -Y(S35) λ and
            //                           ∂Y(Cl35)/∂t = +Y(S35) λ
            rate_eval.dweak_ydot_dYe(S35) -= rho * Y(S35) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cl35) += rho * Y(S35) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(S35) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(S35) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(S35) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(S35) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar36_to_Cl36_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ar36 --> Cl36

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ar36_Cl36_meta, j_Ar36_Cl36_rhoy, j_Ar36_Cl36_temp, j_Ar36_Cl36_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ar36_to_Cl36_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar36_to_Cl36_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ar36 and Cl36
            // this is the derivative of ∂Y(Ar36)/∂t = -Y(Ar36) λ and
            //                           ∂Y(Cl36)/∂t = +Y(Ar36) λ
            rate_eval.dweak_ydot_dYe(Ar36) -= rho * Y(Ar36) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cl36) += rho * Y(Ar36) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ar36) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ar36) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ar36) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ar36) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl36_to_Ar36_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cl36 --> Ar36

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cl36_Ar36_meta, j_Cl36_Ar36_rhoy, j_Cl36_Ar36_temp, j_Cl36_Ar36_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cl36_to_Ar36_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl36_to_Ar36_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cl36 and Ar36
            // this is the derivative of ∂Y(Cl36)/∂t = -Y(Cl36) λ and
            //                           ∂Y(Ar36)/∂t = +Y(Cl36) λ
            rate_eval.dweak_ydot_dYe(Cl36) -= rho * Y(Cl36) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ar36) += rho * Y(Cl36) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cl36) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cl36) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cl36) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cl36) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl36_to_S36_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cl36 --> S36

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cl36_S36_meta, j_Cl36_S36_rhoy, j_Cl36_S36_temp, j_Cl36_S36_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cl36_to_S36_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl36_to_S36_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cl36 and S36
            // this is the derivative of ∂Y(Cl36)/∂t = -Y(Cl36) λ and
            //                           ∂Y(S36)/∂t = +Y(Cl36) λ
            rate_eval.dweak_ydot_dYe(Cl36) -= rho * Y(Cl36) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(S36) += rho * Y(Cl36) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cl36) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cl36) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cl36) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cl36) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S36_to_Cl36_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // S36 --> Cl36

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_S36_Cl36_meta, j_S36_Cl36_rhoy, j_S36_Cl36_temp, j_S36_Cl36_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_S36_to_Cl36_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S36_to_Cl36_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to S36 and Cl36
            // this is the derivative of ∂Y(S36)/∂t = -Y(S36) λ and
            //                           ∂Y(Cl36)/∂t = +Y(S36) λ
            rate_eval.dweak_ydot_dYe(S36) -= rho * Y(S36) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cl36) += rho * Y(S36) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(S36) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(S36) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(S36) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(S36) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar37_to_Cl37_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ar37 --> Cl37

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ar37_Cl37_meta, j_Ar37_Cl37_rhoy, j_Ar37_Cl37_temp, j_Ar37_Cl37_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ar37_to_Cl37_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar37_to_Cl37_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ar37 and Cl37
            // this is the derivative of ∂Y(Ar37)/∂t = -Y(Ar37) λ and
            //                           ∂Y(Cl37)/∂t = +Y(Ar37) λ
            rate_eval.dweak_ydot_dYe(Ar37) -= rho * Y(Ar37) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cl37) += rho * Y(Ar37) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ar37) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ar37) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ar37) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ar37) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar37_to_K37_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ar37 --> K37

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ar37_K37_meta, j_Ar37_K37_rhoy, j_Ar37_K37_temp, j_Ar37_K37_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ar37_to_K37_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar37_to_K37_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ar37 and K37
            // this is the derivative of ∂Y(Ar37)/∂t = -Y(Ar37) λ and
            //                           ∂Y(K37)/∂t = +Y(Ar37) λ
            rate_eval.dweak_ydot_dYe(Ar37) -= rho * Y(Ar37) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(K37) += rho * Y(Ar37) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ar37) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ar37) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ar37) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ar37) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cl37_to_Ar37_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cl37 --> Ar37

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cl37_Ar37_meta, j_Cl37_Ar37_rhoy, j_Cl37_Ar37_temp, j_Cl37_Ar37_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cl37_to_Ar37_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cl37_to_Ar37_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cl37 and Ar37
            // this is the derivative of ∂Y(Cl37)/∂t = -Y(Cl37) λ and
            //                           ∂Y(Ar37)/∂t = +Y(Cl37) λ
            rate_eval.dweak_ydot_dYe(Cl37) -= rho * Y(Cl37) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ar37) += rho * Y(Cl37) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cl37) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cl37) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cl37) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cl37) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K37_to_Ar37_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // K37 --> Ar37

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_K37_Ar37_meta, j_K37_Ar37_rhoy, j_K37_Ar37_temp, j_K37_Ar37_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_K37_to_Ar37_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K37_to_Ar37_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to K37 and Ar37
            // this is the derivative of ∂Y(K37)/∂t = -Y(K37) λ and
            //                           ∂Y(Ar37)/∂t = +Y(K37) λ
            rate_eval.dweak_ydot_dYe(K37) -= rho * Y(K37) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ar37) += rho * Y(K37) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(K37) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(K37) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(K37) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(K37) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar38_to_K38_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ar38 --> K38

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ar38_K38_meta, j_Ar38_K38_rhoy, j_Ar38_K38_temp, j_Ar38_K38_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ar38_to_K38_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar38_to_K38_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ar38 and K38
            // this is the derivative of ∂Y(Ar38)/∂t = -Y(Ar38) λ and
            //                           ∂Y(K38)/∂t = +Y(Ar38) λ
            rate_eval.dweak_ydot_dYe(Ar38) -= rho * Y(Ar38) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(K38) += rho * Y(Ar38) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ar38) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ar38) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ar38) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ar38) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K38_to_Ar38_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // K38 --> Ar38

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_K38_Ar38_meta, j_K38_Ar38_rhoy, j_K38_Ar38_temp, j_K38_Ar38_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_K38_to_Ar38_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K38_to_Ar38_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to K38 and Ar38
            // this is the derivative of ∂Y(K38)/∂t = -Y(K38) λ and
            //                           ∂Y(Ar38)/∂t = +Y(K38) λ
            rate_eval.dweak_ydot_dYe(K38) -= rho * Y(K38) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ar38) += rho * Y(K38) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(K38) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(K38) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(K38) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(K38) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar39_to_K39_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ar39 --> K39

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ar39_K39_meta, j_Ar39_K39_rhoy, j_Ar39_K39_temp, j_Ar39_K39_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ar39_to_K39_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar39_to_K39_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ar39 and K39
            // this is the derivative of ∂Y(Ar39)/∂t = -Y(Ar39) λ and
            //                           ∂Y(K39)/∂t = +Y(Ar39) λ
            rate_eval.dweak_ydot_dYe(Ar39) -= rho * Y(Ar39) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(K39) += rho * Y(Ar39) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ar39) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ar39) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ar39) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ar39) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K39_to_Ar39_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // K39 --> Ar39

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_K39_Ar39_meta, j_K39_Ar39_rhoy, j_K39_Ar39_temp, j_K39_Ar39_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_K39_to_Ar39_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K39_to_Ar39_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to K39 and Ar39
            // this is the derivative of ∂Y(K39)/∂t = -Y(K39) λ and
            //                           ∂Y(Ar39)/∂t = +Y(K39) λ
            rate_eval.dweak_ydot_dYe(K39) -= rho * Y(K39) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ar39) += rho * Y(K39) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(K39) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(K39) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(K39) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(K39) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar40_to_K40_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ar40 --> K40

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ar40_K40_meta, j_Ar40_K40_rhoy, j_Ar40_K40_temp, j_Ar40_K40_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ar40_to_K40_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar40_to_K40_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ar40 and K40
            // this is the derivative of ∂Y(Ar40)/∂t = -Y(Ar40) λ and
            //                           ∂Y(K40)/∂t = +Y(Ar40) λ
            rate_eval.dweak_ydot_dYe(Ar40) -= rho * Y(Ar40) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(K40) += rho * Y(Ar40) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ar40) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ar40) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ar40) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ar40) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca40_to_K40_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ca40 --> K40

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ca40_K40_meta, j_Ca40_K40_rhoy, j_Ca40_K40_temp, j_Ca40_K40_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ca40_to_K40_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca40_to_K40_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ca40 and K40
            // this is the derivative of ∂Y(Ca40)/∂t = -Y(Ca40) λ and
            //                           ∂Y(K40)/∂t = +Y(Ca40) λ
            rate_eval.dweak_ydot_dYe(Ca40) -= rho * Y(Ca40) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(K40) += rho * Y(Ca40) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ca40) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ca40) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ca40) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ca40) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K40_to_Ar40_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // K40 --> Ar40

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_K40_Ar40_meta, j_K40_Ar40_rhoy, j_K40_Ar40_temp, j_K40_Ar40_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_K40_to_Ar40_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K40_to_Ar40_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to K40 and Ar40
            // this is the derivative of ∂Y(K40)/∂t = -Y(K40) λ and
            //                           ∂Y(Ar40)/∂t = +Y(K40) λ
            rate_eval.dweak_ydot_dYe(K40) -= rho * Y(K40) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ar40) += rho * Y(K40) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(K40) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(K40) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(K40) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(K40) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K40_to_Ca40_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // K40 --> Ca40

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_K40_Ca40_meta, j_K40_Ca40_rhoy, j_K40_Ca40_temp, j_K40_Ca40_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_K40_to_Ca40_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K40_to_Ca40_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to K40 and Ca40
            // this is the derivative of ∂Y(K40)/∂t = -Y(K40) λ and
            //                           ∂Y(Ca40)/∂t = +Y(K40) λ
            rate_eval.dweak_ydot_dYe(K40) -= rho * Y(K40) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ca40) += rho * Y(K40) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(K40) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(K40) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(K40) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(K40) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca41_to_K41_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ca41 --> K41

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ca41_K41_meta, j_Ca41_K41_rhoy, j_Ca41_K41_temp, j_Ca41_K41_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ca41_to_K41_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca41_to_K41_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ca41 and K41
            // this is the derivative of ∂Y(Ca41)/∂t = -Y(Ca41) λ and
            //                           ∂Y(K41)/∂t = +Y(Ca41) λ
            rate_eval.dweak_ydot_dYe(Ca41) -= rho * Y(Ca41) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(K41) += rho * Y(Ca41) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ca41) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ca41) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ca41) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ca41) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_K41_to_Ca41_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // K41 --> Ca41

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_K41_Ca41_meta, j_K41_Ca41_rhoy, j_K41_Ca41_temp, j_K41_Ca41_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_K41_to_Ca41_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_K41_to_Ca41_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to K41 and Ca41
            // this is the derivative of ∂Y(K41)/∂t = -Y(K41) λ and
            //                           ∂Y(Ca41)/∂t = +Y(K41) λ
            rate_eval.dweak_ydot_dYe(K41) -= rho * Y(K41) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ca41) += rho * Y(K41) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(K41) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(K41) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(K41) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(K41) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca43_to_Sc43_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ca43 --> Sc43

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ca43_Sc43_meta, j_Ca43_Sc43_rhoy, j_Ca43_Sc43_temp, j_Ca43_Sc43_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ca43_to_Sc43_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca43_to_Sc43_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ca43 and Sc43
            // this is the derivative of ∂Y(Ca43)/∂t = -Y(Ca43) λ and
            //                           ∂Y(Sc43)/∂t = +Y(Ca43) λ
            rate_eval.dweak_ydot_dYe(Ca43) -= rho * Y(Ca43) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Sc43) += rho * Y(Ca43) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ca43) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ca43) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ca43) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ca43) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc43_to_Ca43_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Sc43 --> Ca43

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Sc43_Ca43_meta, j_Sc43_Ca43_rhoy, j_Sc43_Ca43_temp, j_Sc43_Ca43_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Sc43_to_Ca43_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc43_to_Ca43_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Sc43 and Ca43
            // this is the derivative of ∂Y(Sc43)/∂t = -Y(Sc43) λ and
            //                           ∂Y(Ca43)/∂t = +Y(Sc43) λ
            rate_eval.dweak_ydot_dYe(Sc43) -= rho * Y(Sc43) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ca43) += rho * Y(Sc43) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Sc43) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Sc43) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Sc43) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Sc43) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca44_to_Sc44_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ca44 --> Sc44

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ca44_Sc44_meta, j_Ca44_Sc44_rhoy, j_Ca44_Sc44_temp, j_Ca44_Sc44_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ca44_to_Sc44_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca44_to_Sc44_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ca44 and Sc44
            // this is the derivative of ∂Y(Ca44)/∂t = -Y(Ca44) λ and
            //                           ∂Y(Sc44)/∂t = +Y(Ca44) λ
            rate_eval.dweak_ydot_dYe(Ca44) -= rho * Y(Ca44) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Sc44) += rho * Y(Ca44) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ca44) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ca44) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ca44) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ca44) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc44_to_Ca44_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Sc44 --> Ca44

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Sc44_Ca44_meta, j_Sc44_Ca44_rhoy, j_Sc44_Ca44_temp, j_Sc44_Ca44_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Sc44_to_Ca44_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc44_to_Ca44_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Sc44 and Ca44
            // this is the derivative of ∂Y(Sc44)/∂t = -Y(Sc44) λ and
            //                           ∂Y(Ca44)/∂t = +Y(Sc44) λ
            rate_eval.dweak_ydot_dYe(Sc44) -= rho * Y(Sc44) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ca44) += rho * Y(Sc44) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Sc44) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Sc44) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Sc44) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Sc44) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc44_to_Ti44_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Sc44 --> Ti44

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Sc44_Ti44_meta, j_Sc44_Ti44_rhoy, j_Sc44_Ti44_temp, j_Sc44_Ti44_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Sc44_to_Ti44_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc44_to_Ti44_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Sc44 and Ti44
            // this is the derivative of ∂Y(Sc44)/∂t = -Y(Sc44) λ and
            //                           ∂Y(Ti44)/∂t = +Y(Sc44) λ
            rate_eval.dweak_ydot_dYe(Sc44) -= rho * Y(Sc44) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ti44) += rho * Y(Sc44) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Sc44) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Sc44) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Sc44) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Sc44) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti44_to_Sc44_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ti44 --> Sc44

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ti44_Sc44_meta, j_Ti44_Sc44_rhoy, j_Ti44_Sc44_temp, j_Ti44_Sc44_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ti44_to_Sc44_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti44_to_Sc44_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ti44 and Sc44
            // this is the derivative of ∂Y(Ti44)/∂t = -Y(Ti44) λ and
            //                           ∂Y(Sc44)/∂t = +Y(Ti44) λ
            rate_eval.dweak_ydot_dYe(Ti44) -= rho * Y(Ti44) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Sc44) += rho * Y(Ti44) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ti44) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ti44) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ti44) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ti44) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co53_to_Fe53_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Co53 --> Fe53

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Co53_Fe53_meta, j_Co53_Fe53_rhoy, j_Co53_Fe53_temp, j_Co53_Fe53_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Co53_to_Fe53_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co53_to_Fe53_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Co53 and Fe53
            // this is the derivative of ∂Y(Co53)/∂t = -Y(Co53) λ and
            //                           ∂Y(Fe53)/∂t = +Y(Co53) λ
            rate_eval.dweak_ydot_dYe(Co53) -= rho * Y(Co53) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Fe53) += rho * Y(Co53) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Co53) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Co53) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Co53) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Co53) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe53_to_Co53_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Fe53 --> Co53

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Fe53_Co53_meta, j_Fe53_Co53_rhoy, j_Fe53_Co53_temp, j_Fe53_Co53_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Fe53_to_Co53_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe53_to_Co53_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Fe53 and Co53
            // this is the derivative of ∂Y(Fe53)/∂t = -Y(Fe53) λ and
            //                           ∂Y(Co53)/∂t = +Y(Fe53) λ
            rate_eval.dweak_ydot_dYe(Fe53) -= rho * Y(Fe53) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Co53) += rho * Y(Fe53) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Fe53) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Fe53) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Fe53) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Fe53) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu57_to_Ni57_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu57 --> Ni57

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu57_Ni57_meta, j_Cu57_Ni57_rhoy, j_Cu57_Ni57_temp, j_Cu57_Ni57_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu57_to_Ni57_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu57_to_Ni57_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cu57 and Ni57
            // this is the derivative of ∂Y(Cu57)/∂t = -Y(Cu57) λ and
            //                           ∂Y(Ni57)/∂t = +Y(Cu57) λ
            rate_eval.dweak_ydot_dYe(Cu57) -= rho * Y(Cu57) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ni57) += rho * Y(Cu57) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cu57) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cu57) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cu57) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu57) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni57_to_Cu57_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ni57 --> Cu57

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ni57_Cu57_meta, j_Ni57_Cu57_rhoy, j_Ni57_Cu57_temp, j_Ni57_Cu57_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ni57_to_Cu57_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni57_to_Cu57_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ni57 and Cu57
            // this is the derivative of ∂Y(Ni57)/∂t = -Y(Ni57) λ and
            //                           ∂Y(Cu57)/∂t = +Y(Ni57) λ
            rate_eval.dweak_ydot_dYe(Ni57) -= rho * Y(Ni57) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cu57) += rho * Y(Ni57) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ni57) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ni57) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ni57) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ni57) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca45_to_Sc45_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ca45 --> Sc45

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ca45_Sc45_meta, j_Ca45_Sc45_rhoy, j_Ca45_Sc45_temp, j_Ca45_Sc45_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ca45_to_Sc45_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca45_to_Sc45_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ca45 and Sc45
            // this is the derivative of ∂Y(Ca45)/∂t = -Y(Ca45) λ and
            //                           ∂Y(Sc45)/∂t = +Y(Ca45) λ
            rate_eval.dweak_ydot_dYe(Ca45) -= rho * Y(Ca45) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Sc45) += rho * Y(Ca45) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ca45) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ca45) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ca45) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ca45) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc45_to_Ca45_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Sc45 --> Ca45

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Sc45_Ca45_meta, j_Sc45_Ca45_rhoy, j_Sc45_Ca45_temp, j_Sc45_Ca45_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Sc45_to_Ca45_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc45_to_Ca45_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Sc45 and Ca45
            // this is the derivative of ∂Y(Sc45)/∂t = -Y(Sc45) λ and
            //                           ∂Y(Ca45)/∂t = +Y(Sc45) λ
            rate_eval.dweak_ydot_dYe(Sc45) -= rho * Y(Sc45) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ca45) += rho * Y(Sc45) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Sc45) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Sc45) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Sc45) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Sc45) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc45_to_Ti45_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Sc45 --> Ti45

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Sc45_Ti45_meta, j_Sc45_Ti45_rhoy, j_Sc45_Ti45_temp, j_Sc45_Ti45_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Sc45_to_Ti45_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc45_to_Ti45_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Sc45 and Ti45
            // this is the derivative of ∂Y(Sc45)/∂t = -Y(Sc45) λ and
            //                           ∂Y(Ti45)/∂t = +Y(Sc45) λ
            rate_eval.dweak_ydot_dYe(Sc45) -= rho * Y(Sc45) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ti45) += rho * Y(Sc45) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Sc45) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Sc45) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Sc45) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Sc45) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti45_to_Sc45_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ti45 --> Sc45

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ti45_Sc45_meta, j_Ti45_Sc45_rhoy, j_Ti45_Sc45_temp, j_Ti45_Sc45_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ti45_to_Sc45_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti45_to_Sc45_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ti45 and Sc45
            // this is the derivative of ∂Y(Ti45)/∂t = -Y(Ti45) λ and
            //                           ∂Y(Sc45)/∂t = +Y(Ti45) λ
            rate_eval.dweak_ydot_dYe(Ti45) -= rho * Y(Ti45) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Sc45) += rho * Y(Ti45) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ti45) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ti45) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ti45) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ti45) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca46_to_Sc46_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ca46 --> Sc46

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ca46_Sc46_meta, j_Ca46_Sc46_rhoy, j_Ca46_Sc46_temp, j_Ca46_Sc46_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ca46_to_Sc46_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca46_to_Sc46_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ca46 and Sc46
            // this is the derivative of ∂Y(Ca46)/∂t = -Y(Ca46) λ and
            //                           ∂Y(Sc46)/∂t = +Y(Ca46) λ
            rate_eval.dweak_ydot_dYe(Ca46) -= rho * Y(Ca46) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Sc46) += rho * Y(Ca46) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ca46) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ca46) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ca46) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ca46) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc46_to_Ca46_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Sc46 --> Ca46

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Sc46_Ca46_meta, j_Sc46_Ca46_rhoy, j_Sc46_Ca46_temp, j_Sc46_Ca46_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Sc46_to_Ca46_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc46_to_Ca46_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Sc46 and Ca46
            // this is the derivative of ∂Y(Sc46)/∂t = -Y(Sc46) λ and
            //                           ∂Y(Ca46)/∂t = +Y(Sc46) λ
            rate_eval.dweak_ydot_dYe(Sc46) -= rho * Y(Sc46) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ca46) += rho * Y(Sc46) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Sc46) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Sc46) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Sc46) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Sc46) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc46_to_Ti46_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Sc46 --> Ti46

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Sc46_Ti46_meta, j_Sc46_Ti46_rhoy, j_Sc46_Ti46_temp, j_Sc46_Ti46_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Sc46_to_Ti46_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc46_to_Ti46_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Sc46 and Ti46
            // this is the derivative of ∂Y(Sc46)/∂t = -Y(Sc46) λ and
            //                           ∂Y(Ti46)/∂t = +Y(Sc46) λ
            rate_eval.dweak_ydot_dYe(Sc46) -= rho * Y(Sc46) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ti46) += rho * Y(Sc46) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Sc46) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Sc46) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Sc46) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Sc46) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti46_to_Sc46_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ti46 --> Sc46

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ti46_Sc46_meta, j_Ti46_Sc46_rhoy, j_Ti46_Sc46_temp, j_Ti46_Sc46_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ti46_to_Sc46_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti46_to_Sc46_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ti46 and Sc46
            // this is the derivative of ∂Y(Ti46)/∂t = -Y(Ti46) λ and
            //                           ∂Y(Sc46)/∂t = +Y(Ti46) λ
            rate_eval.dweak_ydot_dYe(Ti46) -= rho * Y(Ti46) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Sc46) += rho * Y(Ti46) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ti46) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ti46) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ti46) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ti46) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti46_to_V46_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ti46 --> V46

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ti46_V46_meta, j_Ti46_V46_rhoy, j_Ti46_V46_temp, j_Ti46_V46_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ti46_to_V46_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti46_to_V46_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ti46 and V46
            // this is the derivative of ∂Y(Ti46)/∂t = -Y(Ti46) λ and
            //                           ∂Y(V46)/∂t = +Y(Ti46) λ
            rate_eval.dweak_ydot_dYe(Ti46) -= rho * Y(Ti46) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(V46) += rho * Y(Ti46) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ti46) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ti46) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ti46) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ti46) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V46_to_Ti46_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // V46 --> Ti46

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_V46_Ti46_meta, j_V46_Ti46_rhoy, j_V46_Ti46_temp, j_V46_Ti46_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_V46_to_Ti46_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V46_to_Ti46_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to V46 and Ti46
            // this is the derivative of ∂Y(V46)/∂t = -Y(V46) λ and
            //                           ∂Y(Ti46)/∂t = +Y(V46) λ
            rate_eval.dweak_ydot_dYe(V46) -= rho * Y(V46) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ti46) += rho * Y(V46) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(V46) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(V46) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(V46) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(V46) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca47_to_Sc47_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ca47 --> Sc47

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ca47_Sc47_meta, j_Ca47_Sc47_rhoy, j_Ca47_Sc47_temp, j_Ca47_Sc47_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ca47_to_Sc47_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca47_to_Sc47_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ca47 and Sc47
            // this is the derivative of ∂Y(Ca47)/∂t = -Y(Ca47) λ and
            //                           ∂Y(Sc47)/∂t = +Y(Ca47) λ
            rate_eval.dweak_ydot_dYe(Ca47) -= rho * Y(Ca47) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Sc47) += rho * Y(Ca47) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ca47) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ca47) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ca47) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ca47) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc47_to_Ca47_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Sc47 --> Ca47

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Sc47_Ca47_meta, j_Sc47_Ca47_rhoy, j_Sc47_Ca47_temp, j_Sc47_Ca47_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Sc47_to_Ca47_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc47_to_Ca47_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Sc47 and Ca47
            // this is the derivative of ∂Y(Sc47)/∂t = -Y(Sc47) λ and
            //                           ∂Y(Ca47)/∂t = +Y(Sc47) λ
            rate_eval.dweak_ydot_dYe(Sc47) -= rho * Y(Sc47) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ca47) += rho * Y(Sc47) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Sc47) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Sc47) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Sc47) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Sc47) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc47_to_Ti47_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Sc47 --> Ti47

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Sc47_Ti47_meta, j_Sc47_Ti47_rhoy, j_Sc47_Ti47_temp, j_Sc47_Ti47_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Sc47_to_Ti47_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc47_to_Ti47_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Sc47 and Ti47
            // this is the derivative of ∂Y(Sc47)/∂t = -Y(Sc47) λ and
            //                           ∂Y(Ti47)/∂t = +Y(Sc47) λ
            rate_eval.dweak_ydot_dYe(Sc47) -= rho * Y(Sc47) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ti47) += rho * Y(Sc47) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Sc47) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Sc47) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Sc47) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Sc47) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti47_to_Sc47_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ti47 --> Sc47

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ti47_Sc47_meta, j_Ti47_Sc47_rhoy, j_Ti47_Sc47_temp, j_Ti47_Sc47_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ti47_to_Sc47_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti47_to_Sc47_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ti47 and Sc47
            // this is the derivative of ∂Y(Ti47)/∂t = -Y(Ti47) λ and
            //                           ∂Y(Sc47)/∂t = +Y(Ti47) λ
            rate_eval.dweak_ydot_dYe(Ti47) -= rho * Y(Ti47) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Sc47) += rho * Y(Ti47) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ti47) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ti47) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ti47) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ti47) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti47_to_V47_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ti47 --> V47

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ti47_V47_meta, j_Ti47_V47_rhoy, j_Ti47_V47_temp, j_Ti47_V47_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ti47_to_V47_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti47_to_V47_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ti47 and V47
            // this is the derivative of ∂Y(Ti47)/∂t = -Y(Ti47) λ and
            //                           ∂Y(V47)/∂t = +Y(Ti47) λ
            rate_eval.dweak_ydot_dYe(Ti47) -= rho * Y(Ti47) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(V47) += rho * Y(Ti47) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ti47) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ti47) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ti47) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ti47) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V47_to_Ti47_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // V47 --> Ti47

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_V47_Ti47_meta, j_V47_Ti47_rhoy, j_V47_Ti47_temp, j_V47_Ti47_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_V47_to_Ti47_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V47_to_Ti47_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to V47 and Ti47
            // this is the derivative of ∂Y(V47)/∂t = -Y(V47) λ and
            //                           ∂Y(Ti47)/∂t = +Y(V47) λ
            rate_eval.dweak_ydot_dYe(V47) -= rho * Y(V47) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ti47) += rho * Y(V47) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(V47) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(V47) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(V47) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(V47) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca48_to_Sc48_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ca48 --> Sc48

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ca48_Sc48_meta, j_Ca48_Sc48_rhoy, j_Ca48_Sc48_temp, j_Ca48_Sc48_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ca48_to_Sc48_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca48_to_Sc48_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ca48 and Sc48
            // this is the derivative of ∂Y(Ca48)/∂t = -Y(Ca48) λ and
            //                           ∂Y(Sc48)/∂t = +Y(Ca48) λ
            rate_eval.dweak_ydot_dYe(Ca48) -= rho * Y(Ca48) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Sc48) += rho * Y(Ca48) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ca48) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ca48) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ca48) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ca48) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr48_to_V48_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cr48 --> V48

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cr48_V48_meta, j_Cr48_V48_rhoy, j_Cr48_V48_temp, j_Cr48_V48_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cr48_to_V48_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr48_to_V48_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cr48 and V48
            // this is the derivative of ∂Y(Cr48)/∂t = -Y(Cr48) λ and
            //                           ∂Y(V48)/∂t = +Y(Cr48) λ
            rate_eval.dweak_ydot_dYe(Cr48) -= rho * Y(Cr48) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(V48) += rho * Y(Cr48) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cr48) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cr48) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cr48) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cr48) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc48_to_Ca48_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Sc48 --> Ca48

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Sc48_Ca48_meta, j_Sc48_Ca48_rhoy, j_Sc48_Ca48_temp, j_Sc48_Ca48_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Sc48_to_Ca48_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc48_to_Ca48_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Sc48 and Ca48
            // this is the derivative of ∂Y(Sc48)/∂t = -Y(Sc48) λ and
            //                           ∂Y(Ca48)/∂t = +Y(Sc48) λ
            rate_eval.dweak_ydot_dYe(Sc48) -= rho * Y(Sc48) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ca48) += rho * Y(Sc48) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Sc48) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Sc48) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Sc48) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Sc48) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc48_to_Ti48_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Sc48 --> Ti48

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Sc48_Ti48_meta, j_Sc48_Ti48_rhoy, j_Sc48_Ti48_temp, j_Sc48_Ti48_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Sc48_to_Ti48_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc48_to_Ti48_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Sc48 and Ti48
            // this is the derivative of ∂Y(Sc48)/∂t = -Y(Sc48) λ and
            //                           ∂Y(Ti48)/∂t = +Y(Sc48) λ
            rate_eval.dweak_ydot_dYe(Sc48) -= rho * Y(Sc48) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ti48) += rho * Y(Sc48) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Sc48) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Sc48) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Sc48) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Sc48) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti48_to_Sc48_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ti48 --> Sc48

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ti48_Sc48_meta, j_Ti48_Sc48_rhoy, j_Ti48_Sc48_temp, j_Ti48_Sc48_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ti48_to_Sc48_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti48_to_Sc48_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ti48 and Sc48
            // this is the derivative of ∂Y(Ti48)/∂t = -Y(Ti48) λ and
            //                           ∂Y(Sc48)/∂t = +Y(Ti48) λ
            rate_eval.dweak_ydot_dYe(Ti48) -= rho * Y(Ti48) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Sc48) += rho * Y(Ti48) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ti48) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ti48) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ti48) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ti48) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti48_to_V48_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ti48 --> V48

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ti48_V48_meta, j_Ti48_V48_rhoy, j_Ti48_V48_temp, j_Ti48_V48_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ti48_to_V48_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti48_to_V48_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ti48 and V48
            // this is the derivative of ∂Y(Ti48)/∂t = -Y(Ti48) λ and
            //                           ∂Y(V48)/∂t = +Y(Ti48) λ
            rate_eval.dweak_ydot_dYe(Ti48) -= rho * Y(Ti48) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(V48) += rho * Y(Ti48) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ti48) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ti48) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ti48) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ti48) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V48_to_Cr48_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // V48 --> Cr48

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_V48_Cr48_meta, j_V48_Cr48_rhoy, j_V48_Cr48_temp, j_V48_Cr48_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_V48_to_Cr48_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V48_to_Cr48_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to V48 and Cr48
            // this is the derivative of ∂Y(V48)/∂t = -Y(V48) λ and
            //                           ∂Y(Cr48)/∂t = +Y(V48) λ
            rate_eval.dweak_ydot_dYe(V48) -= rho * Y(V48) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cr48) += rho * Y(V48) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(V48) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(V48) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(V48) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(V48) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V48_to_Ti48_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // V48 --> Ti48

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_V48_Ti48_meta, j_V48_Ti48_rhoy, j_V48_Ti48_temp, j_V48_Ti48_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_V48_to_Ti48_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V48_to_Ti48_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to V48 and Ti48
            // this is the derivative of ∂Y(V48)/∂t = -Y(V48) λ and
            //                           ∂Y(Ti48)/∂t = +Y(V48) λ
            rate_eval.dweak_ydot_dYe(V48) -= rho * Y(V48) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ti48) += rho * Y(V48) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(V48) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(V48) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(V48) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(V48) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr49_to_V49_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cr49 --> V49

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cr49_V49_meta, j_Cr49_V49_rhoy, j_Cr49_V49_temp, j_Cr49_V49_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cr49_to_V49_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr49_to_V49_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cr49 and V49
            // this is the derivative of ∂Y(Cr49)/∂t = -Y(Cr49) λ and
            //                           ∂Y(V49)/∂t = +Y(Cr49) λ
            rate_eval.dweak_ydot_dYe(Cr49) -= rho * Y(Cr49) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(V49) += rho * Y(Cr49) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cr49) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cr49) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cr49) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cr49) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc49_to_Ti49_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Sc49 --> Ti49

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Sc49_Ti49_meta, j_Sc49_Ti49_rhoy, j_Sc49_Ti49_temp, j_Sc49_Ti49_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Sc49_to_Ti49_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc49_to_Ti49_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Sc49 and Ti49
            // this is the derivative of ∂Y(Sc49)/∂t = -Y(Sc49) λ and
            //                           ∂Y(Ti49)/∂t = +Y(Sc49) λ
            rate_eval.dweak_ydot_dYe(Sc49) -= rho * Y(Sc49) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ti49) += rho * Y(Sc49) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Sc49) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Sc49) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Sc49) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Sc49) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti49_to_Sc49_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ti49 --> Sc49

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ti49_Sc49_meta, j_Ti49_Sc49_rhoy, j_Ti49_Sc49_temp, j_Ti49_Sc49_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ti49_to_Sc49_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti49_to_Sc49_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ti49 and Sc49
            // this is the derivative of ∂Y(Ti49)/∂t = -Y(Ti49) λ and
            //                           ∂Y(Sc49)/∂t = +Y(Ti49) λ
            rate_eval.dweak_ydot_dYe(Ti49) -= rho * Y(Ti49) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Sc49) += rho * Y(Ti49) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ti49) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ti49) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ti49) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ti49) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti49_to_V49_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ti49 --> V49

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ti49_V49_meta, j_Ti49_V49_rhoy, j_Ti49_V49_temp, j_Ti49_V49_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ti49_to_V49_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti49_to_V49_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ti49 and V49
            // this is the derivative of ∂Y(Ti49)/∂t = -Y(Ti49) λ and
            //                           ∂Y(V49)/∂t = +Y(Ti49) λ
            rate_eval.dweak_ydot_dYe(Ti49) -= rho * Y(Ti49) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(V49) += rho * Y(Ti49) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ti49) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ti49) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ti49) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ti49) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V49_to_Cr49_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // V49 --> Cr49

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_V49_Cr49_meta, j_V49_Cr49_rhoy, j_V49_Cr49_temp, j_V49_Cr49_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_V49_to_Cr49_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V49_to_Cr49_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to V49 and Cr49
            // this is the derivative of ∂Y(V49)/∂t = -Y(V49) λ and
            //                           ∂Y(Cr49)/∂t = +Y(V49) λ
            rate_eval.dweak_ydot_dYe(V49) -= rho * Y(V49) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cr49) += rho * Y(V49) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(V49) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(V49) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(V49) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(V49) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V49_to_Ti49_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // V49 --> Ti49

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_V49_Ti49_meta, j_V49_Ti49_rhoy, j_V49_Ti49_temp, j_V49_Ti49_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_V49_to_Ti49_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V49_to_Ti49_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to V49 and Ti49
            // this is the derivative of ∂Y(V49)/∂t = -Y(V49) λ and
            //                           ∂Y(Ti49)/∂t = +Y(V49) λ
            rate_eval.dweak_ydot_dYe(V49) -= rho * Y(V49) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ti49) += rho * Y(V49) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(V49) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(V49) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(V49) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(V49) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr50_to_Mn50_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cr50 --> Mn50

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cr50_Mn50_meta, j_Cr50_Mn50_rhoy, j_Cr50_Mn50_temp, j_Cr50_Mn50_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cr50_to_Mn50_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr50_to_Mn50_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cr50 and Mn50
            // this is the derivative of ∂Y(Cr50)/∂t = -Y(Cr50) λ and
            //                           ∂Y(Mn50)/∂t = +Y(Cr50) λ
            rate_eval.dweak_ydot_dYe(Cr50) -= rho * Y(Cr50) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Mn50) += rho * Y(Cr50) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cr50) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cr50) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cr50) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cr50) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr50_to_V50_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cr50 --> V50

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cr50_V50_meta, j_Cr50_V50_rhoy, j_Cr50_V50_temp, j_Cr50_V50_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cr50_to_V50_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr50_to_V50_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cr50 and V50
            // this is the derivative of ∂Y(Cr50)/∂t = -Y(Cr50) λ and
            //                           ∂Y(V50)/∂t = +Y(Cr50) λ
            rate_eval.dweak_ydot_dYe(Cr50) -= rho * Y(Cr50) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(V50) += rho * Y(Cr50) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cr50) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cr50) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cr50) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cr50) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn50_to_Cr50_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Mn50 --> Cr50

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Mn50_Cr50_meta, j_Mn50_Cr50_rhoy, j_Mn50_Cr50_temp, j_Mn50_Cr50_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Mn50_to_Cr50_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn50_to_Cr50_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Mn50 and Cr50
            // this is the derivative of ∂Y(Mn50)/∂t = -Y(Mn50) λ and
            //                           ∂Y(Cr50)/∂t = +Y(Mn50) λ
            rate_eval.dweak_ydot_dYe(Mn50) -= rho * Y(Mn50) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cr50) += rho * Y(Mn50) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Mn50) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Mn50) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Mn50) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Mn50) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti50_to_V50_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ti50 --> V50

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ti50_V50_meta, j_Ti50_V50_rhoy, j_Ti50_V50_temp, j_Ti50_V50_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ti50_to_V50_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti50_to_V50_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ti50 and V50
            // this is the derivative of ∂Y(Ti50)/∂t = -Y(Ti50) λ and
            //                           ∂Y(V50)/∂t = +Y(Ti50) λ
            rate_eval.dweak_ydot_dYe(Ti50) -= rho * Y(Ti50) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(V50) += rho * Y(Ti50) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ti50) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ti50) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ti50) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ti50) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V50_to_Cr50_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // V50 --> Cr50

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_V50_Cr50_meta, j_V50_Cr50_rhoy, j_V50_Cr50_temp, j_V50_Cr50_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_V50_to_Cr50_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V50_to_Cr50_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to V50 and Cr50
            // this is the derivative of ∂Y(V50)/∂t = -Y(V50) λ and
            //                           ∂Y(Cr50)/∂t = +Y(V50) λ
            rate_eval.dweak_ydot_dYe(V50) -= rho * Y(V50) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cr50) += rho * Y(V50) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(V50) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(V50) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(V50) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(V50) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V50_to_Ti50_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // V50 --> Ti50

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_V50_Ti50_meta, j_V50_Ti50_rhoy, j_V50_Ti50_temp, j_V50_Ti50_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_V50_to_Ti50_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V50_to_Ti50_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to V50 and Ti50
            // this is the derivative of ∂Y(V50)/∂t = -Y(V50) λ and
            //                           ∂Y(Ti50)/∂t = +Y(V50) λ
            rate_eval.dweak_ydot_dYe(V50) -= rho * Y(V50) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ti50) += rho * Y(V50) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(V50) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(V50) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(V50) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(V50) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr51_to_Mn51_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cr51 --> Mn51

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cr51_Mn51_meta, j_Cr51_Mn51_rhoy, j_Cr51_Mn51_temp, j_Cr51_Mn51_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cr51_to_Mn51_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr51_to_Mn51_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cr51 and Mn51
            // this is the derivative of ∂Y(Cr51)/∂t = -Y(Cr51) λ and
            //                           ∂Y(Mn51)/∂t = +Y(Cr51) λ
            rate_eval.dweak_ydot_dYe(Cr51) -= rho * Y(Cr51) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Mn51) += rho * Y(Cr51) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cr51) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cr51) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cr51) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cr51) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr51_to_V51_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cr51 --> V51

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cr51_V51_meta, j_Cr51_V51_rhoy, j_Cr51_V51_temp, j_Cr51_V51_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cr51_to_V51_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr51_to_V51_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cr51 and V51
            // this is the derivative of ∂Y(Cr51)/∂t = -Y(Cr51) λ and
            //                           ∂Y(V51)/∂t = +Y(Cr51) λ
            rate_eval.dweak_ydot_dYe(Cr51) -= rho * Y(Cr51) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(V51) += rho * Y(Cr51) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cr51) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cr51) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cr51) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cr51) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn51_to_Cr51_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Mn51 --> Cr51

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Mn51_Cr51_meta, j_Mn51_Cr51_rhoy, j_Mn51_Cr51_temp, j_Mn51_Cr51_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Mn51_to_Cr51_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn51_to_Cr51_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Mn51 and Cr51
            // this is the derivative of ∂Y(Mn51)/∂t = -Y(Mn51) λ and
            //                           ∂Y(Cr51)/∂t = +Y(Mn51) λ
            rate_eval.dweak_ydot_dYe(Mn51) -= rho * Y(Mn51) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cr51) += rho * Y(Mn51) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Mn51) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Mn51) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Mn51) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Mn51) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti51_to_V51_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ti51 --> V51

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ti51_V51_meta, j_Ti51_V51_rhoy, j_Ti51_V51_temp, j_Ti51_V51_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ti51_to_V51_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti51_to_V51_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ti51 and V51
            // this is the derivative of ∂Y(Ti51)/∂t = -Y(Ti51) λ and
            //                           ∂Y(V51)/∂t = +Y(Ti51) λ
            rate_eval.dweak_ydot_dYe(Ti51) -= rho * Y(Ti51) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(V51) += rho * Y(Ti51) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ti51) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ti51) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ti51) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ti51) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V51_to_Cr51_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // V51 --> Cr51

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_V51_Cr51_meta, j_V51_Cr51_rhoy, j_V51_Cr51_temp, j_V51_Cr51_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_V51_to_Cr51_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V51_to_Cr51_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to V51 and Cr51
            // this is the derivative of ∂Y(V51)/∂t = -Y(V51) λ and
            //                           ∂Y(Cr51)/∂t = +Y(V51) λ
            rate_eval.dweak_ydot_dYe(V51) -= rho * Y(V51) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cr51) += rho * Y(V51) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(V51) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(V51) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(V51) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(V51) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V51_to_Ti51_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // V51 --> Ti51

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_V51_Ti51_meta, j_V51_Ti51_rhoy, j_V51_Ti51_temp, j_V51_Ti51_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_V51_to_Ti51_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V51_to_Ti51_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to V51 and Ti51
            // this is the derivative of ∂Y(V51)/∂t = -Y(V51) λ and
            //                           ∂Y(Ti51)/∂t = +Y(V51) λ
            rate_eval.dweak_ydot_dYe(V51) -= rho * Y(V51) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ti51) += rho * Y(V51) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(V51) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(V51) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(V51) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(V51) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr52_to_Mn52_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cr52 --> Mn52

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cr52_Mn52_meta, j_Cr52_Mn52_rhoy, j_Cr52_Mn52_temp, j_Cr52_Mn52_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cr52_to_Mn52_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr52_to_Mn52_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cr52 and Mn52
            // this is the derivative of ∂Y(Cr52)/∂t = -Y(Cr52) λ and
            //                           ∂Y(Mn52)/∂t = +Y(Cr52) λ
            rate_eval.dweak_ydot_dYe(Cr52) -= rho * Y(Cr52) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Mn52) += rho * Y(Cr52) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cr52) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cr52) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cr52) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cr52) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr52_to_V52_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cr52 --> V52

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cr52_V52_meta, j_Cr52_V52_rhoy, j_Cr52_V52_temp, j_Cr52_V52_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cr52_to_V52_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr52_to_V52_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cr52 and V52
            // this is the derivative of ∂Y(Cr52)/∂t = -Y(Cr52) λ and
            //                           ∂Y(V52)/∂t = +Y(Cr52) λ
            rate_eval.dweak_ydot_dYe(Cr52) -= rho * Y(Cr52) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(V52) += rho * Y(Cr52) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cr52) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cr52) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cr52) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cr52) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe52_to_Mn52_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Fe52 --> Mn52

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Fe52_Mn52_meta, j_Fe52_Mn52_rhoy, j_Fe52_Mn52_temp, j_Fe52_Mn52_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Fe52_to_Mn52_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe52_to_Mn52_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Fe52 and Mn52
            // this is the derivative of ∂Y(Fe52)/∂t = -Y(Fe52) λ and
            //                           ∂Y(Mn52)/∂t = +Y(Fe52) λ
            rate_eval.dweak_ydot_dYe(Fe52) -= rho * Y(Fe52) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Mn52) += rho * Y(Fe52) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Fe52) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Fe52) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Fe52) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Fe52) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn52_to_Cr52_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Mn52 --> Cr52

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Mn52_Cr52_meta, j_Mn52_Cr52_rhoy, j_Mn52_Cr52_temp, j_Mn52_Cr52_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Mn52_to_Cr52_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn52_to_Cr52_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Mn52 and Cr52
            // this is the derivative of ∂Y(Mn52)/∂t = -Y(Mn52) λ and
            //                           ∂Y(Cr52)/∂t = +Y(Mn52) λ
            rate_eval.dweak_ydot_dYe(Mn52) -= rho * Y(Mn52) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cr52) += rho * Y(Mn52) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Mn52) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Mn52) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Mn52) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Mn52) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn52_to_Fe52_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Mn52 --> Fe52

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Mn52_Fe52_meta, j_Mn52_Fe52_rhoy, j_Mn52_Fe52_temp, j_Mn52_Fe52_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Mn52_to_Fe52_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn52_to_Fe52_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Mn52 and Fe52
            // this is the derivative of ∂Y(Mn52)/∂t = -Y(Mn52) λ and
            //                           ∂Y(Fe52)/∂t = +Y(Mn52) λ
            rate_eval.dweak_ydot_dYe(Mn52) -= rho * Y(Mn52) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Fe52) += rho * Y(Mn52) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Mn52) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Mn52) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Mn52) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Mn52) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V52_to_Cr52_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // V52 --> Cr52

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_V52_Cr52_meta, j_V52_Cr52_rhoy, j_V52_Cr52_temp, j_V52_Cr52_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_V52_to_Cr52_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V52_to_Cr52_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to V52 and Cr52
            // this is the derivative of ∂Y(V52)/∂t = -Y(V52) λ and
            //                           ∂Y(Cr52)/∂t = +Y(V52) λ
            rate_eval.dweak_ydot_dYe(V52) -= rho * Y(V52) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cr52) += rho * Y(V52) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(V52) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(V52) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(V52) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(V52) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr53_to_Mn53_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cr53 --> Mn53

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cr53_Mn53_meta, j_Cr53_Mn53_rhoy, j_Cr53_Mn53_temp, j_Cr53_Mn53_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cr53_to_Mn53_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr53_to_Mn53_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cr53 and Mn53
            // this is the derivative of ∂Y(Cr53)/∂t = -Y(Cr53) λ and
            //                           ∂Y(Mn53)/∂t = +Y(Cr53) λ
            rate_eval.dweak_ydot_dYe(Cr53) -= rho * Y(Cr53) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Mn53) += rho * Y(Cr53) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cr53) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cr53) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cr53) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cr53) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe53_to_Mn53_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Fe53 --> Mn53

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Fe53_Mn53_meta, j_Fe53_Mn53_rhoy, j_Fe53_Mn53_temp, j_Fe53_Mn53_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Fe53_to_Mn53_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe53_to_Mn53_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Fe53 and Mn53
            // this is the derivative of ∂Y(Fe53)/∂t = -Y(Fe53) λ and
            //                           ∂Y(Mn53)/∂t = +Y(Fe53) λ
            rate_eval.dweak_ydot_dYe(Fe53) -= rho * Y(Fe53) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Mn53) += rho * Y(Fe53) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Fe53) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Fe53) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Fe53) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Fe53) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn53_to_Cr53_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Mn53 --> Cr53

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Mn53_Cr53_meta, j_Mn53_Cr53_rhoy, j_Mn53_Cr53_temp, j_Mn53_Cr53_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Mn53_to_Cr53_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn53_to_Cr53_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Mn53 and Cr53
            // this is the derivative of ∂Y(Mn53)/∂t = -Y(Mn53) λ and
            //                           ∂Y(Cr53)/∂t = +Y(Mn53) λ
            rate_eval.dweak_ydot_dYe(Mn53) -= rho * Y(Mn53) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cr53) += rho * Y(Mn53) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Mn53) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Mn53) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Mn53) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Mn53) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn53_to_Fe53_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Mn53 --> Fe53

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Mn53_Fe53_meta, j_Mn53_Fe53_rhoy, j_Mn53_Fe53_temp, j_Mn53_Fe53_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Mn53_to_Fe53_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn53_to_Fe53_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Mn53 and Fe53
            // this is the derivative of ∂Y(Mn53)/∂t = -Y(Mn53) λ and
            //                           ∂Y(Fe53)/∂t = +Y(Mn53) λ
            rate_eval.dweak_ydot_dYe(Mn53) -= rho * Y(Mn53) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Fe53) += rho * Y(Mn53) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Mn53) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Mn53) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Mn53) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Mn53) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co54_to_Fe54_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Co54 --> Fe54

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Co54_Fe54_meta, j_Co54_Fe54_rhoy, j_Co54_Fe54_temp, j_Co54_Fe54_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Co54_to_Fe54_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co54_to_Fe54_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Co54 and Fe54
            // this is the derivative of ∂Y(Co54)/∂t = -Y(Co54) λ and
            //                           ∂Y(Fe54)/∂t = +Y(Co54) λ
            rate_eval.dweak_ydot_dYe(Co54) -= rho * Y(Co54) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Fe54) += rho * Y(Co54) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Co54) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Co54) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Co54) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Co54) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr54_to_Mn54_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cr54 --> Mn54

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cr54_Mn54_meta, j_Cr54_Mn54_rhoy, j_Cr54_Mn54_temp, j_Cr54_Mn54_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cr54_to_Mn54_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr54_to_Mn54_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cr54 and Mn54
            // this is the derivative of ∂Y(Cr54)/∂t = -Y(Cr54) λ and
            //                           ∂Y(Mn54)/∂t = +Y(Cr54) λ
            rate_eval.dweak_ydot_dYe(Cr54) -= rho * Y(Cr54) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Mn54) += rho * Y(Cr54) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cr54) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cr54) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cr54) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cr54) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe54_to_Co54_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Fe54 --> Co54

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Fe54_Co54_meta, j_Fe54_Co54_rhoy, j_Fe54_Co54_temp, j_Fe54_Co54_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Fe54_to_Co54_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe54_to_Co54_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Fe54 and Co54
            // this is the derivative of ∂Y(Fe54)/∂t = -Y(Fe54) λ and
            //                           ∂Y(Co54)/∂t = +Y(Fe54) λ
            rate_eval.dweak_ydot_dYe(Fe54) -= rho * Y(Fe54) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Co54) += rho * Y(Fe54) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Fe54) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Fe54) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Fe54) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Fe54) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe54_to_Mn54_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Fe54 --> Mn54

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Fe54_Mn54_meta, j_Fe54_Mn54_rhoy, j_Fe54_Mn54_temp, j_Fe54_Mn54_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Fe54_to_Mn54_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe54_to_Mn54_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Fe54 and Mn54
            // this is the derivative of ∂Y(Fe54)/∂t = -Y(Fe54) λ and
            //                           ∂Y(Mn54)/∂t = +Y(Fe54) λ
            rate_eval.dweak_ydot_dYe(Fe54) -= rho * Y(Fe54) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Mn54) += rho * Y(Fe54) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Fe54) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Fe54) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Fe54) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Fe54) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn54_to_Cr54_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Mn54 --> Cr54

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Mn54_Cr54_meta, j_Mn54_Cr54_rhoy, j_Mn54_Cr54_temp, j_Mn54_Cr54_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Mn54_to_Cr54_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn54_to_Cr54_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Mn54 and Cr54
            // this is the derivative of ∂Y(Mn54)/∂t = -Y(Mn54) λ and
            //                           ∂Y(Cr54)/∂t = +Y(Mn54) λ
            rate_eval.dweak_ydot_dYe(Mn54) -= rho * Y(Mn54) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cr54) += rho * Y(Mn54) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Mn54) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Mn54) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Mn54) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Mn54) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn54_to_Fe54_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Mn54 --> Fe54

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Mn54_Fe54_meta, j_Mn54_Fe54_rhoy, j_Mn54_Fe54_temp, j_Mn54_Fe54_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Mn54_to_Fe54_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn54_to_Fe54_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Mn54 and Fe54
            // this is the derivative of ∂Y(Mn54)/∂t = -Y(Mn54) λ and
            //                           ∂Y(Fe54)/∂t = +Y(Mn54) λ
            rate_eval.dweak_ydot_dYe(Mn54) -= rho * Y(Mn54) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Fe54) += rho * Y(Mn54) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Mn54) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Mn54) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Mn54) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Mn54) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co55_to_Fe55_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Co55 --> Fe55

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Co55_Fe55_meta, j_Co55_Fe55_rhoy, j_Co55_Fe55_temp, j_Co55_Fe55_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Co55_to_Fe55_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co55_to_Fe55_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Co55 and Fe55
            // this is the derivative of ∂Y(Co55)/∂t = -Y(Co55) λ and
            //                           ∂Y(Fe55)/∂t = +Y(Co55) λ
            rate_eval.dweak_ydot_dYe(Co55) -= rho * Y(Co55) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Fe55) += rho * Y(Co55) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Co55) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Co55) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Co55) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Co55) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe55_to_Co55_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Fe55 --> Co55

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Fe55_Co55_meta, j_Fe55_Co55_rhoy, j_Fe55_Co55_temp, j_Fe55_Co55_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Fe55_to_Co55_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe55_to_Co55_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Fe55 and Co55
            // this is the derivative of ∂Y(Fe55)/∂t = -Y(Fe55) λ and
            //                           ∂Y(Co55)/∂t = +Y(Fe55) λ
            rate_eval.dweak_ydot_dYe(Fe55) -= rho * Y(Fe55) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Co55) += rho * Y(Fe55) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Fe55) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Fe55) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Fe55) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Fe55) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe55_to_Mn55_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Fe55 --> Mn55

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Fe55_Mn55_meta, j_Fe55_Mn55_rhoy, j_Fe55_Mn55_temp, j_Fe55_Mn55_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Fe55_to_Mn55_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe55_to_Mn55_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Fe55 and Mn55
            // this is the derivative of ∂Y(Fe55)/∂t = -Y(Fe55) λ and
            //                           ∂Y(Mn55)/∂t = +Y(Fe55) λ
            rate_eval.dweak_ydot_dYe(Fe55) -= rho * Y(Fe55) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Mn55) += rho * Y(Fe55) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Fe55) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Fe55) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Fe55) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Fe55) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn55_to_Fe55_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Mn55 --> Fe55

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Mn55_Fe55_meta, j_Mn55_Fe55_rhoy, j_Mn55_Fe55_temp, j_Mn55_Fe55_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Mn55_to_Fe55_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn55_to_Fe55_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Mn55 and Fe55
            // this is the derivative of ∂Y(Mn55)/∂t = -Y(Mn55) λ and
            //                           ∂Y(Fe55)/∂t = +Y(Mn55) λ
            rate_eval.dweak_ydot_dYe(Mn55) -= rho * Y(Mn55) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Fe55) += rho * Y(Mn55) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Mn55) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Mn55) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Mn55) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Mn55) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co56_to_Fe56_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Co56 --> Fe56

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Co56_Fe56_meta, j_Co56_Fe56_rhoy, j_Co56_Fe56_temp, j_Co56_Fe56_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Co56_to_Fe56_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co56_to_Fe56_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Co56 and Fe56
            // this is the derivative of ∂Y(Co56)/∂t = -Y(Co56) λ and
            //                           ∂Y(Fe56)/∂t = +Y(Co56) λ
            rate_eval.dweak_ydot_dYe(Co56) -= rho * Y(Co56) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Fe56) += rho * Y(Co56) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Co56) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Co56) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Co56) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Co56) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co56_to_Ni56_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Co56 --> Ni56

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Co56_Ni56_meta, j_Co56_Ni56_rhoy, j_Co56_Ni56_temp, j_Co56_Ni56_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Co56_to_Ni56_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co56_to_Ni56_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Co56 and Ni56
            // this is the derivative of ∂Y(Co56)/∂t = -Y(Co56) λ and
            //                           ∂Y(Ni56)/∂t = +Y(Co56) λ
            rate_eval.dweak_ydot_dYe(Co56) -= rho * Y(Co56) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ni56) += rho * Y(Co56) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Co56) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Co56) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Co56) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Co56) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe56_to_Co56_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Fe56 --> Co56

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Fe56_Co56_meta, j_Fe56_Co56_rhoy, j_Fe56_Co56_temp, j_Fe56_Co56_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Fe56_to_Co56_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe56_to_Co56_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Fe56 and Co56
            // this is the derivative of ∂Y(Fe56)/∂t = -Y(Fe56) λ and
            //                           ∂Y(Co56)/∂t = +Y(Fe56) λ
            rate_eval.dweak_ydot_dYe(Fe56) -= rho * Y(Fe56) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Co56) += rho * Y(Fe56) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Fe56) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Fe56) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Fe56) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Fe56) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni56_to_Co56_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ni56 --> Co56

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ni56_Co56_meta, j_Ni56_Co56_rhoy, j_Ni56_Co56_temp, j_Ni56_Co56_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ni56_to_Co56_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni56_to_Co56_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ni56 and Co56
            // this is the derivative of ∂Y(Ni56)/∂t = -Y(Ni56) λ and
            //                           ∂Y(Co56)/∂t = +Y(Ni56) λ
            rate_eval.dweak_ydot_dYe(Ni56) -= rho * Y(Ni56) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Co56) += rho * Y(Ni56) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ni56) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ni56) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ni56) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ni56) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co57_to_Fe57_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Co57 --> Fe57

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Co57_Fe57_meta, j_Co57_Fe57_rhoy, j_Co57_Fe57_temp, j_Co57_Fe57_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Co57_to_Fe57_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co57_to_Fe57_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Co57 and Fe57
            // this is the derivative of ∂Y(Co57)/∂t = -Y(Co57) λ and
            //                           ∂Y(Fe57)/∂t = +Y(Co57) λ
            rate_eval.dweak_ydot_dYe(Co57) -= rho * Y(Co57) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Fe57) += rho * Y(Co57) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Co57) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Co57) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Co57) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Co57) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co57_to_Ni57_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Co57 --> Ni57

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Co57_Ni57_meta, j_Co57_Ni57_rhoy, j_Co57_Ni57_temp, j_Co57_Ni57_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Co57_to_Ni57_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co57_to_Ni57_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Co57 and Ni57
            // this is the derivative of ∂Y(Co57)/∂t = -Y(Co57) λ and
            //                           ∂Y(Ni57)/∂t = +Y(Co57) λ
            rate_eval.dweak_ydot_dYe(Co57) -= rho * Y(Co57) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ni57) += rho * Y(Co57) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Co57) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Co57) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Co57) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Co57) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe57_to_Co57_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Fe57 --> Co57

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Fe57_Co57_meta, j_Fe57_Co57_rhoy, j_Fe57_Co57_temp, j_Fe57_Co57_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Fe57_to_Co57_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe57_to_Co57_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Fe57 and Co57
            // this is the derivative of ∂Y(Fe57)/∂t = -Y(Fe57) λ and
            //                           ∂Y(Co57)/∂t = +Y(Fe57) λ
            rate_eval.dweak_ydot_dYe(Fe57) -= rho * Y(Fe57) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Co57) += rho * Y(Fe57) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Fe57) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Fe57) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Fe57) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Fe57) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni57_to_Co57_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ni57 --> Co57

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ni57_Co57_meta, j_Ni57_Co57_rhoy, j_Ni57_Co57_temp, j_Ni57_Co57_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ni57_to_Co57_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni57_to_Co57_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ni57 and Co57
            // this is the derivative of ∂Y(Ni57)/∂t = -Y(Ni57) λ and
            //                           ∂Y(Co57)/∂t = +Y(Ni57) λ
            rate_eval.dweak_ydot_dYe(Ni57) -= rho * Y(Ni57) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Co57) += rho * Y(Ni57) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ni57) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ni57) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ni57) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ni57) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co58_to_Fe58_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Co58 --> Fe58

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Co58_Fe58_meta, j_Co58_Fe58_rhoy, j_Co58_Fe58_temp, j_Co58_Fe58_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Co58_to_Fe58_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co58_to_Fe58_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Co58 and Fe58
            // this is the derivative of ∂Y(Co58)/∂t = -Y(Co58) λ and
            //                           ∂Y(Fe58)/∂t = +Y(Co58) λ
            rate_eval.dweak_ydot_dYe(Co58) -= rho * Y(Co58) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Fe58) += rho * Y(Co58) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Co58) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Co58) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Co58) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Co58) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co58_to_Ni58_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Co58 --> Ni58

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Co58_Ni58_meta, j_Co58_Ni58_rhoy, j_Co58_Ni58_temp, j_Co58_Ni58_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Co58_to_Ni58_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co58_to_Ni58_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Co58 and Ni58
            // this is the derivative of ∂Y(Co58)/∂t = -Y(Co58) λ and
            //                           ∂Y(Ni58)/∂t = +Y(Co58) λ
            rate_eval.dweak_ydot_dYe(Co58) -= rho * Y(Co58) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ni58) += rho * Y(Co58) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Co58) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Co58) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Co58) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Co58) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu58_to_Ni58_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu58 --> Ni58

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu58_Ni58_meta, j_Cu58_Ni58_rhoy, j_Cu58_Ni58_temp, j_Cu58_Ni58_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu58_to_Ni58_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu58_to_Ni58_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cu58 and Ni58
            // this is the derivative of ∂Y(Cu58)/∂t = -Y(Cu58) λ and
            //                           ∂Y(Ni58)/∂t = +Y(Cu58) λ
            rate_eval.dweak_ydot_dYe(Cu58) -= rho * Y(Cu58) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ni58) += rho * Y(Cu58) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cu58) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cu58) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cu58) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu58) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe58_to_Co58_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Fe58 --> Co58

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Fe58_Co58_meta, j_Fe58_Co58_rhoy, j_Fe58_Co58_temp, j_Fe58_Co58_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Fe58_to_Co58_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe58_to_Co58_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Fe58 and Co58
            // this is the derivative of ∂Y(Fe58)/∂t = -Y(Fe58) λ and
            //                           ∂Y(Co58)/∂t = +Y(Fe58) λ
            rate_eval.dweak_ydot_dYe(Fe58) -= rho * Y(Fe58) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Co58) += rho * Y(Fe58) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Fe58) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Fe58) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Fe58) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Fe58) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni58_to_Co58_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ni58 --> Co58

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ni58_Co58_meta, j_Ni58_Co58_rhoy, j_Ni58_Co58_temp, j_Ni58_Co58_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ni58_to_Co58_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni58_to_Co58_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ni58 and Co58
            // this is the derivative of ∂Y(Ni58)/∂t = -Y(Ni58) λ and
            //                           ∂Y(Co58)/∂t = +Y(Ni58) λ
            rate_eval.dweak_ydot_dYe(Ni58) -= rho * Y(Ni58) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Co58) += rho * Y(Ni58) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ni58) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ni58) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ni58) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ni58) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni58_to_Cu58_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ni58 --> Cu58

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ni58_Cu58_meta, j_Ni58_Cu58_rhoy, j_Ni58_Cu58_temp, j_Ni58_Cu58_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ni58_to_Cu58_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni58_to_Cu58_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ni58 and Cu58
            // this is the derivative of ∂Y(Ni58)/∂t = -Y(Ni58) λ and
            //                           ∂Y(Cu58)/∂t = +Y(Ni58) λ
            rate_eval.dweak_ydot_dYe(Ni58) -= rho * Y(Ni58) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cu58) += rho * Y(Ni58) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ni58) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ni58) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ni58) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ni58) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co59_to_Ni59_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Co59 --> Ni59

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Co59_Ni59_meta, j_Co59_Ni59_rhoy, j_Co59_Ni59_temp, j_Co59_Ni59_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Co59_to_Ni59_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co59_to_Ni59_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Co59 and Ni59
            // this is the derivative of ∂Y(Co59)/∂t = -Y(Co59) λ and
            //                           ∂Y(Ni59)/∂t = +Y(Co59) λ
            rate_eval.dweak_ydot_dYe(Co59) -= rho * Y(Co59) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ni59) += rho * Y(Co59) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Co59) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Co59) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Co59) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Co59) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu59_to_Ni59_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu59 --> Ni59

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu59_Ni59_meta, j_Cu59_Ni59_rhoy, j_Cu59_Ni59_temp, j_Cu59_Ni59_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu59_to_Ni59_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu59_to_Ni59_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cu59 and Ni59
            // this is the derivative of ∂Y(Cu59)/∂t = -Y(Cu59) λ and
            //                           ∂Y(Ni59)/∂t = +Y(Cu59) λ
            rate_eval.dweak_ydot_dYe(Cu59) -= rho * Y(Cu59) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ni59) += rho * Y(Cu59) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cu59) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cu59) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cu59) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu59) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni59_to_Co59_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ni59 --> Co59

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ni59_Co59_meta, j_Ni59_Co59_rhoy, j_Ni59_Co59_temp, j_Ni59_Co59_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ni59_to_Co59_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni59_to_Co59_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ni59 and Co59
            // this is the derivative of ∂Y(Ni59)/∂t = -Y(Ni59) λ and
            //                           ∂Y(Co59)/∂t = +Y(Ni59) λ
            rate_eval.dweak_ydot_dYe(Ni59) -= rho * Y(Ni59) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Co59) += rho * Y(Ni59) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ni59) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ni59) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ni59) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ni59) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni59_to_Cu59_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ni59 --> Cu59

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ni59_Cu59_meta, j_Ni59_Cu59_rhoy, j_Ni59_Cu59_temp, j_Ni59_Cu59_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ni59_to_Cu59_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni59_to_Cu59_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ni59 and Cu59
            // this is the derivative of ∂Y(Ni59)/∂t = -Y(Ni59) λ and
            //                           ∂Y(Cu59)/∂t = +Y(Ni59) λ
            rate_eval.dweak_ydot_dYe(Ni59) -= rho * Y(Ni59) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cu59) += rho * Y(Ni59) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ni59) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ni59) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ni59) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ni59) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu60_to_Ni60_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu60 --> Ni60

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu60_Ni60_meta, j_Cu60_Ni60_rhoy, j_Cu60_Ni60_temp, j_Cu60_Ni60_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu60_to_Ni60_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu60_to_Ni60_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cu60 and Ni60
            // this is the derivative of ∂Y(Cu60)/∂t = -Y(Cu60) λ and
            //                           ∂Y(Ni60)/∂t = +Y(Cu60) λ
            rate_eval.dweak_ydot_dYe(Cu60) -= rho * Y(Cu60) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ni60) += rho * Y(Cu60) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cu60) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cu60) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cu60) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu60) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu60_to_Zn60_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu60 --> Zn60

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu60_Zn60_meta, j_Cu60_Zn60_rhoy, j_Cu60_Zn60_temp, j_Cu60_Zn60_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu60_to_Zn60_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu60_to_Zn60_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cu60 and Zn60
            // this is the derivative of ∂Y(Cu60)/∂t = -Y(Cu60) λ and
            //                           ∂Y(Zn60)/∂t = +Y(Cu60) λ
            rate_eval.dweak_ydot_dYe(Cu60) -= rho * Y(Cu60) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Zn60) += rho * Y(Cu60) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cu60) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cu60) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cu60) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu60) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni60_to_Cu60_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ni60 --> Cu60

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ni60_Cu60_meta, j_Ni60_Cu60_rhoy, j_Ni60_Cu60_temp, j_Ni60_Cu60_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ni60_to_Cu60_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni60_to_Cu60_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ni60 and Cu60
            // this is the derivative of ∂Y(Ni60)/∂t = -Y(Ni60) λ and
            //                           ∂Y(Cu60)/∂t = +Y(Ni60) λ
            rate_eval.dweak_ydot_dYe(Ni60) -= rho * Y(Ni60) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cu60) += rho * Y(Ni60) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ni60) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ni60) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ni60) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ni60) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Zn60_to_Cu60_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Zn60 --> Cu60

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Zn60_Cu60_meta, j_Zn60_Cu60_rhoy, j_Zn60_Cu60_temp, j_Zn60_Cu60_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Zn60_to_Cu60_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn60_to_Cu60_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Zn60 and Cu60
            // this is the derivative of ∂Y(Zn60)/∂t = -Y(Zn60) λ and
            //                           ∂Y(Cu60)/∂t = +Y(Zn60) λ
            rate_eval.dweak_ydot_dYe(Zn60) -= rho * Y(Zn60) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cu60) += rho * Y(Zn60) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Zn60) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Zn60) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Zn60) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Zn60) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu61_to_Ni61_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu61 --> Ni61

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu61_Ni61_meta, j_Cu61_Ni61_rhoy, j_Cu61_Ni61_temp, j_Cu61_Ni61_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu61_to_Ni61_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu61_to_Ni61_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cu61 and Ni61
            // this is the derivative of ∂Y(Cu61)/∂t = -Y(Cu61) λ and
            //                           ∂Y(Ni61)/∂t = +Y(Cu61) λ
            rate_eval.dweak_ydot_dYe(Cu61) -= rho * Y(Cu61) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ni61) += rho * Y(Cu61) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cu61) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cu61) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cu61) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu61) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu61_to_Zn61_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu61 --> Zn61

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu61_Zn61_meta, j_Cu61_Zn61_rhoy, j_Cu61_Zn61_temp, j_Cu61_Zn61_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu61_to_Zn61_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu61_to_Zn61_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cu61 and Zn61
            // this is the derivative of ∂Y(Cu61)/∂t = -Y(Cu61) λ and
            //                           ∂Y(Zn61)/∂t = +Y(Cu61) λ
            rate_eval.dweak_ydot_dYe(Cu61) -= rho * Y(Cu61) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Zn61) += rho * Y(Cu61) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cu61) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cu61) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cu61) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu61) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni61_to_Cu61_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ni61 --> Cu61

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ni61_Cu61_meta, j_Ni61_Cu61_rhoy, j_Ni61_Cu61_temp, j_Ni61_Cu61_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ni61_to_Cu61_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni61_to_Cu61_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ni61 and Cu61
            // this is the derivative of ∂Y(Ni61)/∂t = -Y(Ni61) λ and
            //                           ∂Y(Cu61)/∂t = +Y(Ni61) λ
            rate_eval.dweak_ydot_dYe(Ni61) -= rho * Y(Ni61) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cu61) += rho * Y(Ni61) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ni61) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ni61) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ni61) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ni61) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Zn61_to_Cu61_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Zn61 --> Cu61

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Zn61_Cu61_meta, j_Zn61_Cu61_rhoy, j_Zn61_Cu61_temp, j_Zn61_Cu61_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Zn61_to_Cu61_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn61_to_Cu61_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Zn61 and Cu61
            // this is the derivative of ∂Y(Zn61)/∂t = -Y(Zn61) λ and
            //                           ∂Y(Cu61)/∂t = +Y(Zn61) λ
            rate_eval.dweak_ydot_dYe(Zn61) -= rho * Y(Zn61) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cu61) += rho * Y(Zn61) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Zn61) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Zn61) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Zn61) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Zn61) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu62_to_Ni62_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu62 --> Ni62

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu62_Ni62_meta, j_Cu62_Ni62_rhoy, j_Cu62_Ni62_temp, j_Cu62_Ni62_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu62_to_Ni62_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu62_to_Ni62_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cu62 and Ni62
            // this is the derivative of ∂Y(Cu62)/∂t = -Y(Cu62) λ and
            //                           ∂Y(Ni62)/∂t = +Y(Cu62) λ
            rate_eval.dweak_ydot_dYe(Cu62) -= rho * Y(Cu62) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ni62) += rho * Y(Cu62) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cu62) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cu62) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cu62) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu62) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu62_to_Zn62_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu62 --> Zn62

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu62_Zn62_meta, j_Cu62_Zn62_rhoy, j_Cu62_Zn62_temp, j_Cu62_Zn62_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu62_to_Zn62_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu62_to_Zn62_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cu62 and Zn62
            // this is the derivative of ∂Y(Cu62)/∂t = -Y(Cu62) λ and
            //                           ∂Y(Zn62)/∂t = +Y(Cu62) λ
            rate_eval.dweak_ydot_dYe(Cu62) -= rho * Y(Cu62) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Zn62) += rho * Y(Cu62) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cu62) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cu62) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cu62) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu62) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ga62_to_Zn62_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ga62 --> Zn62

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ga62_Zn62_meta, j_Ga62_Zn62_rhoy, j_Ga62_Zn62_temp, j_Ga62_Zn62_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ga62_to_Zn62_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ga62_to_Zn62_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ga62 and Zn62
            // this is the derivative of ∂Y(Ga62)/∂t = -Y(Ga62) λ and
            //                           ∂Y(Zn62)/∂t = +Y(Ga62) λ
            rate_eval.dweak_ydot_dYe(Ga62) -= rho * Y(Ga62) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Zn62) += rho * Y(Ga62) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ga62) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ga62) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ga62) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ga62) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni62_to_Cu62_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ni62 --> Cu62

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ni62_Cu62_meta, j_Ni62_Cu62_rhoy, j_Ni62_Cu62_temp, j_Ni62_Cu62_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ni62_to_Cu62_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni62_to_Cu62_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ni62 and Cu62
            // this is the derivative of ∂Y(Ni62)/∂t = -Y(Ni62) λ and
            //                           ∂Y(Cu62)/∂t = +Y(Ni62) λ
            rate_eval.dweak_ydot_dYe(Ni62) -= rho * Y(Ni62) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cu62) += rho * Y(Ni62) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ni62) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ni62) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ni62) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ni62) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Zn62_to_Cu62_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Zn62 --> Cu62

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Zn62_Cu62_meta, j_Zn62_Cu62_rhoy, j_Zn62_Cu62_temp, j_Zn62_Cu62_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Zn62_to_Cu62_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn62_to_Cu62_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Zn62 and Cu62
            // this is the derivative of ∂Y(Zn62)/∂t = -Y(Zn62) λ and
            //                           ∂Y(Cu62)/∂t = +Y(Zn62) λ
            rate_eval.dweak_ydot_dYe(Zn62) -= rho * Y(Zn62) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cu62) += rho * Y(Zn62) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Zn62) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Zn62) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Zn62) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Zn62) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Zn62_to_Ga62_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Zn62 --> Ga62

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Zn62_Ga62_meta, j_Zn62_Ga62_rhoy, j_Zn62_Ga62_temp, j_Zn62_Ga62_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Zn62_to_Ga62_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn62_to_Ga62_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Zn62 and Ga62
            // this is the derivative of ∂Y(Zn62)/∂t = -Y(Zn62) λ and
            //                           ∂Y(Ga62)/∂t = +Y(Zn62) λ
            rate_eval.dweak_ydot_dYe(Zn62) -= rho * Y(Zn62) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ga62) += rho * Y(Zn62) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Zn62) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Zn62) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Zn62) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Zn62) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu63_to_Ni63_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu63 --> Ni63

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu63_Ni63_meta, j_Cu63_Ni63_rhoy, j_Cu63_Ni63_temp, j_Cu63_Ni63_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu63_to_Ni63_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu63_to_Ni63_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cu63 and Ni63
            // this is the derivative of ∂Y(Cu63)/∂t = -Y(Cu63) λ and
            //                           ∂Y(Ni63)/∂t = +Y(Cu63) λ
            rate_eval.dweak_ydot_dYe(Cu63) -= rho * Y(Cu63) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ni63) += rho * Y(Cu63) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cu63) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cu63) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cu63) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu63) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu63_to_Zn63_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu63 --> Zn63

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu63_Zn63_meta, j_Cu63_Zn63_rhoy, j_Cu63_Zn63_temp, j_Cu63_Zn63_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu63_to_Zn63_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu63_to_Zn63_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cu63 and Zn63
            // this is the derivative of ∂Y(Cu63)/∂t = -Y(Cu63) λ and
            //                           ∂Y(Zn63)/∂t = +Y(Cu63) λ
            rate_eval.dweak_ydot_dYe(Cu63) -= rho * Y(Cu63) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Zn63) += rho * Y(Cu63) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cu63) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cu63) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cu63) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu63) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ga63_to_Zn63_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ga63 --> Zn63

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ga63_Zn63_meta, j_Ga63_Zn63_rhoy, j_Ga63_Zn63_temp, j_Ga63_Zn63_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ga63_to_Zn63_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ga63_to_Zn63_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ga63 and Zn63
            // this is the derivative of ∂Y(Ga63)/∂t = -Y(Ga63) λ and
            //                           ∂Y(Zn63)/∂t = +Y(Ga63) λ
            rate_eval.dweak_ydot_dYe(Ga63) -= rho * Y(Ga63) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Zn63) += rho * Y(Ga63) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ga63) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ga63) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ga63) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ga63) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni63_to_Cu63_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ni63 --> Cu63

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ni63_Cu63_meta, j_Ni63_Cu63_rhoy, j_Ni63_Cu63_temp, j_Ni63_Cu63_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ni63_to_Cu63_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni63_to_Cu63_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ni63 and Cu63
            // this is the derivative of ∂Y(Ni63)/∂t = -Y(Ni63) λ and
            //                           ∂Y(Cu63)/∂t = +Y(Ni63) λ
            rate_eval.dweak_ydot_dYe(Ni63) -= rho * Y(Ni63) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cu63) += rho * Y(Ni63) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ni63) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ni63) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ni63) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ni63) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Zn63_to_Cu63_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Zn63 --> Cu63

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Zn63_Cu63_meta, j_Zn63_Cu63_rhoy, j_Zn63_Cu63_temp, j_Zn63_Cu63_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Zn63_to_Cu63_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn63_to_Cu63_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Zn63 and Cu63
            // this is the derivative of ∂Y(Zn63)/∂t = -Y(Zn63) λ and
            //                           ∂Y(Cu63)/∂t = +Y(Zn63) λ
            rate_eval.dweak_ydot_dYe(Zn63) -= rho * Y(Zn63) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cu63) += rho * Y(Zn63) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Zn63) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Zn63) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Zn63) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Zn63) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Zn63_to_Ga63_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Zn63 --> Ga63

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Zn63_Ga63_meta, j_Zn63_Ga63_rhoy, j_Zn63_Ga63_temp, j_Zn63_Ga63_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Zn63_to_Ga63_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn63_to_Ga63_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Zn63 and Ga63
            // this is the derivative of ∂Y(Zn63)/∂t = -Y(Zn63) λ and
            //                           ∂Y(Ga63)/∂t = +Y(Zn63) λ
            rate_eval.dweak_ydot_dYe(Zn63) -= rho * Y(Zn63) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ga63) += rho * Y(Zn63) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Zn63) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Zn63) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Zn63) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Zn63) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu64_to_Ni64_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu64 --> Ni64

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu64_Ni64_meta, j_Cu64_Ni64_rhoy, j_Cu64_Ni64_temp, j_Cu64_Ni64_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu64_to_Ni64_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu64_to_Ni64_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cu64 and Ni64
            // this is the derivative of ∂Y(Cu64)/∂t = -Y(Cu64) λ and
            //                           ∂Y(Ni64)/∂t = +Y(Cu64) λ
            rate_eval.dweak_ydot_dYe(Cu64) -= rho * Y(Cu64) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ni64) += rho * Y(Cu64) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cu64) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cu64) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cu64) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu64) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu64_to_Zn64_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu64 --> Zn64

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu64_Zn64_meta, j_Cu64_Zn64_rhoy, j_Cu64_Zn64_temp, j_Cu64_Zn64_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu64_to_Zn64_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu64_to_Zn64_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cu64 and Zn64
            // this is the derivative of ∂Y(Cu64)/∂t = -Y(Cu64) λ and
            //                           ∂Y(Zn64)/∂t = +Y(Cu64) λ
            rate_eval.dweak_ydot_dYe(Cu64) -= rho * Y(Cu64) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Zn64) += rho * Y(Cu64) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cu64) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cu64) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cu64) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu64) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ga64_to_Ge64_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ga64 --> Ge64

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ga64_Ge64_meta, j_Ga64_Ge64_rhoy, j_Ga64_Ge64_temp, j_Ga64_Ge64_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ga64_to_Ge64_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ga64_to_Ge64_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ga64 and Ge64
            // this is the derivative of ∂Y(Ga64)/∂t = -Y(Ga64) λ and
            //                           ∂Y(Ge64)/∂t = +Y(Ga64) λ
            rate_eval.dweak_ydot_dYe(Ga64) -= rho * Y(Ga64) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ge64) += rho * Y(Ga64) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ga64) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ga64) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ga64) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ga64) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ga64_to_Zn64_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ga64 --> Zn64

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ga64_Zn64_meta, j_Ga64_Zn64_rhoy, j_Ga64_Zn64_temp, j_Ga64_Zn64_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ga64_to_Zn64_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ga64_to_Zn64_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ga64 and Zn64
            // this is the derivative of ∂Y(Ga64)/∂t = -Y(Ga64) λ and
            //                           ∂Y(Zn64)/∂t = +Y(Ga64) λ
            rate_eval.dweak_ydot_dYe(Ga64) -= rho * Y(Ga64) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Zn64) += rho * Y(Ga64) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ga64) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ga64) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ga64) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ga64) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ge64_to_Ga64_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ge64 --> Ga64

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ge64_Ga64_meta, j_Ge64_Ga64_rhoy, j_Ge64_Ga64_temp, j_Ge64_Ga64_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ge64_to_Ga64_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ge64_to_Ga64_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ge64 and Ga64
            // this is the derivative of ∂Y(Ge64)/∂t = -Y(Ge64) λ and
            //                           ∂Y(Ga64)/∂t = +Y(Ge64) λ
            rate_eval.dweak_ydot_dYe(Ge64) -= rho * Y(Ge64) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ga64) += rho * Y(Ge64) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ge64) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ge64) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ge64) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ge64) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni64_to_Cu64_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Ni64 --> Cu64

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Ni64_Cu64_meta, j_Ni64_Cu64_rhoy, j_Ni64_Cu64_temp, j_Ni64_Cu64_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Ni64_to_Cu64_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni64_to_Cu64_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Ni64 and Cu64
            // this is the derivative of ∂Y(Ni64)/∂t = -Y(Ni64) λ and
            //                           ∂Y(Cu64)/∂t = +Y(Ni64) λ
            rate_eval.dweak_ydot_dYe(Ni64) -= rho * Y(Ni64) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cu64) += rho * Y(Ni64) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Ni64) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Ni64) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Ni64) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Ni64) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Zn64_to_Cu64_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Zn64 --> Cu64

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Zn64_Cu64_meta, j_Zn64_Cu64_rhoy, j_Zn64_Cu64_temp, j_Zn64_Cu64_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Zn64_to_Cu64_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn64_to_Cu64_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Zn64 and Cu64
            // this is the derivative of ∂Y(Zn64)/∂t = -Y(Zn64) λ and
            //                           ∂Y(Cu64)/∂t = +Y(Zn64) λ
            rate_eval.dweak_ydot_dYe(Zn64) -= rho * Y(Zn64) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cu64) += rho * Y(Zn64) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Zn64) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Zn64) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Zn64) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Zn64) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Zn64_to_Ga64_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Zn64 --> Ga64

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Zn64_Ga64_meta, j_Zn64_Ga64_rhoy, j_Zn64_Ga64_temp, j_Zn64_Ga64_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Zn64_to_Ga64_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn64_to_Ga64_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Zn64 and Ga64
            // this is the derivative of ∂Y(Zn64)/∂t = -Y(Zn64) λ and
            //                           ∂Y(Ga64)/∂t = +Y(Zn64) λ
            rate_eval.dweak_ydot_dYe(Zn64) -= rho * Y(Zn64) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Ga64) += rho * Y(Zn64) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Zn64) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Zn64) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Zn64) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Zn64) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cu65_to_Zn65_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Cu65 --> Zn65

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Cu65_Zn65_meta, j_Cu65_Zn65_rhoy, j_Cu65_Zn65_temp, j_Cu65_Zn65_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Cu65_to_Zn65_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cu65_to_Zn65_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Cu65 and Zn65
            // this is the derivative of ∂Y(Cu65)/∂t = -Y(Cu65) λ and
            //                           ∂Y(Zn65)/∂t = +Y(Cu65) λ
            rate_eval.dweak_ydot_dYe(Cu65) -= rho * Y(Cu65) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Zn65) += rho * Y(Cu65) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Cu65) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Cu65) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Cu65) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Cu65) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Zn65_to_Cu65_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Zn65 --> Cu65

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Zn65_Cu65_meta, j_Zn65_Cu65_rhoy, j_Zn65_Cu65_temp, j_Zn65_Cu65_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Zn65_to_Cu65_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Zn65_to_Cu65_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Zn65 and Cu65
            // this is the derivative of ∂Y(Zn65)/∂t = -Y(Zn65) λ and
            //                           ∂Y(Cu65)/∂t = +Y(Zn65) λ
            rate_eval.dweak_ydot_dYe(Zn65) -= rho * Y(Zn65) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Cu65) += rho * Y(Zn65) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Zn65) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Zn65) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Zn65) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Zn65) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_n_to_p_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // n --> p

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_n_p_meta, j_n_p_rhoy, j_n_p_temp, j_n_p_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_n_to_p_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_n_to_p_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to n and p
            // this is the derivative of ∂Y(n)/∂t = -Y(n) λ and
            //                           ∂Y(p)/∂t = +Y(n) λ
            rate_eval.dweak_ydot_dYe(N) -= rho * Y(N) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(H1) += rho * Y(N) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(N) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(N) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(N) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(N) * (table_values.enu + table_values.gamma);
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_p_to_n_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // p --> n

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_p_n_meta, j_p_n_rhoy, j_p_n_temp, j_p_n_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_p_to_n_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_p_to_n_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to p and n
            // this is the derivative of ∂Y(p)/∂t = -Y(p) λ and
            //                           ∂Y(n)/∂t = +Y(p) λ
            rate_eval.dweak_ydot_dYe(H1) -= rho * Y(H1) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(N) += rho * Y(H1) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(H1) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(H1) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(H1) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(H1) * (table_values.enu + table_values.gamma);
    }


    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void
    fill_rates([[maybe_unused]] const amrex::Real temp,
               [[maybe_unused]] const amrex::Real rho,
               [[maybe_unused]] const amrex::Real rhoy,
               [[maybe_unused]] const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
               [[maybe_unused]] T& rate_eval)
    {

        // Calculate tabular weak rates

        [[maybe_unused]] const amrex::Real log_temp = std::log10(temp);
        [[maybe_unused]] const amrex::Real log_rhoy = std::log10(rhoy);

        rate_eval.enuc_weak = 0.0_rt;

        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            for (int n = 1; n <= NumSpec; ++n) {
                rate_eval.dweak_ydot_dYe(n) = 0.0_rt;
                rate_eval.denuc_weak_dY(n) = 0.0_rt;
            }
            rate_eval.denuc_weak_dT = 0.0_rt;
            rate_eval.denuc_weak_dYe = 0.0_rt;

        }

        {
            // Na21_to_Ne21_electron_capture_weaktab
            rate_Na21_to_Ne21_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ne21_to_Na21_beta_neg_weaktab
            rate_Ne21_to_Na21_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Na22_to_Ne22_electron_capture_weaktab
            rate_Na22_to_Ne22_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ne22_to_Na22_beta_neg_weaktab
            rate_Ne22_to_Na22_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mg23_to_Na23_electron_capture_weaktab
            rate_Mg23_to_Na23_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Na23_to_Mg23_beta_neg_weaktab
            rate_Na23_to_Mg23_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Al25_to_Mg25_electron_capture_weaktab
            rate_Al25_to_Mg25_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mg25_to_Al25_beta_neg_weaktab
            rate_Mg25_to_Al25_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Al26_to_Mg26_electron_capture_weaktab
            rate_Al26_to_Mg26_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mg26_to_Al26_beta_neg_weaktab
            rate_Mg26_to_Al26_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // P29_to_Si29_electron_capture_weaktab
            rate_P29_to_Si29_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Si29_to_P29_beta_neg_weaktab
            rate_Si29_to_P29_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // P30_to_Si30_electron_capture_weaktab
            rate_P30_to_Si30_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Si30_to_P30_beta_neg_weaktab
            rate_Si30_to_P30_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // P31_to_Si31_electron_capture_weaktab
            rate_P31_to_Si31_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Si31_to_P31_beta_neg_weaktab
            rate_Si31_to_P31_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // P32_to_S32_beta_neg_weaktab
            rate_P32_to_S32_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // P32_to_Si32_electron_capture_weaktab
            rate_P32_to_Si32_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // S32_to_P32_electron_capture_weaktab
            rate_S32_to_P32_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Si32_to_P32_beta_neg_weaktab
            rate_Si32_to_P32_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cl33_to_S33_electron_capture_weaktab
            rate_Cl33_to_S33_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // P33_to_S33_beta_neg_weaktab
            rate_P33_to_S33_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // S33_to_Cl33_beta_neg_weaktab
            rate_S33_to_Cl33_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // S33_to_P33_electron_capture_weaktab
            rate_S33_to_P33_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cl34_to_S34_electron_capture_weaktab
            rate_Cl34_to_S34_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // S34_to_Cl34_beta_neg_weaktab
            rate_S34_to_Cl34_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cl35_to_S35_electron_capture_weaktab
            rate_Cl35_to_S35_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // S35_to_Cl35_beta_neg_weaktab
            rate_S35_to_Cl35_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ar36_to_Cl36_electron_capture_weaktab
            rate_Ar36_to_Cl36_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cl36_to_Ar36_beta_neg_weaktab
            rate_Cl36_to_Ar36_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cl36_to_S36_electron_capture_weaktab
            rate_Cl36_to_S36_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // S36_to_Cl36_beta_neg_weaktab
            rate_S36_to_Cl36_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ar37_to_Cl37_electron_capture_weaktab
            rate_Ar37_to_Cl37_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ar37_to_K37_beta_neg_weaktab
            rate_Ar37_to_K37_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cl37_to_Ar37_beta_neg_weaktab
            rate_Cl37_to_Ar37_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // K37_to_Ar37_electron_capture_weaktab
            rate_K37_to_Ar37_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ar38_to_K38_beta_neg_weaktab
            rate_Ar38_to_K38_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // K38_to_Ar38_electron_capture_weaktab
            rate_K38_to_Ar38_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ar39_to_K39_beta_neg_weaktab
            rate_Ar39_to_K39_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // K39_to_Ar39_electron_capture_weaktab
            rate_K39_to_Ar39_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ar40_to_K40_beta_neg_weaktab
            rate_Ar40_to_K40_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ca40_to_K40_electron_capture_weaktab
            rate_Ca40_to_K40_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // K40_to_Ar40_electron_capture_weaktab
            rate_K40_to_Ar40_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // K40_to_Ca40_beta_neg_weaktab
            rate_K40_to_Ca40_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ca41_to_K41_electron_capture_weaktab
            rate_Ca41_to_K41_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // K41_to_Ca41_beta_neg_weaktab
            rate_K41_to_Ca41_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ca43_to_Sc43_beta_neg_weaktab
            rate_Ca43_to_Sc43_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc43_to_Ca43_electron_capture_weaktab
            rate_Sc43_to_Ca43_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ca44_to_Sc44_beta_neg_weaktab
            rate_Ca44_to_Sc44_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc44_to_Ca44_electron_capture_weaktab
            rate_Sc44_to_Ca44_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc44_to_Ti44_beta_neg_weaktab
            rate_Sc44_to_Ti44_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti44_to_Sc44_electron_capture_weaktab
            rate_Ti44_to_Sc44_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co53_to_Fe53_electron_capture_weaktab
            rate_Co53_to_Fe53_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe53_to_Co53_beta_neg_weaktab
            rate_Fe53_to_Co53_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu57_to_Ni57_electron_capture_weaktab
            rate_Cu57_to_Ni57_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni57_to_Cu57_beta_neg_weaktab
            rate_Ni57_to_Cu57_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ca45_to_Sc45_beta_neg_weaktab
            rate_Ca45_to_Sc45_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc45_to_Ca45_electron_capture_weaktab
            rate_Sc45_to_Ca45_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc45_to_Ti45_beta_neg_weaktab
            rate_Sc45_to_Ti45_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti45_to_Sc45_electron_capture_weaktab
            rate_Ti45_to_Sc45_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ca46_to_Sc46_beta_neg_weaktab
            rate_Ca46_to_Sc46_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc46_to_Ca46_electron_capture_weaktab
            rate_Sc46_to_Ca46_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc46_to_Ti46_beta_neg_weaktab
            rate_Sc46_to_Ti46_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti46_to_Sc46_electron_capture_weaktab
            rate_Ti46_to_Sc46_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti46_to_V46_beta_neg_weaktab
            rate_Ti46_to_V46_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V46_to_Ti46_electron_capture_weaktab
            rate_V46_to_Ti46_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ca47_to_Sc47_beta_neg_weaktab
            rate_Ca47_to_Sc47_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc47_to_Ca47_electron_capture_weaktab
            rate_Sc47_to_Ca47_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc47_to_Ti47_beta_neg_weaktab
            rate_Sc47_to_Ti47_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti47_to_Sc47_electron_capture_weaktab
            rate_Ti47_to_Sc47_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti47_to_V47_beta_neg_weaktab
            rate_Ti47_to_V47_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V47_to_Ti47_electron_capture_weaktab
            rate_V47_to_Ti47_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ca48_to_Sc48_beta_neg_weaktab
            rate_Ca48_to_Sc48_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr48_to_V48_electron_capture_weaktab
            rate_Cr48_to_V48_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc48_to_Ca48_electron_capture_weaktab
            rate_Sc48_to_Ca48_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc48_to_Ti48_beta_neg_weaktab
            rate_Sc48_to_Ti48_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti48_to_Sc48_electron_capture_weaktab
            rate_Ti48_to_Sc48_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti48_to_V48_beta_neg_weaktab
            rate_Ti48_to_V48_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V48_to_Cr48_beta_neg_weaktab
            rate_V48_to_Cr48_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V48_to_Ti48_electron_capture_weaktab
            rate_V48_to_Ti48_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr49_to_V49_electron_capture_weaktab
            rate_Cr49_to_V49_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Sc49_to_Ti49_beta_neg_weaktab
            rate_Sc49_to_Ti49_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti49_to_Sc49_electron_capture_weaktab
            rate_Ti49_to_Sc49_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti49_to_V49_beta_neg_weaktab
            rate_Ti49_to_V49_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V49_to_Cr49_beta_neg_weaktab
            rate_V49_to_Cr49_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V49_to_Ti49_electron_capture_weaktab
            rate_V49_to_Ti49_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr50_to_Mn50_beta_neg_weaktab
            rate_Cr50_to_Mn50_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr50_to_V50_electron_capture_weaktab
            rate_Cr50_to_V50_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn50_to_Cr50_electron_capture_weaktab
            rate_Mn50_to_Cr50_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti50_to_V50_beta_neg_weaktab
            rate_Ti50_to_V50_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V50_to_Cr50_beta_neg_weaktab
            rate_V50_to_Cr50_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V50_to_Ti50_electron_capture_weaktab
            rate_V50_to_Ti50_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr51_to_Mn51_beta_neg_weaktab
            rate_Cr51_to_Mn51_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr51_to_V51_electron_capture_weaktab
            rate_Cr51_to_V51_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn51_to_Cr51_electron_capture_weaktab
            rate_Mn51_to_Cr51_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ti51_to_V51_beta_neg_weaktab
            rate_Ti51_to_V51_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V51_to_Cr51_beta_neg_weaktab
            rate_V51_to_Cr51_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V51_to_Ti51_electron_capture_weaktab
            rate_V51_to_Ti51_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr52_to_Mn52_beta_neg_weaktab
            rate_Cr52_to_Mn52_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr52_to_V52_electron_capture_weaktab
            rate_Cr52_to_V52_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe52_to_Mn52_electron_capture_weaktab
            rate_Fe52_to_Mn52_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn52_to_Cr52_electron_capture_weaktab
            rate_Mn52_to_Cr52_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn52_to_Fe52_beta_neg_weaktab
            rate_Mn52_to_Fe52_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // V52_to_Cr52_beta_neg_weaktab
            rate_V52_to_Cr52_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr53_to_Mn53_beta_neg_weaktab
            rate_Cr53_to_Mn53_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe53_to_Mn53_electron_capture_weaktab
            rate_Fe53_to_Mn53_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn53_to_Cr53_electron_capture_weaktab
            rate_Mn53_to_Cr53_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn53_to_Fe53_beta_neg_weaktab
            rate_Mn53_to_Fe53_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co54_to_Fe54_electron_capture_weaktab
            rate_Co54_to_Fe54_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cr54_to_Mn54_beta_neg_weaktab
            rate_Cr54_to_Mn54_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe54_to_Co54_beta_neg_weaktab
            rate_Fe54_to_Co54_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe54_to_Mn54_electron_capture_weaktab
            rate_Fe54_to_Mn54_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn54_to_Cr54_electron_capture_weaktab
            rate_Mn54_to_Cr54_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn54_to_Fe54_beta_neg_weaktab
            rate_Mn54_to_Fe54_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co55_to_Fe55_electron_capture_weaktab
            rate_Co55_to_Fe55_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe55_to_Co55_beta_neg_weaktab
            rate_Fe55_to_Co55_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe55_to_Mn55_electron_capture_weaktab
            rate_Fe55_to_Mn55_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Mn55_to_Fe55_beta_neg_weaktab
            rate_Mn55_to_Fe55_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co56_to_Fe56_electron_capture_weaktab
            rate_Co56_to_Fe56_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co56_to_Ni56_beta_neg_weaktab
            rate_Co56_to_Ni56_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe56_to_Co56_beta_neg_weaktab
            rate_Fe56_to_Co56_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni56_to_Co56_electron_capture_weaktab
            rate_Ni56_to_Co56_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co57_to_Fe57_electron_capture_weaktab
            rate_Co57_to_Fe57_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co57_to_Ni57_beta_neg_weaktab
            rate_Co57_to_Ni57_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe57_to_Co57_beta_neg_weaktab
            rate_Fe57_to_Co57_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni57_to_Co57_electron_capture_weaktab
            rate_Ni57_to_Co57_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co58_to_Fe58_electron_capture_weaktab
            rate_Co58_to_Fe58_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co58_to_Ni58_beta_neg_weaktab
            rate_Co58_to_Ni58_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu58_to_Ni58_electron_capture_weaktab
            rate_Cu58_to_Ni58_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Fe58_to_Co58_beta_neg_weaktab
            rate_Fe58_to_Co58_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni58_to_Co58_electron_capture_weaktab
            rate_Ni58_to_Co58_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni58_to_Cu58_beta_neg_weaktab
            rate_Ni58_to_Cu58_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Co59_to_Ni59_beta_neg_weaktab
            rate_Co59_to_Ni59_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu59_to_Ni59_electron_capture_weaktab
            rate_Cu59_to_Ni59_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni59_to_Co59_electron_capture_weaktab
            rate_Ni59_to_Co59_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni59_to_Cu59_beta_neg_weaktab
            rate_Ni59_to_Cu59_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu60_to_Ni60_electron_capture_weaktab
            rate_Cu60_to_Ni60_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu60_to_Zn60_beta_neg_weaktab
            rate_Cu60_to_Zn60_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni60_to_Cu60_beta_neg_weaktab
            rate_Ni60_to_Cu60_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn60_to_Cu60_electron_capture_weaktab
            rate_Zn60_to_Cu60_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu61_to_Ni61_electron_capture_weaktab
            rate_Cu61_to_Ni61_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu61_to_Zn61_beta_neg_weaktab
            rate_Cu61_to_Zn61_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni61_to_Cu61_beta_neg_weaktab
            rate_Ni61_to_Cu61_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn61_to_Cu61_electron_capture_weaktab
            rate_Zn61_to_Cu61_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu62_to_Ni62_electron_capture_weaktab
            rate_Cu62_to_Ni62_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu62_to_Zn62_beta_neg_weaktab
            rate_Cu62_to_Zn62_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ga62_to_Zn62_electron_capture_weaktab
            rate_Ga62_to_Zn62_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni62_to_Cu62_beta_neg_weaktab
            rate_Ni62_to_Cu62_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn62_to_Cu62_electron_capture_weaktab
            rate_Zn62_to_Cu62_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn62_to_Ga62_beta_neg_weaktab
            rate_Zn62_to_Ga62_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu63_to_Ni63_electron_capture_weaktab
            rate_Cu63_to_Ni63_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu63_to_Zn63_beta_neg_weaktab
            rate_Cu63_to_Zn63_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ga63_to_Zn63_electron_capture_weaktab
            rate_Ga63_to_Zn63_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni63_to_Cu63_beta_neg_weaktab
            rate_Ni63_to_Cu63_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn63_to_Cu63_electron_capture_weaktab
            rate_Zn63_to_Cu63_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn63_to_Ga63_beta_neg_weaktab
            rate_Zn63_to_Ga63_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu64_to_Ni64_electron_capture_weaktab
            rate_Cu64_to_Ni64_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu64_to_Zn64_beta_neg_weaktab
            rate_Cu64_to_Zn64_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ga64_to_Ge64_beta_neg_weaktab
            rate_Ga64_to_Ge64_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ga64_to_Zn64_electron_capture_weaktab
            rate_Ga64_to_Zn64_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ge64_to_Ga64_electron_capture_weaktab
            rate_Ge64_to_Ga64_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni64_to_Cu64_beta_neg_weaktab
            rate_Ni64_to_Cu64_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn64_to_Cu64_electron_capture_weaktab
            rate_Zn64_to_Cu64_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn64_to_Ga64_beta_neg_weaktab
            rate_Zn64_to_Ga64_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Cu65_to_Zn65_beta_neg_weaktab
            rate_Cu65_to_Zn65_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Zn65_to_Cu65_electron_capture_weaktab
            rate_Zn65_to_Cu65_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // n_to_p_beta_neg_weaktab
            rate_n_to_p_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // p_to_n_electron_capture_weaktab
            rate_p_to_n_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

    }

}
#endif
