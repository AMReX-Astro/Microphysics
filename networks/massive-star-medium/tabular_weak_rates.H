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
    void rate_Al27_to_Si27_beta_neg_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Al27 --> Si27

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Al27_Si27_meta, j_Al27_Si27_rhoy, j_Al27_Si27_temp, j_Al27_Si27_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Al27_to_Si27_beta_neg_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Al27_to_Si27_beta_neg_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Al27 and Si27
            // this is the derivative of ∂Y(Al27)/∂t = -Y(Al27) λ and
            //                           ∂Y(Si27)/∂t = +Y(Al27) λ
            rate_eval.dweak_ydot_dYe(Al27) -= rho * Y(Al27) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Si27) += rho * Y(Al27) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Al27) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Al27) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Al27) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Al27) * (table_values.enu + table_values.gamma);
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
    void rate_Si27_to_Al27_electron_capture_weaktab(T& rate_eval, const amrex::Real temp, const amrex::Real log_temp, const amrex::Real rho, const amrex::Real rhoy, const amrex::Real log_rhoy, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        // Si27 --> Al27

        weak_rate_t table_values{};
        constexpr int do_T_derivatives = std::is_same_v<T, rate_derivs_t>;
        tabular_evaluate<do_T_derivatives>(j_Si27_Al27_meta, j_Si27_Al27_rhoy, j_Si27_Al27_temp, j_Si27_Al27_data,
                                            log_rhoy, rhoy, log_temp, temp, table_values);

        rate_eval.screened_rates(k_Si27_to_Al27_electron_capture_weaktab) = table_values.rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si27_to_Al27_electron_capture_weaktab) = table_values.drate_dT;

            // accumulate ∂/∂Y_e contributions to Si27 and Al27
            // this is the derivative of ∂Y(Si27)/∂t = -Y(Si27) λ and
            //                           ∂Y(Al27)/∂t = +Y(Si27) λ
            rate_eval.dweak_ydot_dYe(Si27) -= rho * Y(Si27) * table_values.drate_drhoye;
            rate_eval.dweak_ydot_dYe(Al27) += rho * Y(Si27) * table_values.drate_drhoye;

            // also accumulate the derivatives of ε_{ν,weak} with respect to T and Y
            rate_eval.denuc_weak_dY(Si27) += C::n_A * (table_values.enu + table_values.gamma);
            rate_eval.denuc_weak_dT += C::n_A * Y(Si27) * table_values.denu_dT;
            // finally the derivatives of ε_{ν,weak} with respect to Ye
            rate_eval.denuc_weak_dYe += C::n_A * rho * Y(Si27) * table_values.denu_drhoye;
        }

        rate_eval.enuc_weak += C::n_A * Y(Si27) * (table_values.enu + table_values.gamma);
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
            // Al27_to_Si27_beta_neg_weaktab
            rate_Al27_to_Si27_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
            // P31_to_Si31_electron_capture_weaktab
            rate_P31_to_Si31_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Si31_to_P31_beta_neg_weaktab
            rate_Si31_to_P31_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
            // Cl35_to_S35_electron_capture_weaktab
            rate_Cl35_to_S35_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // S35_to_Cl35_beta_neg_weaktab
            rate_S35_to_Cl35_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
            // n_to_p_beta_neg_weaktab
            rate_n_to_p_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // p_to_n_electron_capture_weaktab
            rate_p_to_n_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Si27_to_Al27_electron_capture_weaktab
            rate_Si27_to_Al27_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

    }

}
#endif
