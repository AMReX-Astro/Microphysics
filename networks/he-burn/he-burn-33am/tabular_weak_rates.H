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
            // Co57_to_Ni57_beta_neg_weaktab
            rate_Co57_to_Ni57_beta_neg_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
        }

        {
            // Ni57_to_Co57_electron_capture_weaktab
            rate_Ni57_to_Co57_electron_capture_weaktab(rate_eval, temp, log_temp, rho, rhoy, log_rhoy, Y);
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
