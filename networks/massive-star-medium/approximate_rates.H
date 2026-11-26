#ifndef APPROXIMATE_RATES_H
#define APPROXIMATE_RATES_H

#include <AMReX.H>
#include <AMReX_Print.H>

#include <tfactors.H>
#include <actual_network.H>
#include <rate_type.H>

namespace approximate_rates {

    using namespace Rates;
    using namespace Species;

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni56_n_n_to_Ni58_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Ni57_n_to_Ni58_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Ni57_to_n_Ni56_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Ni56_n_to_Ni57_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Ni56_n_to_Ni57_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Ni57_n_to_Ni58_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ni57_to_n_Ni56_reaclib);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Ni56_n_n_to_Ni58_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Ni56_n_n_to_Ni58_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni56_n_n_to_Ni58_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni58_to_Ni56_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Ni57_n_to_Ni58_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Ni57_to_n_Ni56_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Ni58_to_n_Ni57_reaclib);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ni57_to_n_Ni56_reaclib);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Ni58_to_n_Ni57_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Ni57_n_to_Ni58_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Ni58_to_Ni56_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Ni58_to_Ni56_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni58_to_Ni56_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni58_n_n_to_Ni60_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Ni59_n_to_Ni60_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Ni59_to_n_Ni58_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Ni58_n_to_Ni59_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Ni58_n_to_Ni59_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Ni59_n_to_Ni60_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ni59_to_n_Ni58_reaclib);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Ni58_n_n_to_Ni60_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Ni58_n_n_to_Ni60_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni58_n_n_to_Ni60_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ni60_to_Ni58_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Ni59_n_to_Ni60_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Ni59_to_n_Ni58_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Ni60_to_n_Ni59_reaclib);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ni59_to_n_Ni58_reaclib);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Ni60_to_n_Ni59_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Ni59_n_to_Ni60_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Ni60_to_Ni58_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Ni60_to_Ni58_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ni60_to_Ni58_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co57_n_n_to_Co59_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Co58_n_to_Co59_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Co58_to_n_Co57_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Co57_n_to_Co58_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Co57_n_to_Co58_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Co58_n_to_Co59_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Co58_to_n_Co57_reaclib);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Co57_n_n_to_Co59_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Co57_n_n_to_Co59_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co57_n_n_to_Co59_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co59_to_Co57_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Co58_n_to_Co59_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Co58_to_n_Co57_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Co59_to_n_Co58_reaclib);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Co58_to_n_Co57_reaclib);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Co59_to_n_Co58_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Co58_n_to_Co59_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Co59_to_Co57_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Co59_to_Co57_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co59_to_Co57_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe52_n_n_to_Fe54_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Fe53_n_to_Fe54_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Fe53_to_n_Fe52_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Fe52_n_to_Fe53_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Fe52_n_to_Fe53_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Fe53_n_to_Fe54_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Fe53_to_n_Fe52_reaclib);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Fe52_n_n_to_Fe54_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Fe52_n_n_to_Fe54_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe52_n_n_to_Fe54_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe54_to_Fe52_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Fe53_n_to_Fe54_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Fe53_to_n_Fe52_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Fe54_to_n_Fe53_reaclib);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Fe53_to_n_Fe52_reaclib);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Fe54_to_n_Fe53_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Fe53_n_to_Fe54_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Fe54_to_Fe52_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Fe54_to_Fe52_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe54_to_Fe52_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe54_n_n_to_Fe56_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Fe55_n_to_Fe56_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Fe55_to_n_Fe54_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Fe54_n_to_Fe55_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Fe54_n_to_Fe55_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Fe55_n_to_Fe56_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Fe55_to_n_Fe54_reaclib);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Fe54_n_n_to_Fe56_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Fe54_n_n_to_Fe56_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe54_n_n_to_Fe56_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe56_to_Fe54_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Fe55_n_to_Fe56_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Fe55_to_n_Fe54_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Fe56_to_n_Fe55_reaclib);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Fe55_to_n_Fe54_reaclib);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Fe56_to_n_Fe55_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Fe55_n_to_Fe56_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Fe56_to_Fe54_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Fe56_to_Fe54_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe56_to_Fe54_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe56_n_n_to_Fe58_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Fe57_n_to_Fe58_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Fe57_to_n_Fe56_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Fe56_n_to_Fe57_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Fe56_n_to_Fe57_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Fe57_n_to_Fe58_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Fe57_to_n_Fe56_reaclib);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Fe56_n_n_to_Fe58_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Fe56_n_n_to_Fe58_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe56_n_n_to_Fe58_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe58_to_Fe56_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Fe57_n_to_Fe58_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Fe57_to_n_Fe56_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Fe58_to_n_Fe57_reaclib);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Fe57_to_n_Fe56_reaclib);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Fe58_to_n_Fe57_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Fe57_n_to_Fe58_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Fe58_to_Fe56_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Fe58_to_Fe56_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe58_to_Fe56_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn51_n_n_to_Mn53_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Mn52_n_to_Mn53_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Mn52_to_n_Mn51_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Mn51_n_to_Mn52_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Mn51_n_to_Mn52_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Mn52_n_to_Mn53_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Mn52_to_n_Mn51_reaclib);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Mn51_n_n_to_Mn53_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Mn51_n_n_to_Mn53_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn51_n_n_to_Mn53_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn53_to_Mn51_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Mn52_n_to_Mn53_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Mn52_to_n_Mn51_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Mn53_to_n_Mn52_reaclib);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Mn52_to_n_Mn51_reaclib);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Mn53_to_n_Mn52_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Mn52_n_to_Mn53_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Mn53_to_Mn51_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Mn53_to_Mn51_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn53_to_Mn51_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn53_n_n_to_Mn55_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Mn54_n_to_Mn55_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Mn54_to_n_Mn53_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Mn53_n_to_Mn54_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Mn53_n_to_Mn54_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Mn54_n_to_Mn55_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Mn54_to_n_Mn53_reaclib);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Mn53_n_n_to_Mn55_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Mn53_n_n_to_Mn55_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn53_n_n_to_Mn55_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn55_to_Mn53_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Mn54_n_to_Mn55_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Mn54_to_n_Mn53_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Mn55_to_n_Mn54_reaclib);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Mn54_to_n_Mn53_reaclib);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Mn55_to_n_Mn54_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Mn54_n_to_Mn55_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Mn55_to_Mn53_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Mn55_to_Mn53_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn55_to_Mn53_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr48_n_n_to_Cr50_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Cr49_n_to_Cr50_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Cr49_to_n_Cr48_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Cr48_n_to_Cr49_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Cr48_n_to_Cr49_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Cr49_n_to_Cr50_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Cr49_to_n_Cr48_reaclib);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Cr48_n_n_to_Cr50_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Cr48_n_n_to_Cr50_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr48_n_n_to_Cr50_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr50_to_Cr48_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Cr49_n_to_Cr50_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Cr49_to_n_Cr48_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Cr50_to_n_Cr49_reaclib);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Cr49_to_n_Cr48_reaclib);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Cr50_to_n_Cr49_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Cr49_n_to_Cr50_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Cr50_to_Cr48_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Cr50_to_Cr48_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr50_to_Cr48_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr50_n_n_to_Cr52_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Cr51_n_to_Cr52_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Cr51_to_n_Cr50_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Cr50_n_to_Cr51_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Cr50_n_to_Cr51_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Cr51_n_to_Cr52_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Cr51_to_n_Cr50_reaclib);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Cr50_n_n_to_Cr52_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Cr50_n_n_to_Cr52_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr50_n_n_to_Cr52_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr52_to_Cr50_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Cr51_n_to_Cr52_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Cr51_to_n_Cr50_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Cr52_to_n_Cr51_reaclib);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Cr51_to_n_Cr50_reaclib);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Cr52_to_n_Cr51_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Cr51_n_to_Cr52_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Cr52_to_Cr50_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Cr52_to_Cr50_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr52_to_Cr50_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr52_n_n_to_Cr54_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Cr53_n_to_Cr54_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Cr53_to_n_Cr52_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Cr52_n_to_Cr53_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Cr52_n_to_Cr53_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Cr53_n_to_Cr54_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Cr53_to_n_Cr52_reaclib);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Cr52_n_n_to_Cr54_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Cr52_n_n_to_Cr54_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr52_n_n_to_Cr54_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr54_to_Cr52_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Cr53_n_to_Cr54_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Cr53_to_n_Cr52_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Cr54_to_n_Cr53_reaclib);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Cr53_to_n_Cr52_reaclib);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Cr54_to_n_Cr53_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Cr53_n_to_Cr54_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Cr54_to_Cr52_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Cr54_to_Cr52_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr54_to_Cr52_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V47_n_n_to_V49_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_V48_n_to_V49_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_V48_to_n_V47_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_V47_n_to_V48_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_V47_n_to_V48_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_V48_n_to_V49_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_V48_to_n_V47_reaclib);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_V47_n_n_to_V49_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_V47_n_n_to_V49_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V47_n_n_to_V49_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V49_to_V47_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_V48_n_to_V49_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_V48_to_n_V47_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_V49_to_n_V48_reaclib);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_V48_to_n_V47_reaclib);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_V49_to_n_V48_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_V48_n_to_V49_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_V49_to_V47_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_V49_to_V47_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V49_to_V47_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V49_n_n_to_V51_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_V50_n_to_V51_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_V50_to_n_V49_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_V49_n_to_V50_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_V49_n_to_V50_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_V50_n_to_V51_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_V50_to_n_V49_reaclib);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_V49_n_n_to_V51_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_V49_n_n_to_V51_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V49_n_n_to_V51_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V51_to_V49_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_V50_n_to_V51_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_V50_to_n_V49_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_V51_to_n_V50_reaclib);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_V50_to_n_V49_reaclib);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_V51_to_n_V50_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_V50_n_to_V51_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_V51_to_V49_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_V51_to_V49_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V51_to_V49_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti44_n_n_to_Ti46_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Ti45_n_to_Ti46_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Ti45_to_n_Ti44_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Ti44_n_to_Ti45_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Ti44_n_to_Ti45_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Ti45_n_to_Ti46_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ti45_to_n_Ti44_reaclib);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Ti44_n_n_to_Ti46_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Ti44_n_n_to_Ti46_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti44_n_n_to_Ti46_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti46_to_Ti44_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Ti45_n_to_Ti46_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Ti45_to_n_Ti44_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Ti46_to_n_Ti45_reaclib);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ti45_to_n_Ti44_reaclib);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Ti46_to_n_Ti45_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Ti45_n_to_Ti46_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Ti46_to_Ti44_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Ti46_to_Ti44_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti46_to_Ti44_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti46_n_n_to_Ti48_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Ti47_n_to_Ti48_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Ti47_to_n_Ti46_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Ti46_n_to_Ti47_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Ti46_n_to_Ti47_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Ti47_n_to_Ti48_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ti47_to_n_Ti46_reaclib);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Ti46_n_n_to_Ti48_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Ti46_n_n_to_Ti48_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti46_n_n_to_Ti48_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti48_to_Ti46_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Ti47_n_to_Ti48_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Ti47_to_n_Ti46_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Ti48_to_n_Ti47_reaclib);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ti47_to_n_Ti46_reaclib);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Ti48_to_n_Ti47_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Ti47_n_to_Ti48_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Ti48_to_Ti46_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Ti48_to_Ti46_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti48_to_Ti46_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc41_n_n_to_Sc43_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Sc42_n_to_Sc43_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Sc42_to_n_Sc41_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Sc41_n_to_Sc42_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Sc41_n_to_Sc42_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Sc42_n_to_Sc43_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Sc42_to_n_Sc41_reaclib);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Sc41_n_n_to_Sc43_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Sc41_n_n_to_Sc43_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc41_n_n_to_Sc43_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc43_to_Sc41_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Sc42_n_to_Sc43_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Sc42_to_n_Sc41_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Sc43_to_n_Sc42_reaclib);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Sc42_to_n_Sc41_reaclib);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Sc43_to_n_Sc42_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Sc42_n_to_Sc43_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Sc43_to_Sc41_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Sc43_to_Sc41_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc43_to_Sc41_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc43_n_n_to_Sc45_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Sc44_n_to_Sc45_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Sc44_to_n_Sc43_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Sc43_n_to_Sc44_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Sc43_n_to_Sc44_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Sc44_n_to_Sc45_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Sc44_to_n_Sc43_reaclib);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Sc43_n_n_to_Sc45_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Sc43_n_n_to_Sc45_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc43_n_n_to_Sc45_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Sc45_to_Sc43_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Sc44_n_to_Sc45_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Sc44_to_n_Sc43_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Sc45_to_n_Sc44_reaclib);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Sc44_to_n_Sc43_reaclib);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Sc45_to_n_Sc44_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Sc44_n_to_Sc45_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Sc45_to_Sc43_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Sc45_to_Sc43_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Sc45_to_Sc43_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca40_n_n_to_Ca42_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Ca41_n_to_Ca42_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Ca41_to_n_Ca40_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Ca40_n_to_Ca41_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Ca40_n_to_Ca41_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Ca41_n_to_Ca42_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ca41_to_n_Ca40_reaclib);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Ca40_n_n_to_Ca42_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Ca40_n_n_to_Ca42_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca40_n_n_to_Ca42_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca42_to_Ca40_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Ca41_n_to_Ca42_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Ca41_to_n_Ca40_reaclib);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Ca42_to_n_Ca41_reaclib);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ca41_to_n_Ca40_reaclib);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Ca42_to_n_Ca41_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Ca41_n_to_Ca42_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Ca42_to_Ca40_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Ca42_to_Ca40_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca42_to_Ca40_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_C12_C12_to_Mg24_approx(T& rate_eval) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real r_pY = rate_eval.screened_rates(k_Na23_p_to_C12_C12_reaclib);
        amrex::Real r_pa = rate_eval.screened_rates(k_Na23_p_to_He4_Ne20_reaclib);
        amrex::Real r_pg = rate_eval.screened_rates(k_Na23_p_to_Mg24_reaclib);
        amrex::Real dd = 1.0_rt / (r_pY + r_pa + r_pg);
        amrex::Real r_Yg = 0.0_rt;
        amrex::Real r_Yp = rate_eval.screened_rates(k_C12_C12_to_p_Na23_reaclib);
        rate = r_Yg + r_Yp * r_pg * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real drdT_pY = rate_eval.dscreened_rates_dT(k_Na23_p_to_C12_C12_reaclib);
            amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_Na23_p_to_He4_Ne20_reaclib);
            amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_Na23_p_to_Mg24_reaclib);
            amrex::Real drdT_Yg = 0.0_rt;
            amrex::Real drdT_Yp = rate_eval.dscreened_rates_dT(k_C12_C12_to_p_Na23_reaclib);
            drate_dT = drdT_Yg + drdT_Yp * r_pg * dd + r_Yp * drdT_pg * dd - r_Yp * r_pg * dd * dd * (drdT_pY + drdT_pa + drdT_pg);
        }
        rate_eval.screened_rates(k_C12_C12_to_Mg24_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_C12_C12_to_Mg24_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mg24_to_C12_C12_approx(T& rate_eval) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real r_pY = rate_eval.screened_rates(k_Na23_p_to_C12_C12_reaclib);
        amrex::Real r_pa = rate_eval.screened_rates(k_Na23_p_to_He4_Ne20_reaclib);
        amrex::Real r_pg = rate_eval.screened_rates(k_Na23_p_to_Mg24_reaclib);
        amrex::Real dd = 1.0_rt / (r_pY + r_pa + r_pg);
        amrex::Real r_gY = 0.0_rt;
        amrex::Real r_gp = rate_eval.screened_rates(k_Mg24_to_p_Na23_reaclib);
        rate = r_gY + r_pY * r_gp * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real drdT_pY = rate_eval.dscreened_rates_dT(k_Na23_p_to_C12_C12_reaclib);
            amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_Na23_p_to_He4_Ne20_reaclib);
            amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_Na23_p_to_Mg24_reaclib);
            amrex::Real drdT_gY = 0.0_rt;
            amrex::Real drdT_gp = rate_eval.dscreened_rates_dT(k_Mg24_to_p_Na23_reaclib);
            drate_dT = drdT_gY + drdT_pY * r_gp * dd + r_pY * drdT_gp * dd - r_pY * r_gp * dd * dd * (drdT_pY + drdT_pa + drdT_pg);
        }
        rate_eval.screened_rates(k_Mg24_to_C12_C12_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg24_to_C12_C12_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_C12_C12_to_Ne20_He4_approx(T& rate_eval) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real r_pY = rate_eval.screened_rates(k_Na23_p_to_C12_C12_reaclib);
        amrex::Real r_pa = rate_eval.screened_rates(k_Na23_p_to_He4_Ne20_reaclib);
        amrex::Real r_pg = rate_eval.screened_rates(k_Na23_p_to_Mg24_reaclib);
        amrex::Real dd = 1.0_rt / (r_pY + r_pa + r_pg);
        amrex::Real r_Ya = rate_eval.screened_rates(k_C12_C12_to_He4_Ne20_reaclib);
        amrex::Real r_Yp = rate_eval.screened_rates(k_C12_C12_to_p_Na23_reaclib);
        rate = r_Ya + r_Yp * r_pa * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real drdT_pY = rate_eval.dscreened_rates_dT(k_Na23_p_to_C12_C12_reaclib);
            amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_Na23_p_to_He4_Ne20_reaclib);
            amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_Na23_p_to_Mg24_reaclib);
            amrex::Real drdT_Ya = rate_eval.dscreened_rates_dT(k_C12_C12_to_He4_Ne20_reaclib);
            amrex::Real drdT_Yp = rate_eval.dscreened_rates_dT(k_C12_C12_to_p_Na23_reaclib);
            drate_dT = drdT_Ya + drdT_Yp * r_pa * dd + r_Ya * drdT_pa * dd - r_Yp * r_pa * dd * dd * (drdT_pY + drdT_pa + drdT_pg);
        }
        rate_eval.screened_rates(k_C12_C12_to_Ne20_He4_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_C12_C12_to_Ne20_He4_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ne20_He4_to_C12_C12_approx(T& rate_eval) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real r_pY = rate_eval.screened_rates(k_Na23_p_to_C12_C12_reaclib);
        amrex::Real r_pa = rate_eval.screened_rates(k_Na23_p_to_He4_Ne20_reaclib);
        amrex::Real r_pg = rate_eval.screened_rates(k_Na23_p_to_Mg24_reaclib);
        amrex::Real dd = 1.0_rt / (r_pY + r_pa + r_pg);
        amrex::Real r_aY = rate_eval.screened_rates(k_Ne20_He4_to_C12_C12_reaclib);
        amrex::Real r_ap = rate_eval.screened_rates(k_Ne20_He4_to_p_Na23_reaclib);
        rate = r_aY + r_pY * r_ap * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real drdT_pY = rate_eval.dscreened_rates_dT(k_Na23_p_to_C12_C12_reaclib);
            amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_Na23_p_to_He4_Ne20_reaclib);
            amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_Na23_p_to_Mg24_reaclib);
            amrex::Real drdT_aY = rate_eval.dscreened_rates_dT(k_Ne20_He4_to_C12_C12_reaclib);
            amrex::Real drdT_ap = rate_eval.dscreened_rates_dT(k_Ne20_He4_to_p_Na23_reaclib);
            drate_dT = drdT_aY + drdT_pY * r_ap * dd + r_pY * drdT_ap * dd - r_pY * r_ap * dd * dd * (drdT_pY + drdT_pa + drdT_pg);
        }
        rate_eval.screened_rates(k_Ne20_He4_to_C12_C12_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne20_He4_to_C12_C12_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ne20_He4_to_Mg24_approx(T& rate_eval) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real r_pg = rate_eval.screened_rates(k_Na23_p_to_Mg24_reaclib);
        amrex::Real r_pa = rate_eval.screened_rates(k_Na23_p_to_He4_Ne20_reaclib);
        amrex::Real r_pY = rate_eval.screened_rates(k_Na23_p_to_C12_C12_reaclib);
        amrex::Real dd = 1.0_rt / (r_pg + r_pa + r_pY);
        amrex::Real r_ag = rate_eval.screened_rates(k_Ne20_He4_to_Mg24_reaclib);
        amrex::Real r_ap = rate_eval.screened_rates(k_Ne20_He4_to_p_Na23_reaclib);
        rate = r_ag + r_ap * r_pg * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real drdT_ag = rate_eval.dscreened_rates_dT(k_Ne20_He4_to_Mg24_reaclib);
            amrex::Real drdT_ap = rate_eval.dscreened_rates_dT(k_Ne20_He4_to_p_Na23_reaclib);
            amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_Na23_p_to_Mg24_reaclib);
            amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_Na23_p_to_He4_Ne20_reaclib);
            amrex::Real drdT_pY = rate_eval.dscreened_rates_dT(k_Na23_p_to_C12_C12_reaclib);
            drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa + drdT_pY);
        }
        rate_eval.screened_rates(k_Ne20_He4_to_Mg24_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne20_He4_to_Mg24_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mg24_to_Ne20_He4_approx(T& rate_eval) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real r_pg = rate_eval.screened_rates(k_Na23_p_to_Mg24_reaclib);
        amrex::Real r_pa = rate_eval.screened_rates(k_Na23_p_to_He4_Ne20_reaclib);
        amrex::Real r_pY = rate_eval.screened_rates(k_Na23_p_to_C12_C12_reaclib);
        amrex::Real dd = 1.0_rt / (r_pg + r_pa + r_pY);
        amrex::Real r_ga = rate_eval.screened_rates(k_Mg24_to_He4_Ne20_reaclib);
        amrex::Real r_gp = rate_eval.screened_rates(k_Mg24_to_p_Na23_reaclib);
        rate = r_ga + r_gp * r_pa * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real drdT_ga = rate_eval.dscreened_rates_dT(k_Mg24_to_He4_Ne20_reaclib);
            amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_Na23_p_to_He4_Ne20_reaclib);
            amrex::Real drdT_gp = rate_eval.dscreened_rates_dT(k_Mg24_to_p_Na23_reaclib);
            amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_Na23_p_to_Mg24_reaclib);
            amrex::Real drdT_pY = rate_eval.dscreened_rates_dT(k_Na23_p_to_C12_C12_reaclib);
            drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa + drdT_pY);
        }
        rate_eval.screened_rates(k_Mg24_to_Ne20_He4_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg24_to_Ne20_He4_approx) = drate_dT;
        }
    }


}
#endif
