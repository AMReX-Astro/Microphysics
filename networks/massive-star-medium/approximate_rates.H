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
        amrex::Real r1_gn = rate_eval.screened_rates(k_Ni57_to_n_Ni56_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Ni56_n_to_Ni57_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Ni56_n_to_Ni57_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Ni57_n_to_Ni58_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ni57_to_n_Ni56_derived);
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
        amrex::Real r1_gn = rate_eval.screened_rates(k_Ni57_to_n_Ni56_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Ni58_to_n_Ni57_derived);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ni57_to_n_Ni56_derived);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Ni58_to_n_Ni57_derived);
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
    void rate_Co53_n_n_to_Co55_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Co54_n_to_Co55_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Co54_to_n_Co53_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Co53_n_to_Co54_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Co53_n_to_Co54_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Co54_n_to_Co55_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Co54_to_n_Co53_derived);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Co53_n_n_to_Co55_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Co53_n_n_to_Co55_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co53_n_n_to_Co55_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Co55_to_Co53_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Co54_n_to_Co55_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Co54_to_n_Co53_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Co55_to_n_Co54_derived);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Co54_to_n_Co53_derived);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Co55_to_n_Co54_derived);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Co54_n_to_Co55_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Co55_to_Co53_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Co55_to_Co53_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Co55_to_Co53_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe50_n_n_to_Fe52_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Fe51_n_to_Fe52_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Fe51_to_n_Fe50_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Fe50_n_to_Fe51_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Fe50_n_to_Fe51_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Fe51_n_to_Fe52_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Fe51_to_n_Fe50_derived);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Fe50_n_n_to_Fe52_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Fe50_n_n_to_Fe52_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe50_n_n_to_Fe52_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe52_to_Fe50_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Fe51_n_to_Fe52_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Fe51_to_n_Fe50_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Fe52_to_n_Fe51_derived);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Fe51_to_n_Fe50_derived);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Fe52_to_n_Fe51_derived);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Fe51_n_to_Fe52_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Fe52_to_Fe50_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Fe52_to_Fe50_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Fe52_to_Fe50_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Fe52_n_n_to_Fe54_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Fe53_n_to_Fe54_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Fe53_to_n_Fe52_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Fe52_n_to_Fe53_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Fe52_n_to_Fe53_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Fe53_n_to_Fe54_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Fe53_to_n_Fe52_derived);
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
        amrex::Real r1_gn = rate_eval.screened_rates(k_Fe53_to_n_Fe52_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Fe54_to_n_Fe53_derived);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Fe53_to_n_Fe52_derived);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Fe54_to_n_Fe53_derived);
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
        amrex::Real r1_gn = rate_eval.screened_rates(k_Fe55_to_n_Fe54_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Fe54_n_to_Fe55_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Fe54_n_to_Fe55_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Fe55_n_to_Fe56_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Fe55_to_n_Fe54_derived);
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
        amrex::Real r1_gn = rate_eval.screened_rates(k_Fe55_to_n_Fe54_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Fe56_to_n_Fe55_derived);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Fe55_to_n_Fe54_derived);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Fe56_to_n_Fe55_derived);
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
    void rate_Mn49_n_n_to_Mn51_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Mn50_n_to_Mn51_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Mn50_to_n_Mn49_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Mn49_n_to_Mn50_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Mn49_n_to_Mn50_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Mn50_n_to_Mn51_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Mn50_to_n_Mn49_derived);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Mn49_n_n_to_Mn51_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Mn49_n_n_to_Mn51_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn49_n_n_to_Mn51_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn51_to_Mn49_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Mn50_n_to_Mn51_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Mn50_to_n_Mn49_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Mn51_to_n_Mn50_derived);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Mn50_to_n_Mn49_derived);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Mn51_to_n_Mn50_derived);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Mn50_n_to_Mn51_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Mn51_to_Mn49_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Mn51_to_Mn49_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mn51_to_Mn49_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Mn51_n_n_to_Mn53_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Mn52_n_to_Mn53_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Mn52_to_n_Mn51_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Mn51_n_to_Mn52_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Mn51_n_to_Mn52_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Mn52_n_to_Mn53_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Mn52_to_n_Mn51_derived);
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
        amrex::Real r1_gn = rate_eval.screened_rates(k_Mn52_to_n_Mn51_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Mn53_to_n_Mn52_derived);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Mn52_to_n_Mn51_derived);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Mn53_to_n_Mn52_derived);
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
    void rate_Cr46_n_n_to_Cr48_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Cr47_n_to_Cr48_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Cr47_to_n_Cr46_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Cr46_n_to_Cr47_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Cr46_n_to_Cr47_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Cr47_n_to_Cr48_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Cr47_to_n_Cr46_derived);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Cr46_n_n_to_Cr48_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Cr46_n_n_to_Cr48_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr46_n_n_to_Cr48_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr48_to_Cr46_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Cr47_n_to_Cr48_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Cr47_to_n_Cr46_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Cr48_to_n_Cr47_derived);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Cr47_to_n_Cr46_derived);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Cr48_to_n_Cr47_derived);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Cr47_n_to_Cr48_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Cr48_to_Cr46_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Cr48_to_Cr46_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Cr48_to_Cr46_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Cr48_n_n_to_Cr50_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Cr49_n_to_Cr50_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Cr49_to_n_Cr48_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Cr48_n_to_Cr49_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Cr48_n_to_Cr49_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Cr49_n_to_Cr50_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Cr49_to_n_Cr48_derived);
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
        amrex::Real r1_gn = rate_eval.screened_rates(k_Cr49_to_n_Cr48_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Cr50_to_n_Cr49_derived);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Cr49_to_n_Cr48_derived);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Cr50_to_n_Cr49_derived);
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
        amrex::Real r1_gn = rate_eval.screened_rates(k_Cr51_to_n_Cr50_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Cr50_n_to_Cr51_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Cr50_n_to_Cr51_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Cr51_n_to_Cr52_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Cr51_to_n_Cr50_derived);
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
        amrex::Real r1_gn = rate_eval.screened_rates(k_Cr51_to_n_Cr50_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Cr52_to_n_Cr51_derived);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Cr51_to_n_Cr50_derived);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Cr52_to_n_Cr51_derived);
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
    void rate_V45_n_n_to_V47_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_V46_n_to_V47_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_V46_to_n_V45_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_V45_n_to_V46_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_V45_n_to_V46_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_V46_n_to_V47_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_V46_to_n_V45_derived);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_V45_n_n_to_V47_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_V45_n_n_to_V47_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V45_n_n_to_V47_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V47_to_V45_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_V46_n_to_V47_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_V46_to_n_V45_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_V47_to_n_V46_derived);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_V46_to_n_V45_derived);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_V47_to_n_V46_derived);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_V46_n_to_V47_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_V47_to_V45_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_V47_to_V45_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_V47_to_V45_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_V47_n_n_to_V49_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_V48_n_to_V49_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_V48_to_n_V47_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_V47_n_to_V48_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_V47_n_to_V48_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_V48_n_to_V49_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_V48_to_n_V47_derived);
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
        amrex::Real r1_gn = rate_eval.screened_rates(k_V48_to_n_V47_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_V49_to_n_V48_derived);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_V48_to_n_V47_derived);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_V49_to_n_V48_derived);
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
    void rate_Ti42_n_n_to_Ti44_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Ti43_n_to_Ti44_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Ti43_to_n_Ti42_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Ti42_n_to_Ti43_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Ti42_n_to_Ti43_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Ti43_n_to_Ti44_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ti43_to_n_Ti42_derived);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Ti42_n_n_to_Ti44_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Ti42_n_n_to_Ti44_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti42_n_n_to_Ti44_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti44_to_Ti42_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Ti43_n_to_Ti44_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Ti43_to_n_Ti42_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Ti44_to_n_Ti43_derived);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ti43_to_n_Ti42_derived);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Ti44_to_n_Ti43_derived);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Ti43_n_to_Ti44_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Ti44_to_Ti42_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Ti44_to_Ti42_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ti44_to_Ti42_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ti44_n_n_to_Ti46_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Ti45_n_to_Ti46_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Ti45_to_n_Ti44_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Ti44_n_to_Ti45_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Ti44_n_to_Ti45_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Ti45_n_to_Ti46_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ti45_to_n_Ti44_derived);
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
        amrex::Real r1_gn = rate_eval.screened_rates(k_Ti45_to_n_Ti44_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Ti46_to_n_Ti45_derived);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ti45_to_n_Ti44_derived);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Ti46_to_n_Ti45_derived);
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
    void rate_Sc41_n_n_to_Sc43_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Sc42_n_to_Sc43_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Sc42_to_n_Sc41_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Sc41_n_to_Sc42_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Sc41_n_to_Sc42_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Sc42_n_to_Sc43_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Sc42_to_n_Sc41_derived);
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
        amrex::Real r1_gn = rate_eval.screened_rates(k_Sc42_to_n_Sc41_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Sc43_to_n_Sc42_derived);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Sc42_to_n_Sc41_derived);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Sc43_to_n_Sc42_derived);
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
    void rate_Ca38_n_n_to_Ca40_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Ca39_n_to_Ca40_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Ca39_to_n_Ca38_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Ca38_n_to_Ca39_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Ca38_n_to_Ca39_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Ca39_n_to_Ca40_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ca39_to_n_Ca38_derived);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Ca38_n_n_to_Ca40_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Ca38_n_n_to_Ca40_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca38_n_n_to_Ca40_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca40_to_Ca38_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Ca39_n_to_Ca40_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Ca39_to_n_Ca38_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Ca40_to_n_Ca39_derived);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ca39_to_n_Ca38_derived);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Ca40_to_n_Ca39_derived);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Ca39_n_to_Ca40_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Ca40_to_Ca38_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Ca40_to_Ca38_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ca40_to_Ca38_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ca40_n_n_to_Ca42_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Ca41_n_to_Ca42_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Ca41_to_n_Ca40_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Ca40_n_to_Ca41_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Ca40_n_to_Ca41_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Ca41_n_to_Ca42_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ca41_to_n_Ca40_derived);
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
        amrex::Real r1_gn = rate_eval.screened_rates(k_Ca41_to_n_Ca40_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Ca42_to_n_Ca41_derived);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ca41_to_n_Ca40_derived);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Ca42_to_n_Ca41_derived);
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
    void rate_Ar34_n_n_to_Ar36_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Ar35_n_to_Ar36_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Ar35_to_n_Ar34_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Ar34_n_to_Ar35_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Ar34_n_to_Ar35_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Ar35_n_to_Ar36_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ar35_to_n_Ar34_derived);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Ar34_n_n_to_Ar36_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Ar34_n_n_to_Ar36_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar34_n_n_to_Ar36_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar36_to_Ar34_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Ar35_n_to_Ar36_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Ar35_to_n_Ar34_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Ar36_to_n_Ar35_derived);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ar35_to_n_Ar34_derived);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Ar36_to_n_Ar35_derived);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Ar35_n_to_Ar36_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Ar36_to_Ar34_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Ar36_to_Ar34_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar36_to_Ar34_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar36_n_n_to_Ar38_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Ar37_n_to_Ar38_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Ar37_to_n_Ar36_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Ar36_n_to_Ar37_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Ar36_n_to_Ar37_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Ar37_n_to_Ar38_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ar37_to_n_Ar36_derived);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Ar36_n_n_to_Ar38_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Ar36_n_n_to_Ar38_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar36_n_n_to_Ar38_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Ar38_to_Ar36_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Ar37_n_to_Ar38_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Ar37_to_n_Ar36_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Ar38_to_n_Ar37_derived);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ar37_to_n_Ar36_derived);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Ar38_to_n_Ar37_derived);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Ar37_n_to_Ar38_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Ar38_to_Ar36_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Ar38_to_Ar36_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ar38_to_Ar36_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si28_n_n_to_Si30_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Si29_n_to_Si30_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Si29_to_n_Si28_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_Si28_n_to_Si29_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_Si28_n_to_Si29_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Si29_n_to_Si30_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Si29_to_n_Si28_derived);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Si28_n_n_to_Si30_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Si28_n_n_to_Si30_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si28_n_n_to_Si30_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_Si30_to_Si28_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_Si29_n_to_Si30_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_Si29_to_n_Si28_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_Si30_to_n_Si29_derived);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Si29_to_n_Si28_derived);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Si30_to_n_Si29_derived);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_Si29_n_to_Si30_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_Si30_to_Si28_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_Si30_to_Si28_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Si30_to_Si28_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S30_n_n_to_S32_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_S31_n_to_S32_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_S31_to_n_S30_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_S30_n_to_S31_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_S30_n_to_S31_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_S31_n_to_S32_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_S31_to_n_S30_derived);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_S30_n_n_to_S32_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_S30_n_n_to_S32_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S30_n_n_to_S32_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S32_to_S30_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_S31_n_to_S32_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_S31_to_n_S30_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_S32_to_n_S31_derived);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_S31_to_n_S30_derived);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_S32_to_n_S31_derived);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_S31_n_to_S32_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_S32_to_S30_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_S32_to_S30_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S32_to_S30_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S32_n_n_to_S34_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_S33_n_to_S34_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_S33_to_n_S32_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r1_ng = rate_eval.screened_rates(k_S32_n_to_S33_reaclib);
        rate = 2.0_rt * r1_ng * r2_ng * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_S32_n_to_S33_reaclib);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_S33_n_to_S34_reaclib);
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_S33_to_n_S32_derived);
            drate_dT = 2.0_rt * (dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn));
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_S32_n_n_to_S34_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_S32_n_n_to_S34_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S32_n_n_to_S34_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_S34_to_S32_n_n_approx(T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real drate_dYN{};
        amrex::Real Yn = Y(N);
        amrex::Real r2_ng = rate_eval.screened_rates(k_S33_n_to_S34_reaclib);
        amrex::Real r1_gn = rate_eval.screened_rates(k_S33_to_n_S32_derived);
        amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
        amrex::Real r2_gn = rate_eval.screened_rates(k_S34_to_n_S33_derived);
        rate = r1_gn * r2_gn * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_S33_to_n_S32_derived);
            amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_S34_to_n_S33_derived);
            amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_S33_n_to_S34_reaclib);
            drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
            drate_dYN = -rate * dd * rho * r2_ng;
        }
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.drate_S34_to_S32_n_n_approx_dYN = drate_dYN;
        }
        rate_eval.screened_rates(k_S34_to_S32_n_n_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_S34_to_S32_n_n_approx) = drate_dT;
        }
    }

    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void rate_C12_C12_to_Mg24_approx(T& rate_eval) {

        amrex::Real rate{}, drate_dT{};
        amrex::Real r_pY = rate_eval.screened_rates(k_Na23_p_to_C12_C12_derived);
        amrex::Real r_pa = rate_eval.screened_rates(k_Na23_p_to_He4_Ne20_reaclib);
        amrex::Real r_pg = rate_eval.screened_rates(k_Na23_p_to_Mg24_reaclib);
        amrex::Real dd = 1.0_rt / (r_pY + r_pa + r_pg);
        amrex::Real r_Yg = 0.0_rt;
        amrex::Real r_Yp = rate_eval.screened_rates(k_C12_C12_to_p_Na23_reaclib);
        rate = r_Yg + r_Yp * r_pg * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real drdT_pY = rate_eval.dscreened_rates_dT(k_Na23_p_to_C12_C12_derived);
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
        amrex::Real r_pY = rate_eval.screened_rates(k_Na23_p_to_C12_C12_derived);
        amrex::Real r_pa = rate_eval.screened_rates(k_Na23_p_to_He4_Ne20_reaclib);
        amrex::Real r_pg = rate_eval.screened_rates(k_Na23_p_to_Mg24_reaclib);
        amrex::Real dd = 1.0_rt / (r_pY + r_pa + r_pg);
        amrex::Real r_gY = 0.0_rt;
        amrex::Real r_gp = rate_eval.screened_rates(k_Mg24_to_p_Na23_derived);
        rate = r_gY + r_pY * r_gp * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real drdT_pY = rate_eval.dscreened_rates_dT(k_Na23_p_to_C12_C12_derived);
            amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_Na23_p_to_He4_Ne20_reaclib);
            amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_Na23_p_to_Mg24_reaclib);
            amrex::Real drdT_gY = 0.0_rt;
            amrex::Real drdT_gp = rate_eval.dscreened_rates_dT(k_Mg24_to_p_Na23_derived);
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
        amrex::Real r_pY = rate_eval.screened_rates(k_Na23_p_to_C12_C12_derived);
        amrex::Real r_pa = rate_eval.screened_rates(k_Na23_p_to_He4_Ne20_reaclib);
        amrex::Real r_pg = rate_eval.screened_rates(k_Na23_p_to_Mg24_reaclib);
        amrex::Real dd = 1.0_rt / (r_pY + r_pa + r_pg);
        amrex::Real r_Ya = rate_eval.screened_rates(k_C12_C12_to_He4_Ne20_reaclib);
        amrex::Real r_Yp = rate_eval.screened_rates(k_C12_C12_to_p_Na23_reaclib);
        rate = r_Ya + r_Yp * r_pa * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real drdT_pY = rate_eval.dscreened_rates_dT(k_Na23_p_to_C12_C12_derived);
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
        amrex::Real r_pY = rate_eval.screened_rates(k_Na23_p_to_C12_C12_derived);
        amrex::Real r_pa = rate_eval.screened_rates(k_Na23_p_to_He4_Ne20_reaclib);
        amrex::Real r_pg = rate_eval.screened_rates(k_Na23_p_to_Mg24_reaclib);
        amrex::Real dd = 1.0_rt / (r_pY + r_pa + r_pg);
        amrex::Real r_aY = rate_eval.screened_rates(k_Ne20_He4_to_C12_C12_derived);
        amrex::Real r_ap = rate_eval.screened_rates(k_Ne20_He4_to_p_Na23_derived);
        rate = r_aY + r_pY * r_ap * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real drdT_pY = rate_eval.dscreened_rates_dT(k_Na23_p_to_C12_C12_derived);
            amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_Na23_p_to_He4_Ne20_reaclib);
            amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_Na23_p_to_Mg24_reaclib);
            amrex::Real drdT_aY = rate_eval.dscreened_rates_dT(k_Ne20_He4_to_C12_C12_derived);
            amrex::Real drdT_ap = rate_eval.dscreened_rates_dT(k_Ne20_He4_to_p_Na23_derived);
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
        amrex::Real r_pY = rate_eval.screened_rates(k_Na23_p_to_C12_C12_derived);
        amrex::Real dd = 1.0_rt / (r_pg + r_pa + r_pY);
        amrex::Real r_ag = rate_eval.screened_rates(k_Ne20_He4_to_Mg24_reaclib);
        amrex::Real r_ap = rate_eval.screened_rates(k_Ne20_He4_to_p_Na23_derived);
        rate = r_ag + r_ap * r_pg * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real drdT_ag = rate_eval.dscreened_rates_dT(k_Ne20_He4_to_Mg24_reaclib);
            amrex::Real drdT_ap = rate_eval.dscreened_rates_dT(k_Ne20_He4_to_p_Na23_derived);
            amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_Na23_p_to_Mg24_reaclib);
            amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_Na23_p_to_He4_Ne20_reaclib);
            amrex::Real drdT_pY = rate_eval.dscreened_rates_dT(k_Na23_p_to_C12_C12_derived);
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
        amrex::Real r_pY = rate_eval.screened_rates(k_Na23_p_to_C12_C12_derived);
        amrex::Real dd = 1.0_rt / (r_pg + r_pa + r_pY);
        amrex::Real r_ga = rate_eval.screened_rates(k_Mg24_to_He4_Ne20_derived);
        amrex::Real r_gp = rate_eval.screened_rates(k_Mg24_to_p_Na23_derived);
        rate = r_ga + r_gp * r_pa * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real drdT_ga = rate_eval.dscreened_rates_dT(k_Mg24_to_He4_Ne20_derived);
            amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_Na23_p_to_He4_Ne20_reaclib);
            amrex::Real drdT_gp = rate_eval.dscreened_rates_dT(k_Mg24_to_p_Na23_derived);
            amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_Na23_p_to_Mg24_reaclib);
            amrex::Real drdT_pY = rate_eval.dscreened_rates_dT(k_Na23_p_to_C12_C12_derived);
            drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa + drdT_pY);
        }
        rate_eval.screened_rates(k_Mg24_to_Ne20_He4_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Mg24_to_Ne20_He4_approx) = drate_dT;
        }
    }



    template <typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void
    fill_rates([[maybe_unused]] const tf_t& tfactors,
               [[maybe_unused]] const amrex::Real rho,
               [[maybe_unused]] const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
               [[maybe_unused]] T& rate_eval)
    {

        [[maybe_unused]] amrex::Real rate{};
        [[maybe_unused]] amrex::Real drate_dT{};

        {
            // Ni56_n_n_to_Ni58_approx
            rate_Ni56_n_n_to_Ni58_approx(rate_eval, rho, Y);
        }

        {
            // Ni58_to_Ni56_n_n_approx
            rate_Ni58_to_Ni56_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Co53_n_n_to_Co55_approx
            rate_Co53_n_n_to_Co55_approx(rate_eval, rho, Y);
        }

        {
            // Co55_to_Co53_n_n_approx
            rate_Co55_to_Co53_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Fe50_n_n_to_Fe52_approx
            rate_Fe50_n_n_to_Fe52_approx(rate_eval, rho, Y);
        }

        {
            // Fe52_to_Fe50_n_n_approx
            rate_Fe52_to_Fe50_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Fe52_n_n_to_Fe54_approx
            rate_Fe52_n_n_to_Fe54_approx(rate_eval, rho, Y);
        }

        {
            // Fe54_to_Fe52_n_n_approx
            rate_Fe54_to_Fe52_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Fe54_n_n_to_Fe56_approx
            rate_Fe54_n_n_to_Fe56_approx(rate_eval, rho, Y);
        }

        {
            // Fe56_to_Fe54_n_n_approx
            rate_Fe56_to_Fe54_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Mn49_n_n_to_Mn51_approx
            rate_Mn49_n_n_to_Mn51_approx(rate_eval, rho, Y);
        }

        {
            // Mn51_to_Mn49_n_n_approx
            rate_Mn51_to_Mn49_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Mn51_n_n_to_Mn53_approx
            rate_Mn51_n_n_to_Mn53_approx(rate_eval, rho, Y);
        }

        {
            // Mn53_to_Mn51_n_n_approx
            rate_Mn53_to_Mn51_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Cr46_n_n_to_Cr48_approx
            rate_Cr46_n_n_to_Cr48_approx(rate_eval, rho, Y);
        }

        {
            // Cr48_to_Cr46_n_n_approx
            rate_Cr48_to_Cr46_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Cr48_n_n_to_Cr50_approx
            rate_Cr48_n_n_to_Cr50_approx(rate_eval, rho, Y);
        }

        {
            // Cr50_to_Cr48_n_n_approx
            rate_Cr50_to_Cr48_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Cr50_n_n_to_Cr52_approx
            rate_Cr50_n_n_to_Cr52_approx(rate_eval, rho, Y);
        }

        {
            // Cr52_to_Cr50_n_n_approx
            rate_Cr52_to_Cr50_n_n_approx(rate_eval, rho, Y);
        }

        {
            // V45_n_n_to_V47_approx
            rate_V45_n_n_to_V47_approx(rate_eval, rho, Y);
        }

        {
            // V47_to_V45_n_n_approx
            rate_V47_to_V45_n_n_approx(rate_eval, rho, Y);
        }

        {
            // V47_n_n_to_V49_approx
            rate_V47_n_n_to_V49_approx(rate_eval, rho, Y);
        }

        {
            // V49_to_V47_n_n_approx
            rate_V49_to_V47_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Ti42_n_n_to_Ti44_approx
            rate_Ti42_n_n_to_Ti44_approx(rate_eval, rho, Y);
        }

        {
            // Ti44_to_Ti42_n_n_approx
            rate_Ti44_to_Ti42_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Ti44_n_n_to_Ti46_approx
            rate_Ti44_n_n_to_Ti46_approx(rate_eval, rho, Y);
        }

        {
            // Ti46_to_Ti44_n_n_approx
            rate_Ti46_to_Ti44_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Sc41_n_n_to_Sc43_approx
            rate_Sc41_n_n_to_Sc43_approx(rate_eval, rho, Y);
        }

        {
            // Sc43_to_Sc41_n_n_approx
            rate_Sc43_to_Sc41_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Ca38_n_n_to_Ca40_approx
            rate_Ca38_n_n_to_Ca40_approx(rate_eval, rho, Y);
        }

        {
            // Ca40_to_Ca38_n_n_approx
            rate_Ca40_to_Ca38_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Ca40_n_n_to_Ca42_approx
            rate_Ca40_n_n_to_Ca42_approx(rate_eval, rho, Y);
        }

        {
            // Ca42_to_Ca40_n_n_approx
            rate_Ca42_to_Ca40_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Ar34_n_n_to_Ar36_approx
            rate_Ar34_n_n_to_Ar36_approx(rate_eval, rho, Y);
        }

        {
            // Ar36_to_Ar34_n_n_approx
            rate_Ar36_to_Ar34_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Ar36_n_n_to_Ar38_approx
            rate_Ar36_n_n_to_Ar38_approx(rate_eval, rho, Y);
        }

        {
            // Ar38_to_Ar36_n_n_approx
            rate_Ar38_to_Ar36_n_n_approx(rate_eval, rho, Y);
        }

        {
            // Si28_n_n_to_Si30_approx
            rate_Si28_n_n_to_Si30_approx(rate_eval, rho, Y);
        }

        {
            // Si30_to_Si28_n_n_approx
            rate_Si30_to_Si28_n_n_approx(rate_eval, rho, Y);
        }

        {
            // S30_n_n_to_S32_approx
            rate_S30_n_n_to_S32_approx(rate_eval, rho, Y);
        }

        {
            // S32_to_S30_n_n_approx
            rate_S32_to_S30_n_n_approx(rate_eval, rho, Y);
        }

        {
            // S32_n_n_to_S34_approx
            rate_S32_n_n_to_S34_approx(rate_eval, rho, Y);
        }

        {
            // S34_to_S32_n_n_approx
            rate_S34_to_S32_n_n_approx(rate_eval, rho, Y);
        }

        {
            // C12_C12_to_Mg24_approx
            rate_C12_C12_to_Mg24_approx(rate_eval);
        }

        {
            // Mg24_to_C12_C12_approx
            rate_Mg24_to_C12_C12_approx(rate_eval);
        }

        {
            // C12_C12_to_Ne20_He4_approx
            rate_C12_C12_to_Ne20_He4_approx(rate_eval);
        }

        {
            // Ne20_He4_to_C12_C12_approx
            rate_Ne20_He4_to_C12_C12_approx(rate_eval);
        }

        {
            // Ne20_He4_to_Mg24_approx
            rate_Ne20_He4_to_Mg24_approx(rate_eval);
        }

        {
            // Mg24_to_Ne20_He4_approx
            rate_Mg24_to_Ne20_He4_approx(rate_eval);
        }


    }
}
#endif
