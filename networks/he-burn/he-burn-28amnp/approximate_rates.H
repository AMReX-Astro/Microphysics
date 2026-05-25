#ifndef APPROXIMATE_RATES_H
#define APPROXIMATE_RATES_H

#include <AMReX.H>
#include <AMReX_Print.H>

#include <tfactors.H>
#include <actual_network.H>
#include <rate_type.H>

using namespace Rates;
using namespace Species;

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_S32_He4_to_Ar36_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_pg = rate_eval.screened_rates(k_p_Cl35_to_Ar36_reaclib);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_Cl35_to_He4_S32_reaclib);
    amrex::Real r_pY = 0.0_rt;
    amrex::Real dd = 1.0_rt / (r_pg + r_pa + r_pY);
    amrex::Real r_ag = rate_eval.screened_rates(k_He4_S32_to_Ar36_reaclib);
    amrex::Real r_ap = rate_eval.screened_rates(k_He4_S32_to_p_Cl35_derived);
    rate = r_ag + r_ap * r_pg * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ag = rate_eval.dscreened_rates_dT(k_He4_S32_to_Ar36_reaclib);
        amrex::Real drdT_ap = rate_eval.dscreened_rates_dT(k_He4_S32_to_p_Cl35_derived);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_Cl35_to_Ar36_reaclib);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_Cl35_to_He4_S32_reaclib);
        amrex::Real drdT_pY = 0.0_rt;
        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa + drdT_pY);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ar36_to_S32_He4_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_pg = rate_eval.screened_rates(k_p_Cl35_to_Ar36_reaclib);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_Cl35_to_He4_S32_reaclib);
    amrex::Real r_pY = 0.0_rt;
    amrex::Real dd = 1.0_rt / (r_pg + r_pa + r_pY);
    amrex::Real r_ga = rate_eval.screened_rates(k_Ar36_to_He4_S32_derived);
    amrex::Real r_gp = rate_eval.screened_rates(k_Ar36_to_p_Cl35_derived);
    rate = r_ga + r_gp * r_pa * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ga = rate_eval.dscreened_rates_dT(k_Ar36_to_He4_S32_derived);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_Cl35_to_He4_S32_reaclib);
        amrex::Real drdT_gp = rate_eval.dscreened_rates_dT(k_Ar36_to_p_Cl35_derived);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_Cl35_to_Ar36_reaclib);
        amrex::Real drdT_pY = 0.0_rt;
        drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa + drdT_pY);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ar36_He4_to_Ca40_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_pg = rate_eval.screened_rates(k_p_K39_to_Ca40_reaclib);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_K39_to_He4_Ar36_reaclib);
    amrex::Real r_pY = 0.0_rt;
    amrex::Real dd = 1.0_rt / (r_pg + r_pa + r_pY);
    amrex::Real r_ag = rate_eval.screened_rates(k_He4_Ar36_to_Ca40_reaclib);
    amrex::Real r_ap = rate_eval.screened_rates(k_He4_Ar36_to_p_K39_derived);
    rate = r_ag + r_ap * r_pg * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ag = rate_eval.dscreened_rates_dT(k_He4_Ar36_to_Ca40_reaclib);
        amrex::Real drdT_ap = rate_eval.dscreened_rates_dT(k_He4_Ar36_to_p_K39_derived);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_K39_to_Ca40_reaclib);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_K39_to_He4_Ar36_reaclib);
        amrex::Real drdT_pY = 0.0_rt;
        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa + drdT_pY);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ca40_to_Ar36_He4_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_pg = rate_eval.screened_rates(k_p_K39_to_Ca40_reaclib);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_K39_to_He4_Ar36_reaclib);
    amrex::Real r_pY = 0.0_rt;
    amrex::Real dd = 1.0_rt / (r_pg + r_pa + r_pY);
    amrex::Real r_ga = rate_eval.screened_rates(k_Ca40_to_He4_Ar36_derived);
    amrex::Real r_gp = rate_eval.screened_rates(k_Ca40_to_p_K39_derived);
    rate = r_ga + r_gp * r_pa * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ga = rate_eval.dscreened_rates_dT(k_Ca40_to_He4_Ar36_derived);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_K39_to_He4_Ar36_reaclib);
        amrex::Real drdT_gp = rate_eval.dscreened_rates_dT(k_Ca40_to_p_K39_derived);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_K39_to_Ca40_reaclib);
        amrex::Real drdT_pY = 0.0_rt;
        drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa + drdT_pY);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ca40_He4_to_Ti44_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_pg = rate_eval.screened_rates(k_p_Sc43_to_Ti44_reaclib);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_Sc43_to_He4_Ca40_reaclib);
    amrex::Real r_pY = 0.0_rt;
    amrex::Real dd = 1.0_rt / (r_pg + r_pa + r_pY);
    amrex::Real r_ag = rate_eval.screened_rates(k_He4_Ca40_to_Ti44_reaclib);
    amrex::Real r_ap = rate_eval.screened_rates(k_He4_Ca40_to_p_Sc43_derived);
    rate = r_ag + r_ap * r_pg * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ag = rate_eval.dscreened_rates_dT(k_He4_Ca40_to_Ti44_reaclib);
        amrex::Real drdT_ap = rate_eval.dscreened_rates_dT(k_He4_Ca40_to_p_Sc43_derived);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_Sc43_to_Ti44_reaclib);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_Sc43_to_He4_Ca40_reaclib);
        amrex::Real drdT_pY = 0.0_rt;
        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa + drdT_pY);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ti44_to_Ca40_He4_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_pg = rate_eval.screened_rates(k_p_Sc43_to_Ti44_reaclib);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_Sc43_to_He4_Ca40_reaclib);
    amrex::Real r_pY = 0.0_rt;
    amrex::Real dd = 1.0_rt / (r_pg + r_pa + r_pY);
    amrex::Real r_ga = rate_eval.screened_rates(k_Ti44_to_He4_Ca40_derived);
    amrex::Real r_gp = rate_eval.screened_rates(k_Ti44_to_p_Sc43_derived);
    rate = r_ga + r_gp * r_pa * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ga = rate_eval.dscreened_rates_dT(k_Ti44_to_He4_Ca40_derived);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_Sc43_to_He4_Ca40_reaclib);
        amrex::Real drdT_gp = rate_eval.dscreened_rates_dT(k_Ti44_to_p_Sc43_derived);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_Sc43_to_Ti44_reaclib);
        amrex::Real drdT_pY = 0.0_rt;
        drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa + drdT_pY);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ti44_He4_to_Cr48_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_pg = rate_eval.screened_rates(k_p_V47_to_Cr48_reaclib);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_V47_to_He4_Ti44_derived);
    amrex::Real r_pY = 0.0_rt;
    amrex::Real dd = 1.0_rt / (r_pg + r_pa + r_pY);
    amrex::Real r_ag = rate_eval.screened_rates(k_He4_Ti44_to_Cr48_reaclib);
    amrex::Real r_ap = rate_eval.screened_rates(k_He4_Ti44_to_p_V47_reaclib);
    rate = r_ag + r_ap * r_pg * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ag = rate_eval.dscreened_rates_dT(k_He4_Ti44_to_Cr48_reaclib);
        amrex::Real drdT_ap = rate_eval.dscreened_rates_dT(k_He4_Ti44_to_p_V47_reaclib);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_V47_to_Cr48_reaclib);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_V47_to_He4_Ti44_derived);
        amrex::Real drdT_pY = 0.0_rt;
        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa + drdT_pY);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Cr48_to_Ti44_He4_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_pg = rate_eval.screened_rates(k_p_V47_to_Cr48_reaclib);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_V47_to_He4_Ti44_derived);
    amrex::Real r_pY = 0.0_rt;
    amrex::Real dd = 1.0_rt / (r_pg + r_pa + r_pY);
    amrex::Real r_ga = rate_eval.screened_rates(k_Cr48_to_He4_Ti44_derived);
    amrex::Real r_gp = rate_eval.screened_rates(k_Cr48_to_p_V47_derived);
    rate = r_ga + r_gp * r_pa * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ga = rate_eval.dscreened_rates_dT(k_Cr48_to_He4_Ti44_derived);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_V47_to_He4_Ti44_derived);
        amrex::Real drdT_gp = rate_eval.dscreened_rates_dT(k_Cr48_to_p_V47_derived);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_V47_to_Cr48_reaclib);
        amrex::Real drdT_pY = 0.0_rt;
        drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa + drdT_pY);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Fe52_n_n_to_Fe54_approx(const T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real Yn = Y(N);
    amrex::Real r1_ng = rate_eval.screened_rates(k_n_Fe52_to_Fe53_reaclib);
    amrex::Real r2_ng = rate_eval.screened_rates(k_n_Fe53_to_Fe54_reaclib);
    amrex::Real r1_gn = rate_eval.screened_rates(k_Fe53_to_n_Fe52_derived);
    amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
    rate = r1_ng * r2_ng * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_n_Fe52_to_Fe53_reaclib);
        amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_n_Fe53_to_Fe54_reaclib);
        amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Fe53_to_n_Fe52_derived);
        drate_dT = dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Fe54_to_Fe52_n_n_approx(const T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real Yn = Y(N);
    amrex::Real r1_gn = rate_eval.screened_rates(k_Fe53_to_n_Fe52_derived);
    amrex::Real r2_gn = rate_eval.screened_rates(k_Fe54_to_n_Fe53_derived);
    amrex::Real r2_ng = rate_eval.screened_rates(k_n_Fe53_to_Fe54_reaclib);
    amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
    rate = r1_gn * r2_gn * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Fe53_to_n_Fe52_derived);
        amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Fe54_to_n_Fe53_derived);
        amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_n_Fe53_to_Fe54_reaclib);
        drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Fe54_n_n_to_Fe56_approx(const T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real Yn = Y(N);
    amrex::Real r1_ng = rate_eval.screened_rates(k_n_Fe54_to_Fe55_reaclib);
    amrex::Real r2_ng = rate_eval.screened_rates(k_n_Fe55_to_Fe56_reaclib);
    amrex::Real r1_gn = rate_eval.screened_rates(k_Fe55_to_n_Fe54_derived);
    amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
    rate = r1_ng * r2_ng * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_n_Fe54_to_Fe55_reaclib);
        amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_n_Fe55_to_Fe56_reaclib);
        amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Fe55_to_n_Fe54_derived);
        drate_dT = dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Fe56_to_Fe54_n_n_approx(const T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real Yn = Y(N);
    amrex::Real r1_gn = rate_eval.screened_rates(k_Fe55_to_n_Fe54_derived);
    amrex::Real r2_gn = rate_eval.screened_rates(k_Fe56_to_n_Fe55_derived);
    amrex::Real r2_ng = rate_eval.screened_rates(k_n_Fe55_to_Fe56_reaclib);
    amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
    rate = r1_gn * r2_gn * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Fe55_to_n_Fe54_derived);
        amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Fe56_to_n_Fe55_derived);
        amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_n_Fe55_to_Fe56_reaclib);
        drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ni56_n_n_to_Ni58_approx(const T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real Yn = Y(N);
    amrex::Real r1_ng = rate_eval.screened_rates(k_n_Ni56_to_Ni57_reaclib);
    amrex::Real r2_ng = rate_eval.screened_rates(k_n_Ni57_to_Ni58_reaclib);
    amrex::Real r1_gn = rate_eval.screened_rates(k_Ni57_to_n_Ni56_derived);
    amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
    rate = r1_ng * r2_ng * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dr1dT_ng = rate_eval.dscreened_rates_dT(k_n_Ni56_to_Ni57_reaclib);
        amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_n_Ni57_to_Ni58_reaclib);
        amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ni57_to_n_Ni56_derived);
        drate_dT = dr1dT_ng * r2_ng * dd + r1_ng * dr2dT_ng * dd - r1_ng * r2_ng * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ni58_to_Ni56_n_n_approx(const T& rate_eval, const amrex::Real rho, const amrex::Array1D<amrex::Real, 1, NumSpec>& Y, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real Yn = Y(N);
    amrex::Real r1_gn = rate_eval.screened_rates(k_Ni57_to_n_Ni56_derived);
    amrex::Real r2_gn = rate_eval.screened_rates(k_Ni58_to_n_Ni57_derived);
    amrex::Real r2_ng = rate_eval.screened_rates(k_n_Ni57_to_Ni58_reaclib);
    amrex::Real dd = 1.0_rt / (rho * Yn * r2_ng + r1_gn);
    rate = r1_gn * r2_gn * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real dr1dT_gn = rate_eval.dscreened_rates_dT(k_Ni57_to_n_Ni56_derived);
        amrex::Real dr2dT_gn = rate_eval.dscreened_rates_dT(k_Ni58_to_n_Ni57_derived);
        amrex::Real dr2dT_ng = rate_eval.dscreened_rates_dT(k_n_Ni57_to_Ni58_reaclib);
        drate_dT = dr1dT_gn * r2_gn * dd + r1_gn * dr2dT_gn * dd - r1_gn * r2_gn * dd * dd * (rho * Yn * dr2dT_ng + dr1dT_gn);
    }
}



template <int do_T_derivatives, typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void
fill_approx_rates([[maybe_unused]] const tf_t& tfactors,
                  [[maybe_unused]] const amrex::Real rho,
                  [[maybe_unused]] const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
                  [[maybe_unused]] T& rate_eval)
{

    [[maybe_unused]] amrex::Real rate{};
    [[maybe_unused]] amrex::Real drate_dT{};

    rate_S32_He4_to_Ar36_approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_S32_He4_to_Ar36_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_S32_He4_to_Ar36_approx) = drate_dT;
    }

    rate_Ar36_to_S32_He4_approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_Ar36_to_S32_He4_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ar36_to_S32_He4_approx) = drate_dT;
    }

    rate_Ar36_He4_to_Ca40_approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_Ar36_He4_to_Ca40_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ar36_He4_to_Ca40_approx) = drate_dT;
    }

    rate_Ca40_to_Ar36_He4_approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_Ca40_to_Ar36_He4_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ca40_to_Ar36_He4_approx) = drate_dT;
    }

    rate_Ca40_He4_to_Ti44_approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_Ca40_He4_to_Ti44_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ca40_He4_to_Ti44_approx) = drate_dT;
    }

    rate_Ti44_to_Ca40_He4_approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_Ti44_to_Ca40_He4_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ti44_to_Ca40_He4_approx) = drate_dT;
    }

    rate_Ti44_He4_to_Cr48_approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_Ti44_He4_to_Cr48_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ti44_He4_to_Cr48_approx) = drate_dT;
    }

    rate_Cr48_to_Ti44_He4_approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_Cr48_to_Ti44_He4_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Cr48_to_Ti44_He4_approx) = drate_dT;
    }

    rate_Fe52_n_n_to_Fe54_approx<T>(rate_eval, rho, Y, rate, drate_dT);
    rate_eval.screened_rates(k_Fe52_n_n_to_Fe54_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe52_n_n_to_Fe54_approx) = drate_dT;
    }

    rate_Fe54_to_Fe52_n_n_approx<T>(rate_eval, rho, Y, rate, drate_dT);
    rate_eval.screened_rates(k_Fe54_to_Fe52_n_n_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe54_to_Fe52_n_n_approx) = drate_dT;
    }

    rate_Fe54_n_n_to_Fe56_approx<T>(rate_eval, rho, Y, rate, drate_dT);
    rate_eval.screened_rates(k_Fe54_n_n_to_Fe56_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe54_n_n_to_Fe56_approx) = drate_dT;
    }

    rate_Fe56_to_Fe54_n_n_approx<T>(rate_eval, rho, Y, rate, drate_dT);
    rate_eval.screened_rates(k_Fe56_to_Fe54_n_n_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe56_to_Fe54_n_n_approx) = drate_dT;
    }

    rate_Ni56_n_n_to_Ni58_approx<T>(rate_eval, rho, Y, rate, drate_dT);
    rate_eval.screened_rates(k_Ni56_n_n_to_Ni58_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni56_n_n_to_Ni58_approx) = drate_dT;
    }

    rate_Ni58_to_Ni56_n_n_approx<T>(rate_eval, rho, Y, rate, drate_dT);
    rate_eval.screened_rates(k_Ni58_to_Ni56_n_n_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni58_to_Ni56_n_n_approx) = drate_dT;
    }


}

#endif
