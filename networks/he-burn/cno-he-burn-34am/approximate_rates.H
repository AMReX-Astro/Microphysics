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
void rate_Cr48_He4_to_Fe52_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_pg = rate_eval.screened_rates(k_p_Mn51_to_Fe52_reaclib);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_Mn51_to_He4_Cr48_derived);
    amrex::Real r_pY = 0.0_rt;
    amrex::Real dd = 1.0_rt / (r_pg + r_pa + r_pY);
    amrex::Real r_ag = rate_eval.screened_rates(k_He4_Cr48_to_Fe52_reaclib);
    amrex::Real r_ap = rate_eval.screened_rates(k_He4_Cr48_to_p_Mn51_reaclib);
    rate = r_ag + r_ap * r_pg * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ag = rate_eval.dscreened_rates_dT(k_He4_Cr48_to_Fe52_reaclib);
        amrex::Real drdT_ap = rate_eval.dscreened_rates_dT(k_He4_Cr48_to_p_Mn51_reaclib);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_Mn51_to_Fe52_reaclib);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_Mn51_to_He4_Cr48_derived);
        amrex::Real drdT_pY = 0.0_rt;
        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa + drdT_pY);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Fe52_to_Cr48_He4_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_pg = rate_eval.screened_rates(k_p_Mn51_to_Fe52_reaclib);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_Mn51_to_He4_Cr48_derived);
    amrex::Real r_pY = 0.0_rt;
    amrex::Real dd = 1.0_rt / (r_pg + r_pa + r_pY);
    amrex::Real r_ga = rate_eval.screened_rates(k_Fe52_to_He4_Cr48_derived);
    amrex::Real r_gp = rate_eval.screened_rates(k_Fe52_to_p_Mn51_derived);
    rate = r_ga + r_gp * r_pa * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ga = rate_eval.dscreened_rates_dT(k_Fe52_to_He4_Cr48_derived);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_Mn51_to_He4_Cr48_derived);
        amrex::Real drdT_gp = rate_eval.dscreened_rates_dT(k_Fe52_to_p_Mn51_derived);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_Mn51_to_Fe52_reaclib);
        amrex::Real drdT_pY = 0.0_rt;
        drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa + drdT_pY);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Fe52_He4_to_Ni56_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_pg = rate_eval.screened_rates(k_p_Co55_to_Ni56_reaclib);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_Co55_to_He4_Fe52_derived);
    amrex::Real r_pY = 0.0_rt;
    amrex::Real dd = 1.0_rt / (r_pg + r_pa + r_pY);
    amrex::Real r_ag = rate_eval.screened_rates(k_He4_Fe52_to_Ni56_reaclib);
    amrex::Real r_ap = rate_eval.screened_rates(k_He4_Fe52_to_p_Co55_reaclib);
    rate = r_ag + r_ap * r_pg * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ag = rate_eval.dscreened_rates_dT(k_He4_Fe52_to_Ni56_reaclib);
        amrex::Real drdT_ap = rate_eval.dscreened_rates_dT(k_He4_Fe52_to_p_Co55_reaclib);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_Co55_to_Ni56_reaclib);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_Co55_to_He4_Fe52_derived);
        amrex::Real drdT_pY = 0.0_rt;
        drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa + drdT_pY);
    }
}

template <typename T>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void rate_Ni56_to_Fe52_He4_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

    amrex::Real r_pg = rate_eval.screened_rates(k_p_Co55_to_Ni56_reaclib);
    amrex::Real r_pa = rate_eval.screened_rates(k_p_Co55_to_He4_Fe52_derived);
    amrex::Real r_pY = 0.0_rt;
    amrex::Real dd = 1.0_rt / (r_pg + r_pa + r_pY);
    amrex::Real r_ga = rate_eval.screened_rates(k_Ni56_to_He4_Fe52_derived);
    amrex::Real r_gp = rate_eval.screened_rates(k_Ni56_to_p_Co55_derived);
    rate = r_ga + r_gp * r_pa * dd;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        amrex::Real drdT_ga = rate_eval.dscreened_rates_dT(k_Ni56_to_He4_Fe52_derived);
        amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_p_Co55_to_He4_Fe52_derived);
        amrex::Real drdT_gp = rate_eval.dscreened_rates_dT(k_Ni56_to_p_Co55_derived);
        amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_p_Co55_to_Ni56_reaclib);
        amrex::Real drdT_pY = 0.0_rt;
        drate_dT = drdT_ga + drdT_gp * r_pa * dd + r_gp * drdT_pa * dd - r_gp * r_pa * dd * dd * (drdT_pg + drdT_pa + drdT_pY);
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

    rate_Cr48_He4_to_Fe52_approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_Cr48_He4_to_Fe52_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Cr48_He4_to_Fe52_approx) = drate_dT;
    }

    rate_Fe52_to_Cr48_He4_approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_Fe52_to_Cr48_He4_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe52_to_Cr48_He4_approx) = drate_dT;
    }

    rate_Fe52_He4_to_Ni56_approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_Fe52_He4_to_Ni56_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Fe52_He4_to_Ni56_approx) = drate_dT;
    }

    rate_Ni56_to_Fe52_He4_approx<T>(rate_eval, rate, drate_dT);
    rate_eval.screened_rates(k_Ni56_to_Fe52_He4_approx) = rate;
    if constexpr (std::is_same_v<T, rate_derivs_t>) {
        rate_eval.dscreened_rates_dT(k_Ni56_to_Fe52_He4_approx) = drate_dT;
    }


}

#endif
