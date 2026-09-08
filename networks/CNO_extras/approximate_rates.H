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
    void rate_Ne18_He4_to_Mg22_approx(const T& rate_eval, amrex::Real& rate, amrex::Real& drate_dT) {

        amrex::Real r_pg = rate_eval.screened_rates(k_Na21_p_to_Mg22_reaclib);
        amrex::Real r_pa = rate_eval.screened_rates(k_Na21_p_to_He4_Ne18_derived);
        amrex::Real r_pY = 0.0_rt;
        amrex::Real dd = 1.0_rt / (r_pg + r_pa + r_pY);
        amrex::Real r_ag = rate_eval.screened_rates(k_Ne18_He4_to_Mg22_reaclib);
        amrex::Real r_ap = rate_eval.screened_rates(k_Ne18_He4_to_p_Na21_reaclib);
        rate = r_ag + r_ap * r_pg * dd;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            amrex::Real drdT_ag = rate_eval.dscreened_rates_dT(k_Ne18_He4_to_Mg22_reaclib);
            amrex::Real drdT_ap = rate_eval.dscreened_rates_dT(k_Ne18_He4_to_p_Na21_reaclib);
            amrex::Real drdT_pg = rate_eval.dscreened_rates_dT(k_Na21_p_to_Mg22_reaclib);
            amrex::Real drdT_pa = rate_eval.dscreened_rates_dT(k_Na21_p_to_He4_Ne18_derived);
            amrex::Real drdT_pY = 0.0_rt;
            drate_dT = drdT_ag + drdT_ap * r_pg * dd + r_ap * drdT_pg * dd - r_ap * r_pg * dd * dd * (drdT_pg + drdT_pa + drdT_pY);
        }
    }



    template <int do_T_derivatives, typename T>
    AMREX_GPU_DEVICE AMREX_INLINE
    void
    fill_rates([[maybe_unused]] const tf_t& tfactors,
               [[maybe_unused]] const amrex::Real rho,
               [[maybe_unused]] const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
               [[maybe_unused]] T& rate_eval)
    {

        [[maybe_unused]] amrex::Real rate{};
        [[maybe_unused]] amrex::Real drate_dT{};

        rate_Ne18_He4_to_Mg22_approx<T>(rate_eval, rate, drate_dT);
        rate_eval.screened_rates(k_Ne18_He4_to_Mg22_approx) = rate;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dscreened_rates_dT(k_Ne18_He4_to_Mg22_approx) = drate_dT;
        }


    }
}
#endif
