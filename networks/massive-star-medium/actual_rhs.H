#ifndef ACTUAL_RHS_H
#define ACTUAL_RHS_H

#include <AMReX_REAL.H>
#include <AMReX_Array.H>

#include <extern_parameters.H>
#include <actual_network.H>
#include <burn_type.H>
#include <jacobian_utilities.H>
#ifdef SCREENING
#include <screen.H>
#endif
#include <microphysics_autodiff.H>
#ifdef NEUTRINOS
#include <neutrino.H>
#endif
#include <rate_type.H>
#include <reaclib_rates.H>
#include <approximate_rates.H>
#include <modified_rates.H>
#include <branched_rates.H>
#include <beta_limited_rates.H>
#include <temperature_table_rates.H>
#include <tabular_weak_rates.H>
#include <derived_rates.H>
#include <fill_rates.H>

using namespace amrex::literals;

template<class T>
AMREX_GPU_DEVICE AMREX_INLINE
void ener_gener_rate(T const& dydt, amrex::Real& enuc)
{

    // Computes the instantaneous energy generation rate (from the nuclei)

    // This is basically e = m c**2

    enuc = 0.0_rt;

    amrex::constexpr_for<1, NumSpec+1>([&] (auto n)
    {
        enuc += dydt(n) * network::mion<n>();
    });

    enuc *= C::enuc_conv2;
}


#ifdef SCREENING
template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void evaluate_screening(const burn_t& state,
                        const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
                        T& rate_eval) {

    using namespace Rates;

    // Evaluate screening terms

    using number_t = std::conditional_t<do_T_derivatives, autodiff::dual, amrex::Real>;
    number_t temp = state.T;
    if constexpr (do_T_derivatives) {
        // seed the dual number for temperature before calculating anything with it
        autodiff::seed(temp);
    }

    plasma_state_t<number_t> pstate{};
    fill_plasma_state(pstate, temp, state.rho, Y);

    amrex::Real log_scor, dlog_scor_dT;

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 15.0_rt, 33.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_P33) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_P33) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 25.0_rt, 53.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Mn53) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Mn53) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 14.0_rt, 31.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Si31) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Si31) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 22.0_rt, 44.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Ti44) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Ti44) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 22.0_rt, 48.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Ti48) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Ti48) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 20.0_rt, 42.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Ca42) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Ca42) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 21.0_rt, 43.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Sc43) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Sc43) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 26.0_rt, 56.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Fe56) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Fe56) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 18.0_rt, 38.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Ar38) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Ar38) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 24.0_rt, 54.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Cr54) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Cr54) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 22.0_rt, 44.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Ti44) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Ti44) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 27.0_rt, 57.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Co57) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Co57) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 20.0_rt, 42.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Ca42) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Ca42) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 12.0_rt, 24.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Mg24) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Mg24) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(6.0_rt, 12.0_rt, 8.0_rt, 16.0_rt);
        static_assert(scn_fac.z1 == 6.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_C12_O16) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_C12_O16) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 18.0_rt, 38.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Ar38) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Ar38) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 23.0_rt, 47.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_V47) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_V47) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 23.0_rt, 51.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_V51) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_V51) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 21.0_rt, 45.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Sc45) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Sc45) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 27.0_rt, 55.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Co55) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Co55) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(6.0_rt, 12.0_rt, 6.0_rt, 12.0_rt);
        static_assert(scn_fac.z1 == 6.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_C12_C12) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_C12_C12) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 17.0_rt, 33.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Cl33) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Cl33) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 13.0_rt, 27.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Al27) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Al27) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 22.0_rt, 46.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Ti46) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Ti46) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 23.0_rt, 47.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_V47) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_V47) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 21.0_rt, 45.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Sc45) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Sc45) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 28.0_rt, 56.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Ni56) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Ni56) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 15.0_rt, 31.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_P31) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_P31) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 28.0_rt, 58.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Ni58) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Ni58) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 14.0_rt, 29.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Si29) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Si29) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 23.0_rt, 49.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_V49) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_V49) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(8.0_rt, 16.0_rt, 8.0_rt, 16.0_rt);
        static_assert(scn_fac.z1 == 8.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_O16_O16) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_O16_O16) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 15.0_rt, 31.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_P31) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_P31) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 14.0_rt, 29.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Si29) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Si29) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 24.0_rt, 52.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Cr52) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Cr52) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 16.0_rt, 34.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_S34) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_S34) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 18.0_rt, 37.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Ar37) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Ar37) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 17.0_rt, 37.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Cl37) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Cl37) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 20.0_rt, 40.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Ca40) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Ca40) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 18.0_rt, 36.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Ar36) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Ar36) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 15.0_rt, 33.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_P33) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_P33) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 24.0_rt, 52.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Cr52) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Cr52) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 16.0_rt, 34.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_S34) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_S34) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 14.0_rt, 31.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Si31) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Si31) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 25.0_rt, 55.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Mn55) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Mn55) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 21.0_rt, 43.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Sc43) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Sc43) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 26.0_rt, 56.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Fe56) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Fe56) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 24.0_rt, 54.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Cr54) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Cr54) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 25.0_rt, 55.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Mn55) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Mn55) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 10.0_rt, 20.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Ne20) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Ne20) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 27.0_rt, 56.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Co56) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Co56) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 14.0_rt, 28.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Si28) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Si28) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 26.0_rt, 54.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Fe54) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Fe54) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 27.0_rt, 55.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Co55) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Co55) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 13.0_rt, 27.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Al27) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Al27) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 17.0_rt, 35.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Cl35) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Cl35) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 16.0_rt, 33.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_S33) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_S33) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 14.0_rt, 30.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Si30) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Si30) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 8.0_rt, 16.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_O16) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_O16) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 14.0_rt, 30.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Si30) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Si30) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 16.0_rt, 35.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_S35) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_S35) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 12.0_rt, 23.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Mg23) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Mg23) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 18.0_rt, 37.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Ar37) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Ar37) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 27.0_rt, 59.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Co59) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Co59) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 16.0_rt, 35.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_S35) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_S35) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 15.0_rt, 29.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_P29) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_P29) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 24.0_rt, 50.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Cr50) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Cr50) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 16.0_rt, 32.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_S32) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_S32) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 29.0_rt, 59.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Cu59) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Cu59) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 16.0_rt, 32.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_S32) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_S32) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 25.0_rt, 53.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Mn53) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Mn53) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 21.0_rt, 41.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Sc41) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Sc41) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 27.0_rt, 56.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Co56) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Co56) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 14.0_rt, 28.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Si28) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Si28) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 26.0_rt, 54.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Fe54) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Fe54) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 22.0_rt, 48.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Ti48) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Ti48) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 4.0_rt, 8.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Be8) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Be8) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 27.0_rt, 57.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Co57) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Co57) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 19.0_rt, 39.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_K39) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_K39) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 17.0_rt, 35.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Cl35) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Cl35) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 16.0_rt, 33.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_S33) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_S33) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 23.0_rt, 51.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_V51) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_V51) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 6.0_rt, 12.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_C12) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_C12) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 19.0_rt, 39.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_K39) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_K39) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 22.0_rt, 46.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Ti46) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Ti46) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 28.0_rt, 60.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Ni60) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Ni60) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 26.0_rt, 58.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Fe58) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Fe58) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 24.0_rt, 48.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Cr48) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Cr48) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 25.0_rt, 51.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Mn51) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Mn51) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 11.0_rt, 23.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Na23) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Na23) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 23.0_rt, 49.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_V49) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_V49) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 26.0_rt, 52.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Fe52) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Fe52) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 19.0_rt, 37.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_K37) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_K37) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 24.0_rt, 50.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Cr50) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Cr50) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 25.0_rt, 51.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_Mn51) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_Mn51) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(2.0_rt, 4.0_rt, 2.0_rt, 4.0_rt);
        static_assert(scn_fac.z1 == 2.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_He4_He4) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_He4_He4) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(6.0_rt, 12.0_rt, 10.0_rt, 20.0_rt);
        static_assert(scn_fac.z1 == 6.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_C12_Ne20) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_C12_Ne20) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 20.0_rt, 40.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Ca40) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Ca40) = dlog_scor_dT;
        }
    }

    {
        constexpr auto scn_fac = scrn::calculate_screen_factor(1.0_rt, 1.0_rt, 18.0_rt, 36.0_rt);
        static_assert(scn_fac.z1 == 1.0_rt);
        actual_log_screen(pstate, scn_fac, log_scor, dlog_scor_dT);
        rate_eval.log_screen(k_p_Ar36) = log_scor;
        if constexpr (std::is_same_v<T, rate_derivs_t>) {
            rate_eval.dlog_screen_dT(k_p_Ar36) = dlog_scor_dT;
        }
    }

}
#endif


template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void evaluate_rates(const burn_t& state,
                    const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
                    T& rate_eval) {

    using namespace Rates;

    // here Y is consistent with state.xn

    [[maybe_unused]] const amrex::Real rhoy = state.rho * state.y_e;

    const tf_t tfactors = evaluate_tfactors(state.T);

    // Precompute screening terms

#ifdef SCREENING
    evaluate_screening<do_T_derivatives, T>(state, Y, rate_eval);
#endif

    // Fill in different rates

    rate_util::fill_rates(tfactors, state.T, state.rho, rhoy, Y, rate_eval);

    // Apply any parameters to disable rates


}


// backwards-compatible wrapper
template <int do_T_derivatives, typename T>
AMREX_GPU_DEVICE AMREX_INLINE
void evaluate_rates(const burn_t& state,
                    T& rate_eval) {

    // create molar fractions
    amrex::Array1D<amrex::Real, 1, NumSpec> Y;
    for (int n = 1; n <= NumSpec; ++n) {
        Y(n) = state.xn[n-1] * aion_inv[n-1];
    }

    evaluate_rates<do_T_derivatives>(state, Y, rate_eval);
}


#ifdef NSE_NET
AMREX_GPU_DEVICE AMREX_INLINE
void get_ydot_weak(const burn_t& state,
                   amrex::Array1D<amrex::Real, 1, neqs>& ydot_nuc,
                   amrex::Real& enuc_weak,
                   [[maybe_unused]] const amrex::Array1D<amrex::Real, 1, NumSpec>& Y) {
    ///
    /// Calculate Ydots contribute only from weak reactions.
    /// This is used to calculate dyedt and energy generation from
    /// weak reactions for self-consistent NSE
    ///

    using namespace modified_rates;
    using namespace temp_tabular;

    // initialize ydot_nuc to 0

    for (int i = 1; i <= neqs; ++i) {
        ydot_nuc(i) = 0.0_rt;
    }

    rate_t rate_eval;

    [[maybe_unused]] const amrex::Real rhoy = state.rho * state.y_e;

    const tf_t tfactors = evaluate_tfactors(state.T);

    // fill only the weak rates
    rate_util::fill_weak_rates(tfactors, state.T, state.rho, rhoy, Y, rate_eval);

    // assemble Ydot for just the weak rate evolution

    const auto& screened_rates = rate_eval.screened_rates;

    ydot_nuc(N) =
        (-screened_rates(k_n_to_p_beta_neg_weaktab)*Y(N) + screened_rates(k_p_to_n_electron_capture_weaktab)*Y(H1));

    ydot_nuc(H1) =
        (screened_rates(k_n_to_p_beta_neg_weaktab)*Y(N) + -screened_rates(k_p_to_n_electron_capture_weaktab)*Y(H1));

    ydot_nuc(He4) = 0.0_rt;

    ydot_nuc(C12) = 0.0_rt;

    ydot_nuc(O16) = 0.0_rt;

    ydot_nuc(Ne20) = 0.0_rt;

    ydot_nuc(Mg23) = 0.0_rt;

    ydot_nuc(Mg24) = 0.0_rt;

    ydot_nuc(Al27) =
        (screened_rates(k_Si27_to_Al27_electron_capture_weaktab)*Y(Si27) + -screened_rates(k_Al27_to_Si27_beta_neg_weaktab)*Y(Al27));

    ydot_nuc(Si27) =
        (-screened_rates(k_Si27_to_Al27_electron_capture_weaktab)*Y(Si27) + screened_rates(k_Al27_to_Si27_beta_neg_weaktab)*Y(Al27));

    ydot_nuc(Si28) = 0.0_rt;

    ydot_nuc(Si29) =
        (screened_rates(k_P29_to_Si29_electron_capture_weaktab)*Y(P29) + -screened_rates(k_Si29_to_P29_beta_neg_weaktab)*Y(Si29));

    ydot_nuc(Si30) = 0.0_rt;

    ydot_nuc(Si31) =
        (-screened_rates(k_Si31_to_P31_beta_neg_weaktab)*Y(Si31) + screened_rates(k_P31_to_Si31_electron_capture_weaktab)*Y(P31));

    ydot_nuc(P29) =
        (-screened_rates(k_P29_to_Si29_electron_capture_weaktab)*Y(P29) + screened_rates(k_Si29_to_P29_beta_neg_weaktab)*Y(Si29));

    ydot_nuc(P31) =
        (screened_rates(k_Si31_to_P31_beta_neg_weaktab)*Y(Si31) + -screened_rates(k_P31_to_Si31_electron_capture_weaktab)*Y(P31));

    ydot_nuc(P33) =
        (-screened_rates(k_P33_to_S33_beta_neg_weaktab)*Y(P33) + screened_rates(k_S33_to_P33_electron_capture_weaktab)*Y(S33));

    ydot_nuc(S32) = 0.0_rt;

    ydot_nuc(S33) =
        (screened_rates(k_Cl33_to_S33_electron_capture_weaktab)*Y(Cl33) + -screened_rates(k_S33_to_Cl33_beta_neg_weaktab)*Y(S33)) +
        (screened_rates(k_P33_to_S33_beta_neg_weaktab)*Y(P33) + -screened_rates(k_S33_to_P33_electron_capture_weaktab)*Y(S33));

    ydot_nuc(S34) = 0.0_rt;

    ydot_nuc(S35) =
        (-screened_rates(k_S35_to_Cl35_beta_neg_weaktab)*Y(S35) + screened_rates(k_Cl35_to_S35_electron_capture_weaktab)*Y(Cl35));

    ydot_nuc(Cl33) =
        (-screened_rates(k_Cl33_to_S33_electron_capture_weaktab)*Y(Cl33) + screened_rates(k_S33_to_Cl33_beta_neg_weaktab)*Y(S33));

    ydot_nuc(Cl35) =
        (screened_rates(k_S35_to_Cl35_beta_neg_weaktab)*Y(S35) + -screened_rates(k_Cl35_to_S35_electron_capture_weaktab)*Y(Cl35));

    ydot_nuc(Cl37) =
        (screened_rates(k_Ar37_to_Cl37_electron_capture_weaktab)*Y(Ar37) + -screened_rates(k_Cl37_to_Ar37_beta_neg_weaktab)*Y(Cl37));

    ydot_nuc(Ar36) = 0.0_rt;

    ydot_nuc(Ar37) =
        (-screened_rates(k_Ar37_to_Cl37_electron_capture_weaktab)*Y(Ar37) + screened_rates(k_Cl37_to_Ar37_beta_neg_weaktab)*Y(Cl37)) +
        (screened_rates(k_K37_to_Ar37_electron_capture_weaktab)*Y(K37) + -screened_rates(k_Ar37_to_K37_beta_neg_weaktab)*Y(Ar37));

    ydot_nuc(Ar38) = 0.0_rt;

    ydot_nuc(K37) =
        (-screened_rates(k_K37_to_Ar37_electron_capture_weaktab)*Y(K37) + screened_rates(k_Ar37_to_K37_beta_neg_weaktab)*Y(Ar37));

    ydot_nuc(K39) = 0.0_rt;

    ydot_nuc(Ca40) = 0.0_rt;

    ydot_nuc(Ca42) = 0.0_rt;

    ydot_nuc(Sc41) = 0.0_rt;

    ydot_nuc(Sc43) = 0.0_rt;

    ydot_nuc(Sc45) = 0.0_rt;

    ydot_nuc(Ti44) = 0.0_rt;

    ydot_nuc(Ti46) = 0.0_rt;

    ydot_nuc(Ti48) = 0.0_rt;

    ydot_nuc(V47) = 0.0_rt;

    ydot_nuc(V49) = 0.0_rt;

    ydot_nuc(V51) = 0.0_rt;

    ydot_nuc(Cr48) = 0.0_rt;

    ydot_nuc(Cr50) = 0.0_rt;

    ydot_nuc(Cr52) = 0.0_rt;

    ydot_nuc(Cr54) = 0.0_rt;

    ydot_nuc(Mn51) = 0.0_rt;

    ydot_nuc(Mn53) = 0.0_rt;

    ydot_nuc(Mn55) = 0.0_rt;

    ydot_nuc(Fe52) = 0.0_rt;

    ydot_nuc(Fe54) = 0.0_rt;

    ydot_nuc(Fe56) =
        (screened_rates(k_Co56_to_Fe56_electron_capture_weaktab)*Y(Co56) + -screened_rates(k_Fe56_to_Co56_beta_neg_weaktab)*Y(Fe56));

    ydot_nuc(Fe58) = 0.0_rt;

    ydot_nuc(Co55) = 0.0_rt;

    ydot_nuc(Co56) =
        (-screened_rates(k_Co56_to_Fe56_electron_capture_weaktab)*Y(Co56) + screened_rates(k_Fe56_to_Co56_beta_neg_weaktab)*Y(Fe56)) +
        (screened_rates(k_Ni56_to_Co56_electron_capture_weaktab)*Y(Ni56) + -screened_rates(k_Co56_to_Ni56_beta_neg_weaktab)*Y(Co56));

    ydot_nuc(Co57) = 0.0_rt;

    ydot_nuc(Co59) = 0.0_rt;

    ydot_nuc(Ni56) =
        (-screened_rates(k_Ni56_to_Co56_electron_capture_weaktab)*Y(Ni56) + screened_rates(k_Co56_to_Ni56_beta_neg_weaktab)*Y(Co56));

    ydot_nuc(Ni58) = 0.0_rt;

    ydot_nuc(Ni60) = 0.0_rt;

    ydot_nuc(Cu59) = 0.0_rt;

    ydot_nuc(Zn60) = 0.0_rt;


    enuc_weak = rate_eval.enuc_weak;
}
#endif


AMREX_GPU_DEVICE AMREX_INLINE
void rhs_nuc(const burn_t& state,
             amrex::Array1D<amrex::Real, 1, neqs>& ydot_nuc,
             const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
             const amrex::Array1D<amrex::Real, 1, Rates::NumRates>& screened_rates) {

    using namespace Rates;

    ydot_nuc(N) =
        (-screened_rates(k_Mg23_n_to_Mg24_reaclib)*Y(Mg23)*Y(N)*state.rho + screened_rates(k_Mg24_to_n_Mg23_reaclib)*Y(Mg24)) +
        (-screened_rates(k_Si27_n_to_Si28_reaclib)*Y(Si27)*Y(N)*state.rho + screened_rates(k_Si28_to_n_Si27_reaclib)*Y(Si28)) +
        (-screened_rates(k_Si28_n_to_Si29_reaclib)*Y(Si28)*Y(N)*state.rho + screened_rates(k_Si29_to_n_Si28_reaclib)*Y(Si29)) +
        (-screened_rates(k_Si29_n_to_Si30_reaclib)*Y(Si29)*Y(N)*state.rho + screened_rates(k_Si30_to_n_Si29_reaclib)*Y(Si30)) +
        (-screened_rates(k_Si30_n_to_Si31_reaclib)*Y(Si30)*Y(N)*state.rho + screened_rates(k_Si31_to_n_Si30_reaclib)*Y(Si31)) +
        (-screened_rates(k_S32_n_to_S33_reaclib)*Y(S32)*Y(N)*state.rho + screened_rates(k_S33_to_n_S32_reaclib)*Y(S33)) +
        (-screened_rates(k_S33_n_to_S34_reaclib)*Y(S33)*Y(N)*state.rho + screened_rates(k_S34_to_n_S33_reaclib)*Y(S34)) +
        (-screened_rates(k_S34_n_to_S35_reaclib)*Y(S34)*Y(N)*state.rho + screened_rates(k_S35_to_n_S34_reaclib)*Y(S35)) +
        (-screened_rates(k_Ar36_n_to_Ar37_reaclib)*Y(Ar36)*Y(N)*state.rho + screened_rates(k_Ar37_to_n_Ar36_reaclib)*Y(Ar37)) +
        (-screened_rates(k_Ar37_n_to_Ar38_reaclib)*Y(Ar37)*Y(N)*state.rho + screened_rates(k_Ar38_to_n_Ar37_reaclib)*Y(Ar38)) +
        (-screened_rates(k_Co55_n_to_Co56_reaclib)*Y(Co55)*Y(N)*state.rho + screened_rates(k_Co56_to_n_Co55_reaclib)*Y(Co56)) +
        (-screened_rates(k_Co56_n_to_Co57_reaclib)*Y(Co56)*Y(N)*state.rho + screened_rates(k_Co57_to_n_Co56_reaclib)*Y(Co57)) +
        (-screened_rates(k_Mg23_n_to_He4_Ne20_reaclib)*Y(Mg23)*Y(N)*state.rho + screened_rates(k_Ne20_He4_to_n_Mg23_reaclib)*Y(He4)*Y(Ne20)*state.rho) +
        (-screened_rates(k_Mg23_n_to_C12_C12_reaclib)*Y(Mg23)*Y(N)*state.rho + 0.5*screened_rates(k_C12_C12_to_n_Mg23_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho) +
        (-screened_rates(k_Si27_n_to_p_Al27_reaclib)*Y(Si27)*Y(N)*state.rho + screened_rates(k_Al27_p_to_n_Si27_reaclib)*Y(Al27)*Y(H1)*state.rho) +
        (-screened_rates(k_Si27_n_to_He4_Mg24_reaclib)*Y(Si27)*Y(N)*state.rho + screened_rates(k_Mg24_He4_to_n_Si27_reaclib)*Y(He4)*Y(Mg24)*state.rho) +
        (-screened_rates(k_Si27_n_to_C12_O16_reaclib)*Y(Si27)*Y(N)*state.rho + screened_rates(k_O16_C12_to_n_Si27_reaclib)*Y(C12)*Y(O16)*state.rho) +
        (screened_rates(k_Si31_p_to_n_P31_reaclib)*Y(Si31)*Y(H1)*state.rho + -screened_rates(k_P31_n_to_p_Si31_reaclib)*Y(P31)*Y(N)*state.rho) +
        (screened_rates(k_Si31_He4_to_n_S34_reaclib)*Y(He4)*Y(Si31)*state.rho + -screened_rates(k_S34_n_to_He4_Si31_reaclib)*Y(S34)*Y(N)*state.rho) +
        (-screened_rates(k_P29_n_to_p_Si29_reaclib)*Y(P29)*Y(N)*state.rho + screened_rates(k_Si29_p_to_n_P29_reaclib)*Y(Si29)*Y(H1)*state.rho) +
        (-screened_rates(k_S32_n_to_He4_Si29_reaclib)*Y(S32)*Y(N)*state.rho + screened_rates(k_Si29_He4_to_n_S32_reaclib)*Y(He4)*Y(Si29)*state.rho) +
        (-screened_rates(k_S33_n_to_p_P33_reaclib)*Y(S33)*Y(N)*state.rho + screened_rates(k_P33_p_to_n_S33_reaclib)*Y(P33)*Y(H1)*state.rho) +
        (-screened_rates(k_S33_n_to_He4_Si30_reaclib)*Y(S33)*Y(N)*state.rho + screened_rates(k_Si30_He4_to_n_S33_reaclib)*Y(He4)*Y(Si30)*state.rho) +
        (screened_rates(k_S35_He4_to_n_Ar38_reaclib)*Y(He4)*Y(S35)*state.rho + -screened_rates(k_Ar38_n_to_He4_S35_reaclib)*Y(Ar38)*Y(N)*state.rho) +
        (-screened_rates(k_Cl33_n_to_p_S33_reaclib)*Y(Cl33)*Y(N)*state.rho + screened_rates(k_S33_p_to_n_Cl33_reaclib)*Y(S33)*Y(H1)*state.rho) +
        (-screened_rates(k_Cl35_n_to_p_S35_reaclib)*Y(Cl35)*Y(N)*state.rho + screened_rates(k_S35_p_to_n_Cl35_reaclib)*Y(S35)*Y(H1)*state.rho) +
        (-screened_rates(k_Ar36_n_to_He4_S33_reaclib)*Y(Ar36)*Y(N)*state.rho + screened_rates(k_S33_He4_to_n_Ar36_reaclib)*Y(He4)*Y(S33)*state.rho) +
        (-screened_rates(k_Ar37_n_to_p_Cl37_reaclib)*Y(Ar37)*Y(N)*state.rho + screened_rates(k_Cl37_p_to_n_Ar37_reaclib)*Y(Cl37)*Y(H1)*state.rho) +
        (-screened_rates(k_Ar37_n_to_He4_S34_reaclib)*Y(Ar37)*Y(N)*state.rho + screened_rates(k_S34_He4_to_n_Ar37_reaclib)*Y(He4)*Y(S34)*state.rho) +
        (-screened_rates(k_K37_n_to_p_Ar37_reaclib)*Y(K37)*Y(N)*state.rho + screened_rates(k_Ar37_p_to_n_K37_reaclib)*Y(Ar37)*Y(H1)*state.rho) +
        (-screened_rates(k_Ca40_n_to_He4_Ar37_reaclib)*Y(Ca40)*Y(N)*state.rho + screened_rates(k_Ar37_He4_to_n_Ca40_reaclib)*Y(Ar37)*Y(He4)*state.rho) +
        (-screened_rates(k_Co56_n_to_p_Fe56_reaclib)*Y(Co56)*Y(N)*state.rho + screened_rates(k_Fe56_p_to_n_Co56_reaclib)*Y(Fe56)*Y(H1)*state.rho) +
        (-screened_rates(k_Co56_n_to_He4_Mn53_reaclib)*Y(Co56)*Y(N)*state.rho + screened_rates(k_Mn53_He4_to_n_Co56_reaclib)*Y(He4)*Y(Mn53)*state.rho) +
        (-screened_rates(k_Ni56_n_to_p_Co56_reaclib)*Y(Ni56)*Y(N)*state.rho + screened_rates(k_Co56_p_to_n_Ni56_reaclib)*Y(Co56)*Y(H1)*state.rho) +
        (-screened_rates(k_Cu59_n_to_He4_Co56_reaclib)*Y(Cu59)*Y(N)*state.rho + screened_rates(k_Co56_He4_to_n_Cu59_reaclib)*Y(Co56)*Y(He4)*state.rho) +
        (-screened_rates(k_Ni56_n_n_to_Ni58_approx)*Y(Ni56)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + 2.0*screened_rates(k_Ni58_to_Ni56_n_n_approx)*Y(Ni58)) +
        (-screened_rates(k_Ni58_n_n_to_Ni60_approx)*Y(Ni58)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + 2.0*screened_rates(k_Ni60_to_Ni58_n_n_approx)*Y(Ni60)) +
        (-screened_rates(k_Co57_n_n_to_Co59_approx)*Y(Co57)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + 2.0*screened_rates(k_Co59_to_Co57_n_n_approx)*Y(Co59)) +
        (-screened_rates(k_Fe52_n_n_to_Fe54_approx)*Y(Fe52)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + 2.0*screened_rates(k_Fe54_to_Fe52_n_n_approx)*Y(Fe54)) +
        (-screened_rates(k_Fe54_n_n_to_Fe56_approx)*Y(Fe54)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + 2.0*screened_rates(k_Fe56_to_Fe54_n_n_approx)*Y(Fe56)) +
        (-screened_rates(k_Fe56_n_n_to_Fe58_approx)*Y(Fe56)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + 2.0*screened_rates(k_Fe58_to_Fe56_n_n_approx)*Y(Fe58)) +
        (-screened_rates(k_Mn51_n_n_to_Mn53_approx)*Y(Mn51)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + 2.0*screened_rates(k_Mn53_to_Mn51_n_n_approx)*Y(Mn53)) +
        (-screened_rates(k_Mn53_n_n_to_Mn55_approx)*Y(Mn53)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + 2.0*screened_rates(k_Mn55_to_Mn53_n_n_approx)*Y(Mn55)) +
        (-screened_rates(k_Cr48_n_n_to_Cr50_approx)*Y(Cr48)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + 2.0*screened_rates(k_Cr50_to_Cr48_n_n_approx)*Y(Cr50)) +
        (-screened_rates(k_Cr50_n_n_to_Cr52_approx)*Y(Cr50)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + 2.0*screened_rates(k_Cr52_to_Cr50_n_n_approx)*Y(Cr52)) +
        (-screened_rates(k_Cr52_n_n_to_Cr54_approx)*Y(Cr52)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + 2.0*screened_rates(k_Cr54_to_Cr52_n_n_approx)*Y(Cr54)) +
        (-screened_rates(k_V47_n_n_to_V49_approx)*Y(V47)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + 2.0*screened_rates(k_V49_to_V47_n_n_approx)*Y(V49)) +
        (-screened_rates(k_V49_n_n_to_V51_approx)*Y(V49)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + 2.0*screened_rates(k_V51_to_V49_n_n_approx)*Y(V51)) +
        (-screened_rates(k_Ti44_n_n_to_Ti46_approx)*Y(Ti44)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + 2.0*screened_rates(k_Ti46_to_Ti44_n_n_approx)*Y(Ti46)) +
        (-screened_rates(k_Ti46_n_n_to_Ti48_approx)*Y(Ti46)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + 2.0*screened_rates(k_Ti48_to_Ti46_n_n_approx)*Y(Ti48)) +
        (-screened_rates(k_Sc41_n_n_to_Sc43_approx)*Y(Sc41)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + 2.0*screened_rates(k_Sc43_to_Sc41_n_n_approx)*Y(Sc43)) +
        (-screened_rates(k_Sc43_n_n_to_Sc45_approx)*Y(Sc43)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + 2.0*screened_rates(k_Sc45_to_Sc43_n_n_approx)*Y(Sc45)) +
        (-screened_rates(k_Ca40_n_n_to_Ca42_approx)*Y(Ca40)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + 2.0*screened_rates(k_Ca42_to_Ca40_n_n_approx)*Y(Ca42)) +
        (-screened_rates(k_n_to_p_beta_neg_weaktab)*Y(N) + screened_rates(k_p_to_n_electron_capture_weaktab)*Y(H1));

    ydot_nuc(H1) =
        (-screened_rates(k_Al27_p_to_Si28_reaclib)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_Si28_to_p_Al27_reaclib)*Y(Si28)) +
        (-screened_rates(k_Si28_p_to_P29_reaclib)*Y(Si28)*Y(H1)*state.rho + screened_rates(k_P29_to_p_Si28_reaclib)*Y(P29)) +
        (-screened_rates(k_Si30_p_to_P31_reaclib)*Y(Si30)*Y(H1)*state.rho + screened_rates(k_P31_to_p_Si30_reaclib)*Y(P31)) +
        (-screened_rates(k_P31_p_to_S32_reaclib)*Y(P31)*Y(H1)*state.rho + screened_rates(k_S32_to_p_P31_reaclib)*Y(S32)) +
        (-screened_rates(k_P33_p_to_S34_reaclib)*Y(P33)*Y(H1)*state.rho + screened_rates(k_S34_to_p_P33_reaclib)*Y(S34)) +
        (-screened_rates(k_S32_p_to_Cl33_reaclib)*Y(S32)*Y(H1)*state.rho + screened_rates(k_Cl33_to_p_S32_reaclib)*Y(Cl33)) +
        (-screened_rates(k_S34_p_to_Cl35_reaclib)*Y(S34)*Y(H1)*state.rho + screened_rates(k_Cl35_to_p_S34_reaclib)*Y(Cl35)) +
        (-screened_rates(k_Cl35_p_to_Ar36_reaclib)*Y(Cl35)*Y(H1)*state.rho + screened_rates(k_Ar36_to_p_Cl35_reaclib)*Y(Ar36)) +
        (-screened_rates(k_Cl37_p_to_Ar38_reaclib)*Y(Cl37)*Y(H1)*state.rho + screened_rates(k_Ar38_to_p_Cl37_reaclib)*Y(Ar38)) +
        (-screened_rates(k_Ar36_p_to_K37_reaclib)*Y(Ar36)*Y(H1)*state.rho + screened_rates(k_K37_to_p_Ar36_reaclib)*Y(K37)) +
        (-screened_rates(k_Ar38_p_to_K39_reaclib)*Y(Ar38)*Y(H1)*state.rho + screened_rates(k_K39_to_p_Ar38_reaclib)*Y(K39)) +
        (-screened_rates(k_K39_p_to_Ca40_reaclib)*Y(K39)*Y(H1)*state.rho + screened_rates(k_Ca40_to_p_K39_reaclib)*Y(Ca40)) +
        (-screened_rates(k_Ca40_p_to_Sc41_reaclib)*Y(Ca40)*Y(H1)*state.rho + screened_rates(k_Sc41_to_p_Ca40_reaclib)*Y(Sc41)) +
        (-screened_rates(k_Ca42_p_to_Sc43_reaclib)*Y(Ca42)*Y(H1)*state.rho + screened_rates(k_Sc43_to_p_Ca42_reaclib)*Y(Sc43)) +
        (-screened_rates(k_Sc43_p_to_Ti44_reaclib)*Y(Sc43)*Y(H1)*state.rho + screened_rates(k_Ti44_to_p_Sc43_reaclib)*Y(Ti44)) +
        (-screened_rates(k_Sc45_p_to_Ti46_reaclib)*Y(Sc45)*Y(H1)*state.rho + screened_rates(k_Ti46_to_p_Sc45_reaclib)*Y(Ti46)) +
        (-screened_rates(k_Ti46_p_to_V47_reaclib)*Y(Ti46)*Y(H1)*state.rho + screened_rates(k_V47_to_p_Ti46_reaclib)*Y(V47)) +
        (-screened_rates(k_Ti48_p_to_V49_reaclib)*Y(Ti48)*Y(H1)*state.rho + screened_rates(k_V49_to_p_Ti48_reaclib)*Y(V49)) +
        (-screened_rates(k_V47_p_to_Cr48_reaclib)*Y(V47)*Y(H1)*state.rho + screened_rates(k_Cr48_to_p_V47_reaclib)*Y(Cr48)) +
        (-screened_rates(k_V49_p_to_Cr50_reaclib)*Y(V49)*Y(H1)*state.rho + screened_rates(k_Cr50_to_p_V49_reaclib)*Y(Cr50)) +
        (-screened_rates(k_V51_p_to_Cr52_reaclib)*Y(V51)*Y(H1)*state.rho + screened_rates(k_Cr52_to_p_V51_reaclib)*Y(Cr52)) +
        (-screened_rates(k_Cr50_p_to_Mn51_reaclib)*Y(Cr50)*Y(H1)*state.rho + screened_rates(k_Mn51_to_p_Cr50_reaclib)*Y(Mn51)) +
        (-screened_rates(k_Cr52_p_to_Mn53_reaclib)*Y(Cr52)*Y(H1)*state.rho + screened_rates(k_Mn53_to_p_Cr52_reaclib)*Y(Mn53)) +
        (-screened_rates(k_Cr54_p_to_Mn55_reaclib)*Y(Cr54)*Y(H1)*state.rho + screened_rates(k_Mn55_to_p_Cr54_reaclib)*Y(Mn55)) +
        (-screened_rates(k_Mn51_p_to_Fe52_reaclib)*Y(Mn51)*Y(H1)*state.rho + screened_rates(k_Fe52_to_p_Mn51_reaclib)*Y(Fe52)) +
        (-screened_rates(k_Mn53_p_to_Fe54_reaclib)*Y(Mn53)*Y(H1)*state.rho + screened_rates(k_Fe54_to_p_Mn53_reaclib)*Y(Fe54)) +
        (-screened_rates(k_Mn55_p_to_Fe56_reaclib)*Y(Mn55)*Y(H1)*state.rho + screened_rates(k_Fe56_to_p_Mn55_reaclib)*Y(Fe56)) +
        (-screened_rates(k_Fe54_p_to_Co55_reaclib)*Y(Fe54)*Y(H1)*state.rho + screened_rates(k_Co55_to_p_Fe54_reaclib)*Y(Co55)) +
        (-screened_rates(k_Fe56_p_to_Co57_reaclib)*Y(Fe56)*Y(H1)*state.rho + screened_rates(k_Co57_to_p_Fe56_reaclib)*Y(Co57)) +
        (-screened_rates(k_Fe58_p_to_Co59_reaclib)*Y(Fe58)*Y(H1)*state.rho + screened_rates(k_Co59_to_p_Fe58_reaclib)*Y(Co59)) +
        (-screened_rates(k_Co55_p_to_Ni56_reaclib)*Y(Co55)*Y(H1)*state.rho + screened_rates(k_Ni56_to_p_Co55_reaclib)*Y(Ni56)) +
        (-screened_rates(k_Co57_p_to_Ni58_reaclib)*Y(Co57)*Y(H1)*state.rho + screened_rates(k_Ni58_to_p_Co57_reaclib)*Y(Ni58)) +
        (-screened_rates(k_Co59_p_to_Ni60_reaclib)*Y(Co59)*Y(H1)*state.rho + screened_rates(k_Ni60_to_p_Co59_reaclib)*Y(Ni60)) +
        (-screened_rates(k_Ni58_p_to_Cu59_reaclib)*Y(Ni58)*Y(H1)*state.rho + screened_rates(k_Cu59_to_p_Ni58_reaclib)*Y(Cu59)) +
        (-screened_rates(k_Cu59_p_to_Zn60_reaclib)*Y(Cu59)*Y(H1)*state.rho + screened_rates(k_Zn60_to_p_Cu59_reaclib)*Y(Zn60)) +
        (screened_rates(k_O16_C12_to_p_Al27_reaclib)*Y(C12)*Y(O16)*state.rho + -screened_rates(k_Al27_p_to_C12_O16_reaclib)*Y(Al27)*Y(H1)*state.rho) +
        (0.5*screened_rates(k_O16_O16_to_p_P31_reaclib)*amrex::Math::powi<2>(Y(O16))*state.rho + -screened_rates(k_P31_p_to_O16_O16_reaclib)*Y(P31)*Y(H1)*state.rho) +
        (screened_rates(k_Ne20_C12_to_p_P31_reaclib)*Y(C12)*Y(Ne20)*state.rho + -screened_rates(k_P31_p_to_C12_Ne20_reaclib)*Y(P31)*Y(H1)*state.rho) +
        (-screened_rates(k_Al27_p_to_He4_Mg24_reaclib)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_Mg24_He4_to_p_Al27_reaclib)*Y(He4)*Y(Mg24)*state.rho) +
        (screened_rates(k_Al27_He4_to_p_Si30_reaclib)*Y(Al27)*Y(He4)*state.rho + -screened_rates(k_Si30_p_to_He4_Al27_reaclib)*Y(Si30)*Y(H1)*state.rho) +
        (screened_rates(k_Si27_n_to_p_Al27_reaclib)*Y(Si27)*Y(N)*state.rho + -screened_rates(k_Al27_p_to_n_Si27_reaclib)*Y(Al27)*Y(H1)*state.rho) +
        (-screened_rates(k_Si31_p_to_n_P31_reaclib)*Y(Si31)*Y(H1)*state.rho + screened_rates(k_P31_n_to_p_Si31_reaclib)*Y(P31)*Y(N)*state.rho) +
        (screened_rates(k_P29_n_to_p_Si29_reaclib)*Y(P29)*Y(N)*state.rho + -screened_rates(k_Si29_p_to_n_P29_reaclib)*Y(Si29)*Y(H1)*state.rho) +
        (screened_rates(k_P29_He4_to_p_S32_reaclib)*Y(He4)*Y(P29)*state.rho + -screened_rates(k_S32_p_to_He4_P29_reaclib)*Y(S32)*Y(H1)*state.rho) +
        (-screened_rates(k_P31_p_to_He4_Si28_reaclib)*Y(P31)*Y(H1)*state.rho + screened_rates(k_Si28_He4_to_p_P31_reaclib)*Y(He4)*Y(Si28)*state.rho) +
        (screened_rates(k_P31_He4_to_p_S34_reaclib)*Y(He4)*Y(P31)*state.rho + -screened_rates(k_S34_p_to_He4_P31_reaclib)*Y(S34)*Y(H1)*state.rho) +
        (-screened_rates(k_P33_p_to_He4_Si30_reaclib)*Y(P33)*Y(H1)*state.rho + screened_rates(k_Si30_He4_to_p_P33_reaclib)*Y(He4)*Y(Si30)*state.rho) +
        (screened_rates(k_S33_n_to_p_P33_reaclib)*Y(S33)*Y(N)*state.rho + -screened_rates(k_P33_p_to_n_S33_reaclib)*Y(P33)*Y(H1)*state.rho) +
        (screened_rates(k_Cl33_n_to_p_S33_reaclib)*Y(Cl33)*Y(N)*state.rho + -screened_rates(k_S33_p_to_n_Cl33_reaclib)*Y(S33)*Y(H1)*state.rho) +
        (screened_rates(k_Cl33_He4_to_p_Ar36_reaclib)*Y(Cl33)*Y(He4)*state.rho + -screened_rates(k_Ar36_p_to_He4_Cl33_reaclib)*Y(Ar36)*Y(H1)*state.rho) +
        (screened_rates(k_Cl35_n_to_p_S35_reaclib)*Y(Cl35)*Y(N)*state.rho + -screened_rates(k_S35_p_to_n_Cl35_reaclib)*Y(S35)*Y(H1)*state.rho) +
        (-screened_rates(k_Cl35_p_to_He4_S32_reaclib)*Y(Cl35)*Y(H1)*state.rho + screened_rates(k_S32_He4_to_p_Cl35_reaclib)*Y(He4)*Y(S32)*state.rho) +
        (screened_rates(k_Cl35_He4_to_p_Ar38_reaclib)*Y(Cl35)*Y(He4)*state.rho + -screened_rates(k_Ar38_p_to_He4_Cl35_reaclib)*Y(Ar38)*Y(H1)*state.rho) +
        (-screened_rates(k_Cl37_p_to_He4_S34_reaclib)*Y(Cl37)*Y(H1)*state.rho + screened_rates(k_S34_He4_to_p_Cl37_reaclib)*Y(He4)*Y(S34)*state.rho) +
        (screened_rates(k_Ar37_n_to_p_Cl37_reaclib)*Y(Ar37)*Y(N)*state.rho + -screened_rates(k_Cl37_p_to_n_Ar37_reaclib)*Y(Cl37)*Y(H1)*state.rho) +
        (screened_rates(k_K37_n_to_p_Ar37_reaclib)*Y(K37)*Y(N)*state.rho + -screened_rates(k_Ar37_p_to_n_K37_reaclib)*Y(Ar37)*Y(H1)*state.rho) +
        (screened_rates(k_K37_He4_to_p_Ca40_reaclib)*Y(He4)*Y(K37)*state.rho + -screened_rates(k_Ca40_p_to_He4_K37_reaclib)*Y(Ca40)*Y(H1)*state.rho) +
        (-screened_rates(k_K39_p_to_He4_Ar36_reaclib)*Y(K39)*Y(H1)*state.rho + screened_rates(k_Ar36_He4_to_p_K39_reaclib)*Y(Ar36)*Y(He4)*state.rho) +
        (-screened_rates(k_Ca42_p_to_He4_K39_reaclib)*Y(Ca42)*Y(H1)*state.rho + screened_rates(k_K39_He4_to_p_Ca42_reaclib)*Y(He4)*Y(K39)*state.rho) +
        (screened_rates(k_Sc41_He4_to_p_Ti44_reaclib)*Y(He4)*Y(Sc41)*state.rho + -screened_rates(k_Ti44_p_to_He4_Sc41_reaclib)*Y(Ti44)*Y(H1)*state.rho) +
        (-screened_rates(k_Sc43_p_to_He4_Ca40_reaclib)*Y(Sc43)*Y(H1)*state.rho + screened_rates(k_Ca40_He4_to_p_Sc43_reaclib)*Y(Ca40)*Y(He4)*state.rho) +
        (screened_rates(k_Sc43_He4_to_p_Ti46_reaclib)*Y(He4)*Y(Sc43)*state.rho + -screened_rates(k_Ti46_p_to_He4_Sc43_reaclib)*Y(Ti46)*Y(H1)*state.rho) +
        (-screened_rates(k_Sc45_p_to_He4_Ca42_reaclib)*Y(Sc45)*Y(H1)*state.rho + screened_rates(k_Ca42_He4_to_p_Sc45_reaclib)*Y(Ca42)*Y(He4)*state.rho) +
        (screened_rates(k_Sc45_He4_to_p_Ti48_reaclib)*Y(He4)*Y(Sc45)*state.rho + -screened_rates(k_Ti48_p_to_He4_Sc45_reaclib)*Y(Ti48)*Y(H1)*state.rho) +
        (-screened_rates(k_V47_p_to_He4_Ti44_reaclib)*Y(V47)*Y(H1)*state.rho + screened_rates(k_Ti44_He4_to_p_V47_reaclib)*Y(He4)*Y(Ti44)*state.rho) +
        (screened_rates(k_V47_He4_to_p_Cr50_reaclib)*Y(He4)*Y(V47)*state.rho + -screened_rates(k_Cr50_p_to_He4_V47_reaclib)*Y(Cr50)*Y(H1)*state.rho) +
        (-screened_rates(k_V49_p_to_He4_Ti46_reaclib)*Y(V49)*Y(H1)*state.rho + screened_rates(k_Ti46_He4_to_p_V49_reaclib)*Y(He4)*Y(Ti46)*state.rho) +
        (screened_rates(k_V49_He4_to_p_Cr52_reaclib)*Y(He4)*Y(V49)*state.rho + -screened_rates(k_Cr52_p_to_He4_V49_reaclib)*Y(Cr52)*Y(H1)*state.rho) +
        (-screened_rates(k_V51_p_to_He4_Ti48_reaclib)*Y(V51)*Y(H1)*state.rho + screened_rates(k_Ti48_He4_to_p_V51_reaclib)*Y(He4)*Y(Ti48)*state.rho) +
        (screened_rates(k_Cr48_He4_to_p_Mn51_reaclib)*Y(Cr48)*Y(He4)*state.rho + -screened_rates(k_Mn51_p_to_He4_Cr48_reaclib)*Y(Mn51)*Y(H1)*state.rho) +
        (-screened_rates(k_Cr54_p_to_He4_V51_reaclib)*Y(Cr54)*Y(H1)*state.rho + screened_rates(k_V51_He4_to_p_Cr54_reaclib)*Y(He4)*Y(V51)*state.rho) +
        (screened_rates(k_Mn51_He4_to_p_Fe54_reaclib)*Y(He4)*Y(Mn51)*state.rho + -screened_rates(k_Fe54_p_to_He4_Mn51_reaclib)*Y(Fe54)*Y(H1)*state.rho) +
        (-screened_rates(k_Mn53_p_to_He4_Cr50_reaclib)*Y(Mn53)*Y(H1)*state.rho + screened_rates(k_Cr50_He4_to_p_Mn53_reaclib)*Y(Cr50)*Y(He4)*state.rho) +
        (screened_rates(k_Mn53_He4_to_p_Fe56_reaclib)*Y(He4)*Y(Mn53)*state.rho + -screened_rates(k_Fe56_p_to_He4_Mn53_reaclib)*Y(Fe56)*Y(H1)*state.rho) +
        (-screened_rates(k_Mn55_p_to_He4_Cr52_reaclib)*Y(Mn55)*Y(H1)*state.rho + screened_rates(k_Cr52_He4_to_p_Mn55_reaclib)*Y(Cr52)*Y(He4)*state.rho) +
        (screened_rates(k_Fe52_He4_to_p_Co55_reaclib)*Y(Fe52)*Y(He4)*state.rho + -screened_rates(k_Co55_p_to_He4_Fe52_reaclib)*Y(Co55)*Y(H1)*state.rho) +
        (-screened_rates(k_Fe58_p_to_He4_Mn55_reaclib)*Y(Fe58)*Y(H1)*state.rho + screened_rates(k_Mn55_He4_to_p_Fe58_reaclib)*Y(He4)*Y(Mn55)*state.rho) +
        (screened_rates(k_Co55_He4_to_p_Ni58_reaclib)*Y(Co55)*Y(He4)*state.rho + -screened_rates(k_Ni58_p_to_He4_Co55_reaclib)*Y(Ni58)*Y(H1)*state.rho) +
        (screened_rates(k_Co56_n_to_p_Fe56_reaclib)*Y(Co56)*Y(N)*state.rho + -screened_rates(k_Fe56_p_to_n_Co56_reaclib)*Y(Fe56)*Y(H1)*state.rho) +
        (-screened_rates(k_Co57_p_to_He4_Fe54_reaclib)*Y(Co57)*Y(H1)*state.rho + screened_rates(k_Fe54_He4_to_p_Co57_reaclib)*Y(Fe54)*Y(He4)*state.rho) +
        (screened_rates(k_Co57_He4_to_p_Ni60_reaclib)*Y(Co57)*Y(He4)*state.rho + -screened_rates(k_Ni60_p_to_He4_Co57_reaclib)*Y(Ni60)*Y(H1)*state.rho) +
        (-screened_rates(k_Co59_p_to_He4_Fe56_reaclib)*Y(Co59)*Y(H1)*state.rho + screened_rates(k_Fe56_He4_to_p_Co59_reaclib)*Y(Fe56)*Y(He4)*state.rho) +
        (screened_rates(k_Ni56_n_to_p_Co56_reaclib)*Y(Ni56)*Y(N)*state.rho + -screened_rates(k_Co56_p_to_n_Ni56_reaclib)*Y(Co56)*Y(H1)*state.rho) +
        (-screened_rates(k_Cu59_p_to_He4_Ni56_reaclib)*Y(Cu59)*Y(H1)*state.rho + screened_rates(k_Ni56_He4_to_p_Cu59_reaclib)*Y(He4)*Y(Ni56)*state.rho) +
        (screened_rates(k_n_to_p_beta_neg_weaktab)*Y(N) + -screened_rates(k_p_to_n_electron_capture_weaktab)*Y(H1));

    ydot_nuc(He4) =
        (-screened_rates(k_C12_He4_to_O16_reaclib)*Y(C12)*Y(He4)*state.rho + screened_rates(k_O16_to_He4_C12_reaclib)*Y(O16)) +
        (-screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(He4)*Y(O16)*state.rho + screened_rates(k_Ne20_to_He4_O16_reaclib)*Y(Ne20)) +
        (-screened_rates(k_Mg23_He4_to_Si27_reaclib)*Y(He4)*Y(Mg23)*state.rho + screened_rates(k_Si27_to_He4_Mg23_reaclib)*Y(Si27)) +
        (-screened_rates(k_Mg24_He4_to_Si28_reaclib)*Y(He4)*Y(Mg24)*state.rho + screened_rates(k_Si28_to_He4_Mg24_reaclib)*Y(Si28)) +
        (-screened_rates(k_Al27_He4_to_P31_reaclib)*Y(Al27)*Y(He4)*state.rho + screened_rates(k_P31_to_He4_Al27_reaclib)*Y(P31)) +
        (-screened_rates(k_Si28_He4_to_S32_reaclib)*Y(He4)*Y(Si28)*state.rho + screened_rates(k_S32_to_He4_Si28_reaclib)*Y(S32)) +
        (-screened_rates(k_Si29_He4_to_S33_reaclib)*Y(He4)*Y(Si29)*state.rho + screened_rates(k_S33_to_He4_Si29_reaclib)*Y(S33)) +
        (-screened_rates(k_Si30_He4_to_S34_reaclib)*Y(He4)*Y(Si30)*state.rho + screened_rates(k_S34_to_He4_Si30_reaclib)*Y(S34)) +
        (-screened_rates(k_Si31_He4_to_S35_reaclib)*Y(He4)*Y(Si31)*state.rho + screened_rates(k_S35_to_He4_Si31_reaclib)*Y(S35)) +
        (-screened_rates(k_P29_He4_to_Cl33_reaclib)*Y(He4)*Y(P29)*state.rho + screened_rates(k_Cl33_to_He4_P29_reaclib)*Y(Cl33)) +
        (-screened_rates(k_P31_He4_to_Cl35_reaclib)*Y(He4)*Y(P31)*state.rho + screened_rates(k_Cl35_to_He4_P31_reaclib)*Y(Cl35)) +
        (-screened_rates(k_P33_He4_to_Cl37_reaclib)*Y(He4)*Y(P33)*state.rho + screened_rates(k_Cl37_to_He4_P33_reaclib)*Y(Cl37)) +
        (-screened_rates(k_S32_He4_to_Ar36_reaclib)*Y(He4)*Y(S32)*state.rho + screened_rates(k_Ar36_to_He4_S32_reaclib)*Y(Ar36)) +
        (-screened_rates(k_S33_He4_to_Ar37_reaclib)*Y(He4)*Y(S33)*state.rho + screened_rates(k_Ar37_to_He4_S33_reaclib)*Y(Ar37)) +
        (-screened_rates(k_S34_He4_to_Ar38_reaclib)*Y(He4)*Y(S34)*state.rho + screened_rates(k_Ar38_to_He4_S34_reaclib)*Y(Ar38)) +
        (-screened_rates(k_Cl33_He4_to_K37_reaclib)*Y(Cl33)*Y(He4)*state.rho + screened_rates(k_K37_to_He4_Cl33_reaclib)*Y(K37)) +
        (-screened_rates(k_Cl35_He4_to_K39_reaclib)*Y(Cl35)*Y(He4)*state.rho + screened_rates(k_K39_to_He4_Cl35_reaclib)*Y(K39)) +
        (-screened_rates(k_Ar36_He4_to_Ca40_reaclib)*Y(Ar36)*Y(He4)*state.rho + screened_rates(k_Ca40_to_He4_Ar36_reaclib)*Y(Ca40)) +
        (-screened_rates(k_Ar38_He4_to_Ca42_reaclib)*Y(Ar38)*Y(He4)*state.rho + screened_rates(k_Ca42_to_He4_Ar38_reaclib)*Y(Ca42)) +
        (-screened_rates(k_K37_He4_to_Sc41_reaclib)*Y(He4)*Y(K37)*state.rho + screened_rates(k_Sc41_to_He4_K37_reaclib)*Y(Sc41)) +
        (-screened_rates(k_K39_He4_to_Sc43_reaclib)*Y(He4)*Y(K39)*state.rho + screened_rates(k_Sc43_to_He4_K39_reaclib)*Y(Sc43)) +
        (-screened_rates(k_Ca40_He4_to_Ti44_reaclib)*Y(Ca40)*Y(He4)*state.rho + screened_rates(k_Ti44_to_He4_Ca40_reaclib)*Y(Ti44)) +
        (-screened_rates(k_Ca42_He4_to_Ti46_reaclib)*Y(Ca42)*Y(He4)*state.rho + screened_rates(k_Ti46_to_He4_Ca42_reaclib)*Y(Ti46)) +
        (-screened_rates(k_Sc43_He4_to_V47_reaclib)*Y(He4)*Y(Sc43)*state.rho + screened_rates(k_V47_to_He4_Sc43_reaclib)*Y(V47)) +
        (-screened_rates(k_Sc45_He4_to_V49_reaclib)*Y(He4)*Y(Sc45)*state.rho + screened_rates(k_V49_to_He4_Sc45_reaclib)*Y(V49)) +
        (-screened_rates(k_Ti44_He4_to_Cr48_reaclib)*Y(He4)*Y(Ti44)*state.rho + screened_rates(k_Cr48_to_He4_Ti44_reaclib)*Y(Cr48)) +
        (-screened_rates(k_Ti46_He4_to_Cr50_reaclib)*Y(He4)*Y(Ti46)*state.rho + screened_rates(k_Cr50_to_He4_Ti46_reaclib)*Y(Cr50)) +
        (-screened_rates(k_Ti48_He4_to_Cr52_reaclib)*Y(He4)*Y(Ti48)*state.rho + screened_rates(k_Cr52_to_He4_Ti48_reaclib)*Y(Cr52)) +
        (-screened_rates(k_V47_He4_to_Mn51_reaclib)*Y(He4)*Y(V47)*state.rho + screened_rates(k_Mn51_to_He4_V47_reaclib)*Y(Mn51)) +
        (-screened_rates(k_V49_He4_to_Mn53_reaclib)*Y(He4)*Y(V49)*state.rho + screened_rates(k_Mn53_to_He4_V49_reaclib)*Y(Mn53)) +
        (-screened_rates(k_V51_He4_to_Mn55_reaclib)*Y(He4)*Y(V51)*state.rho + screened_rates(k_Mn55_to_He4_V51_reaclib)*Y(Mn55)) +
        (-screened_rates(k_Cr48_He4_to_Fe52_reaclib)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_Fe52_to_He4_Cr48_reaclib)*Y(Fe52)) +
        (-screened_rates(k_Cr50_He4_to_Fe54_reaclib)*Y(Cr50)*Y(He4)*state.rho + screened_rates(k_Fe54_to_He4_Cr50_reaclib)*Y(Fe54)) +
        (-screened_rates(k_Cr52_He4_to_Fe56_reaclib)*Y(Cr52)*Y(He4)*state.rho + screened_rates(k_Fe56_to_He4_Cr52_reaclib)*Y(Fe56)) +
        (-screened_rates(k_Cr54_He4_to_Fe58_reaclib)*Y(Cr54)*Y(He4)*state.rho + screened_rates(k_Fe58_to_He4_Cr54_reaclib)*Y(Fe58)) +
        (-screened_rates(k_Mn51_He4_to_Co55_reaclib)*Y(He4)*Y(Mn51)*state.rho + screened_rates(k_Co55_to_He4_Mn51_reaclib)*Y(Co55)) +
        (-screened_rates(k_Mn53_He4_to_Co57_reaclib)*Y(He4)*Y(Mn53)*state.rho + screened_rates(k_Co57_to_He4_Mn53_reaclib)*Y(Co57)) +
        (-screened_rates(k_Mn55_He4_to_Co59_reaclib)*Y(He4)*Y(Mn55)*state.rho + screened_rates(k_Co59_to_He4_Mn55_reaclib)*Y(Co59)) +
        (-screened_rates(k_Fe52_He4_to_Ni56_reaclib)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_Ni56_to_He4_Fe52_reaclib)*Y(Ni56)) +
        (-screened_rates(k_Fe54_He4_to_Ni58_reaclib)*Y(Fe54)*Y(He4)*state.rho + screened_rates(k_Ni58_to_He4_Fe54_reaclib)*Y(Ni58)) +
        (-screened_rates(k_Fe56_He4_to_Ni60_reaclib)*Y(Fe56)*Y(He4)*state.rho + screened_rates(k_Ni60_to_He4_Fe56_reaclib)*Y(Ni60)) +
        (-screened_rates(k_Co55_He4_to_Cu59_reaclib)*Y(Co55)*Y(He4)*state.rho + screened_rates(k_Cu59_to_He4_Co55_reaclib)*Y(Cu59)) +
        (-screened_rates(k_Ni56_He4_to_Zn60_reaclib)*Y(He4)*Y(Ni56)*state.rho + screened_rates(k_Zn60_to_He4_Ni56_reaclib)*Y(Zn60)) +
        (screened_rates(k_O16_C12_to_He4_Mg24_reaclib)*Y(C12)*Y(O16)*state.rho + -screened_rates(k_Mg24_He4_to_C12_O16_reaclib)*Y(He4)*Y(Mg24)*state.rho) +
        (0.5*screened_rates(k_O16_O16_to_He4_Si28_reaclib)*amrex::Math::powi<2>(Y(O16))*state.rho + -screened_rates(k_Si28_He4_to_O16_O16_reaclib)*Y(He4)*Y(Si28)*state.rho) +
        (screened_rates(k_Ne20_C12_to_He4_Si28_reaclib)*Y(C12)*Y(Ne20)*state.rho + -screened_rates(k_Si28_He4_to_C12_Ne20_reaclib)*Y(He4)*Y(Si28)*state.rho) +
        (screened_rates(k_Mg23_n_to_He4_Ne20_reaclib)*Y(Mg23)*Y(N)*state.rho + -screened_rates(k_Ne20_He4_to_n_Mg23_reaclib)*Y(He4)*Y(Ne20)*state.rho) +
        (screened_rates(k_Al27_p_to_He4_Mg24_reaclib)*Y(Al27)*Y(H1)*state.rho + -screened_rates(k_Mg24_He4_to_p_Al27_reaclib)*Y(He4)*Y(Mg24)*state.rho) +
        (-screened_rates(k_Al27_He4_to_p_Si30_reaclib)*Y(Al27)*Y(He4)*state.rho + screened_rates(k_Si30_p_to_He4_Al27_reaclib)*Y(Si30)*Y(H1)*state.rho) +
        (screened_rates(k_Si27_n_to_He4_Mg24_reaclib)*Y(Si27)*Y(N)*state.rho + -screened_rates(k_Mg24_He4_to_n_Si27_reaclib)*Y(He4)*Y(Mg24)*state.rho) +
        (-screened_rates(k_Si31_He4_to_n_S34_reaclib)*Y(He4)*Y(Si31)*state.rho + screened_rates(k_S34_n_to_He4_Si31_reaclib)*Y(S34)*Y(N)*state.rho) +
        (-screened_rates(k_P29_He4_to_p_S32_reaclib)*Y(He4)*Y(P29)*state.rho + screened_rates(k_S32_p_to_He4_P29_reaclib)*Y(S32)*Y(H1)*state.rho) +
        (screened_rates(k_P31_p_to_He4_Si28_reaclib)*Y(P31)*Y(H1)*state.rho + -screened_rates(k_Si28_He4_to_p_P31_reaclib)*Y(He4)*Y(Si28)*state.rho) +
        (-screened_rates(k_P31_He4_to_p_S34_reaclib)*Y(He4)*Y(P31)*state.rho + screened_rates(k_S34_p_to_He4_P31_reaclib)*Y(S34)*Y(H1)*state.rho) +
        (screened_rates(k_P33_p_to_He4_Si30_reaclib)*Y(P33)*Y(H1)*state.rho + -screened_rates(k_Si30_He4_to_p_P33_reaclib)*Y(He4)*Y(Si30)*state.rho) +
        (screened_rates(k_S32_n_to_He4_Si29_reaclib)*Y(S32)*Y(N)*state.rho + -screened_rates(k_Si29_He4_to_n_S32_reaclib)*Y(He4)*Y(Si29)*state.rho) +
        (screened_rates(k_S33_n_to_He4_Si30_reaclib)*Y(S33)*Y(N)*state.rho + -screened_rates(k_Si30_He4_to_n_S33_reaclib)*Y(He4)*Y(Si30)*state.rho) +
        (-screened_rates(k_S35_He4_to_n_Ar38_reaclib)*Y(He4)*Y(S35)*state.rho + screened_rates(k_Ar38_n_to_He4_S35_reaclib)*Y(Ar38)*Y(N)*state.rho) +
        (-screened_rates(k_Cl33_He4_to_p_Ar36_reaclib)*Y(Cl33)*Y(He4)*state.rho + screened_rates(k_Ar36_p_to_He4_Cl33_reaclib)*Y(Ar36)*Y(H1)*state.rho) +
        (screened_rates(k_Cl35_p_to_He4_S32_reaclib)*Y(Cl35)*Y(H1)*state.rho + -screened_rates(k_S32_He4_to_p_Cl35_reaclib)*Y(He4)*Y(S32)*state.rho) +
        (-screened_rates(k_Cl35_He4_to_p_Ar38_reaclib)*Y(Cl35)*Y(He4)*state.rho + screened_rates(k_Ar38_p_to_He4_Cl35_reaclib)*Y(Ar38)*Y(H1)*state.rho) +
        (screened_rates(k_Cl37_p_to_He4_S34_reaclib)*Y(Cl37)*Y(H1)*state.rho + -screened_rates(k_S34_He4_to_p_Cl37_reaclib)*Y(He4)*Y(S34)*state.rho) +
        (screened_rates(k_Ar36_n_to_He4_S33_reaclib)*Y(Ar36)*Y(N)*state.rho + -screened_rates(k_S33_He4_to_n_Ar36_reaclib)*Y(He4)*Y(S33)*state.rho) +
        (screened_rates(k_Ar37_n_to_He4_S34_reaclib)*Y(Ar37)*Y(N)*state.rho + -screened_rates(k_S34_He4_to_n_Ar37_reaclib)*Y(He4)*Y(S34)*state.rho) +
        (-screened_rates(k_K37_He4_to_p_Ca40_reaclib)*Y(He4)*Y(K37)*state.rho + screened_rates(k_Ca40_p_to_He4_K37_reaclib)*Y(Ca40)*Y(H1)*state.rho) +
        (screened_rates(k_K39_p_to_He4_Ar36_reaclib)*Y(K39)*Y(H1)*state.rho + -screened_rates(k_Ar36_He4_to_p_K39_reaclib)*Y(Ar36)*Y(He4)*state.rho) +
        (screened_rates(k_Ca40_n_to_He4_Ar37_reaclib)*Y(Ca40)*Y(N)*state.rho + -screened_rates(k_Ar37_He4_to_n_Ca40_reaclib)*Y(Ar37)*Y(He4)*state.rho) +
        (screened_rates(k_Ca42_p_to_He4_K39_reaclib)*Y(Ca42)*Y(H1)*state.rho + -screened_rates(k_K39_He4_to_p_Ca42_reaclib)*Y(He4)*Y(K39)*state.rho) +
        (-screened_rates(k_Sc41_He4_to_p_Ti44_reaclib)*Y(He4)*Y(Sc41)*state.rho + screened_rates(k_Ti44_p_to_He4_Sc41_reaclib)*Y(Ti44)*Y(H1)*state.rho) +
        (screened_rates(k_Sc43_p_to_He4_Ca40_reaclib)*Y(Sc43)*Y(H1)*state.rho + -screened_rates(k_Ca40_He4_to_p_Sc43_reaclib)*Y(Ca40)*Y(He4)*state.rho) +
        (-screened_rates(k_Sc43_He4_to_p_Ti46_reaclib)*Y(He4)*Y(Sc43)*state.rho + screened_rates(k_Ti46_p_to_He4_Sc43_reaclib)*Y(Ti46)*Y(H1)*state.rho) +
        (screened_rates(k_Sc45_p_to_He4_Ca42_reaclib)*Y(Sc45)*Y(H1)*state.rho + -screened_rates(k_Ca42_He4_to_p_Sc45_reaclib)*Y(Ca42)*Y(He4)*state.rho) +
        (-screened_rates(k_Sc45_He4_to_p_Ti48_reaclib)*Y(He4)*Y(Sc45)*state.rho + screened_rates(k_Ti48_p_to_He4_Sc45_reaclib)*Y(Ti48)*Y(H1)*state.rho) +
        (screened_rates(k_V47_p_to_He4_Ti44_reaclib)*Y(V47)*Y(H1)*state.rho + -screened_rates(k_Ti44_He4_to_p_V47_reaclib)*Y(He4)*Y(Ti44)*state.rho) +
        (-screened_rates(k_V47_He4_to_p_Cr50_reaclib)*Y(He4)*Y(V47)*state.rho + screened_rates(k_Cr50_p_to_He4_V47_reaclib)*Y(Cr50)*Y(H1)*state.rho) +
        (screened_rates(k_V49_p_to_He4_Ti46_reaclib)*Y(V49)*Y(H1)*state.rho + -screened_rates(k_Ti46_He4_to_p_V49_reaclib)*Y(He4)*Y(Ti46)*state.rho) +
        (-screened_rates(k_V49_He4_to_p_Cr52_reaclib)*Y(He4)*Y(V49)*state.rho + screened_rates(k_Cr52_p_to_He4_V49_reaclib)*Y(Cr52)*Y(H1)*state.rho) +
        (screened_rates(k_V51_p_to_He4_Ti48_reaclib)*Y(V51)*Y(H1)*state.rho + -screened_rates(k_Ti48_He4_to_p_V51_reaclib)*Y(He4)*Y(Ti48)*state.rho) +
        (-screened_rates(k_Cr48_He4_to_p_Mn51_reaclib)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_Mn51_p_to_He4_Cr48_reaclib)*Y(Mn51)*Y(H1)*state.rho) +
        (screened_rates(k_Cr54_p_to_He4_V51_reaclib)*Y(Cr54)*Y(H1)*state.rho + -screened_rates(k_V51_He4_to_p_Cr54_reaclib)*Y(He4)*Y(V51)*state.rho) +
        (-screened_rates(k_Mn51_He4_to_p_Fe54_reaclib)*Y(He4)*Y(Mn51)*state.rho + screened_rates(k_Fe54_p_to_He4_Mn51_reaclib)*Y(Fe54)*Y(H1)*state.rho) +
        (screened_rates(k_Mn53_p_to_He4_Cr50_reaclib)*Y(Mn53)*Y(H1)*state.rho + -screened_rates(k_Cr50_He4_to_p_Mn53_reaclib)*Y(Cr50)*Y(He4)*state.rho) +
        (-screened_rates(k_Mn53_He4_to_p_Fe56_reaclib)*Y(He4)*Y(Mn53)*state.rho + screened_rates(k_Fe56_p_to_He4_Mn53_reaclib)*Y(Fe56)*Y(H1)*state.rho) +
        (screened_rates(k_Mn55_p_to_He4_Cr52_reaclib)*Y(Mn55)*Y(H1)*state.rho + -screened_rates(k_Cr52_He4_to_p_Mn55_reaclib)*Y(Cr52)*Y(He4)*state.rho) +
        (-screened_rates(k_Fe52_He4_to_p_Co55_reaclib)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_Co55_p_to_He4_Fe52_reaclib)*Y(Co55)*Y(H1)*state.rho) +
        (screened_rates(k_Fe58_p_to_He4_Mn55_reaclib)*Y(Fe58)*Y(H1)*state.rho + -screened_rates(k_Mn55_He4_to_p_Fe58_reaclib)*Y(He4)*Y(Mn55)*state.rho) +
        (-screened_rates(k_Co55_He4_to_p_Ni58_reaclib)*Y(Co55)*Y(He4)*state.rho + screened_rates(k_Ni58_p_to_He4_Co55_reaclib)*Y(Ni58)*Y(H1)*state.rho) +
        (screened_rates(k_Co56_n_to_He4_Mn53_reaclib)*Y(Co56)*Y(N)*state.rho + -screened_rates(k_Mn53_He4_to_n_Co56_reaclib)*Y(He4)*Y(Mn53)*state.rho) +
        (screened_rates(k_Co57_p_to_He4_Fe54_reaclib)*Y(Co57)*Y(H1)*state.rho + -screened_rates(k_Fe54_He4_to_p_Co57_reaclib)*Y(Fe54)*Y(He4)*state.rho) +
        (-screened_rates(k_Co57_He4_to_p_Ni60_reaclib)*Y(Co57)*Y(He4)*state.rho + screened_rates(k_Ni60_p_to_He4_Co57_reaclib)*Y(Ni60)*Y(H1)*state.rho) +
        (screened_rates(k_Co59_p_to_He4_Fe56_reaclib)*Y(Co59)*Y(H1)*state.rho + -screened_rates(k_Fe56_He4_to_p_Co59_reaclib)*Y(Fe56)*Y(He4)*state.rho) +
        (screened_rates(k_Cu59_n_to_He4_Co56_reaclib)*Y(Cu59)*Y(N)*state.rho + -screened_rates(k_Co56_He4_to_n_Cu59_reaclib)*Y(Co56)*Y(He4)*state.rho) +
        (screened_rates(k_Cu59_p_to_He4_Ni56_reaclib)*Y(Cu59)*Y(H1)*state.rho + -screened_rates(k_Ni56_He4_to_p_Cu59_reaclib)*Y(He4)*Y(Ni56)*state.rho) +
        (-0.5*screened_rates(k_He4_He4_He4_to_C12_reaclib)*amrex::Math::powi<3>(Y(He4))*amrex::Math::powi<2>(state.rho) + 3.0*screened_rates(k_C12_to_3He4_reaclib)*Y(C12)) +
        (0.5*screened_rates(k_C12_C12_to_Ne20_He4_approx)*amrex::Math::powi<2>(Y(C12))*state.rho + -screened_rates(k_Ne20_He4_to_C12_C12_approx)*Y(He4)*Y(Ne20)*state.rho) +
        (-screened_rates(k_Ne20_He4_to_Mg24_approx)*Y(He4)*Y(Ne20)*state.rho + screened_rates(k_Mg24_to_Ne20_He4_approx)*Y(Mg24));

    ydot_nuc(C12) =
        (-screened_rates(k_C12_He4_to_O16_reaclib)*Y(C12)*Y(He4)*state.rho + screened_rates(k_O16_to_He4_C12_reaclib)*Y(O16)) +
        (-screened_rates(k_O16_C12_to_p_Al27_reaclib)*Y(C12)*Y(O16)*state.rho + screened_rates(k_Al27_p_to_C12_O16_reaclib)*Y(Al27)*Y(H1)*state.rho) +
        (-screened_rates(k_O16_C12_to_He4_Mg24_reaclib)*Y(C12)*Y(O16)*state.rho + screened_rates(k_Mg24_He4_to_C12_O16_reaclib)*Y(He4)*Y(Mg24)*state.rho) +
        (-screened_rates(k_Ne20_C12_to_p_P31_reaclib)*Y(C12)*Y(Ne20)*state.rho + screened_rates(k_P31_p_to_C12_Ne20_reaclib)*Y(P31)*Y(H1)*state.rho) +
        (-screened_rates(k_Ne20_C12_to_He4_Si28_reaclib)*Y(C12)*Y(Ne20)*state.rho + screened_rates(k_Si28_He4_to_C12_Ne20_reaclib)*Y(He4)*Y(Si28)*state.rho) +
        (2.0*screened_rates(k_Mg23_n_to_C12_C12_reaclib)*Y(Mg23)*Y(N)*state.rho + -screened_rates(k_C12_C12_to_n_Mg23_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho) +
        (screened_rates(k_Si27_n_to_C12_O16_reaclib)*Y(Si27)*Y(N)*state.rho + -screened_rates(k_O16_C12_to_n_Si27_reaclib)*Y(C12)*Y(O16)*state.rho) +
        (0.16666666666666667*screened_rates(k_He4_He4_He4_to_C12_reaclib)*amrex::Math::powi<3>(Y(He4))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_C12_to_3He4_reaclib)*Y(C12)) +
        (-screened_rates(k_C12_C12_to_Mg24_approx)*amrex::Math::powi<2>(Y(C12))*state.rho + 2.0*screened_rates(k_Mg24_to_C12_C12_approx)*Y(Mg24)) +
        (-screened_rates(k_C12_C12_to_Ne20_He4_approx)*amrex::Math::powi<2>(Y(C12))*state.rho + 2.0*screened_rates(k_Ne20_He4_to_C12_C12_approx)*Y(He4)*Y(Ne20)*state.rho);

    ydot_nuc(O16) =
        (screened_rates(k_C12_He4_to_O16_reaclib)*Y(C12)*Y(He4)*state.rho + -screened_rates(k_O16_to_He4_C12_reaclib)*Y(O16)) +
        (-screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(He4)*Y(O16)*state.rho + screened_rates(k_Ne20_to_He4_O16_reaclib)*Y(Ne20)) +
        (-screened_rates(k_O16_C12_to_p_Al27_reaclib)*Y(C12)*Y(O16)*state.rho + screened_rates(k_Al27_p_to_C12_O16_reaclib)*Y(Al27)*Y(H1)*state.rho) +
        (-screened_rates(k_O16_C12_to_He4_Mg24_reaclib)*Y(C12)*Y(O16)*state.rho + screened_rates(k_Mg24_He4_to_C12_O16_reaclib)*Y(He4)*Y(Mg24)*state.rho) +
        (-screened_rates(k_O16_O16_to_p_P31_reaclib)*amrex::Math::powi<2>(Y(O16))*state.rho + 2.0*screened_rates(k_P31_p_to_O16_O16_reaclib)*Y(P31)*Y(H1)*state.rho) +
        (-screened_rates(k_O16_O16_to_He4_Si28_reaclib)*amrex::Math::powi<2>(Y(O16))*state.rho + 2.0*screened_rates(k_Si28_He4_to_O16_O16_reaclib)*Y(He4)*Y(Si28)*state.rho) +
        (screened_rates(k_Si27_n_to_C12_O16_reaclib)*Y(Si27)*Y(N)*state.rho + -screened_rates(k_O16_C12_to_n_Si27_reaclib)*Y(C12)*Y(O16)*state.rho);

    ydot_nuc(Ne20) =
        (screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(He4)*Y(O16)*state.rho + -screened_rates(k_Ne20_to_He4_O16_reaclib)*Y(Ne20)) +
        (-screened_rates(k_Ne20_C12_to_p_P31_reaclib)*Y(C12)*Y(Ne20)*state.rho + screened_rates(k_P31_p_to_C12_Ne20_reaclib)*Y(P31)*Y(H1)*state.rho) +
        (-screened_rates(k_Ne20_C12_to_He4_Si28_reaclib)*Y(C12)*Y(Ne20)*state.rho + screened_rates(k_Si28_He4_to_C12_Ne20_reaclib)*Y(He4)*Y(Si28)*state.rho) +
        (screened_rates(k_Mg23_n_to_He4_Ne20_reaclib)*Y(Mg23)*Y(N)*state.rho + -screened_rates(k_Ne20_He4_to_n_Mg23_reaclib)*Y(He4)*Y(Ne20)*state.rho) +
        (0.5*screened_rates(k_C12_C12_to_Ne20_He4_approx)*amrex::Math::powi<2>(Y(C12))*state.rho + -screened_rates(k_Ne20_He4_to_C12_C12_approx)*Y(He4)*Y(Ne20)*state.rho) +
        (-screened_rates(k_Ne20_He4_to_Mg24_approx)*Y(He4)*Y(Ne20)*state.rho + screened_rates(k_Mg24_to_Ne20_He4_approx)*Y(Mg24));

    ydot_nuc(Mg23) =
        (-screened_rates(k_Mg23_n_to_Mg24_reaclib)*Y(Mg23)*Y(N)*state.rho + screened_rates(k_Mg24_to_n_Mg23_reaclib)*Y(Mg24)) +
        (-screened_rates(k_Mg23_He4_to_Si27_reaclib)*Y(He4)*Y(Mg23)*state.rho + screened_rates(k_Si27_to_He4_Mg23_reaclib)*Y(Si27)) +
        (-screened_rates(k_Mg23_n_to_He4_Ne20_reaclib)*Y(Mg23)*Y(N)*state.rho + screened_rates(k_Ne20_He4_to_n_Mg23_reaclib)*Y(He4)*Y(Ne20)*state.rho) +
        (-screened_rates(k_Mg23_n_to_C12_C12_reaclib)*Y(Mg23)*Y(N)*state.rho + 0.5*screened_rates(k_C12_C12_to_n_Mg23_reaclib)*amrex::Math::powi<2>(Y(C12))*state.rho);

    ydot_nuc(Mg24) =
        (screened_rates(k_Mg23_n_to_Mg24_reaclib)*Y(Mg23)*Y(N)*state.rho + -screened_rates(k_Mg24_to_n_Mg23_reaclib)*Y(Mg24)) +
        (-screened_rates(k_Mg24_He4_to_Si28_reaclib)*Y(He4)*Y(Mg24)*state.rho + screened_rates(k_Si28_to_He4_Mg24_reaclib)*Y(Si28)) +
        (screened_rates(k_O16_C12_to_He4_Mg24_reaclib)*Y(C12)*Y(O16)*state.rho + -screened_rates(k_Mg24_He4_to_C12_O16_reaclib)*Y(He4)*Y(Mg24)*state.rho) +
        (screened_rates(k_Al27_p_to_He4_Mg24_reaclib)*Y(Al27)*Y(H1)*state.rho + -screened_rates(k_Mg24_He4_to_p_Al27_reaclib)*Y(He4)*Y(Mg24)*state.rho) +
        (screened_rates(k_Si27_n_to_He4_Mg24_reaclib)*Y(Si27)*Y(N)*state.rho + -screened_rates(k_Mg24_He4_to_n_Si27_reaclib)*Y(He4)*Y(Mg24)*state.rho) +
        (0.5*screened_rates(k_C12_C12_to_Mg24_approx)*amrex::Math::powi<2>(Y(C12))*state.rho + -screened_rates(k_Mg24_to_C12_C12_approx)*Y(Mg24)) +
        (screened_rates(k_Ne20_He4_to_Mg24_approx)*Y(He4)*Y(Ne20)*state.rho + -screened_rates(k_Mg24_to_Ne20_He4_approx)*Y(Mg24));

    ydot_nuc(Al27) =
        (-screened_rates(k_Al27_p_to_Si28_reaclib)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_Si28_to_p_Al27_reaclib)*Y(Si28)) +
        (-screened_rates(k_Al27_He4_to_P31_reaclib)*Y(Al27)*Y(He4)*state.rho + screened_rates(k_P31_to_He4_Al27_reaclib)*Y(P31)) +
        (screened_rates(k_O16_C12_to_p_Al27_reaclib)*Y(C12)*Y(O16)*state.rho + -screened_rates(k_Al27_p_to_C12_O16_reaclib)*Y(Al27)*Y(H1)*state.rho) +
        (-screened_rates(k_Al27_p_to_He4_Mg24_reaclib)*Y(Al27)*Y(H1)*state.rho + screened_rates(k_Mg24_He4_to_p_Al27_reaclib)*Y(He4)*Y(Mg24)*state.rho) +
        (-screened_rates(k_Al27_He4_to_p_Si30_reaclib)*Y(Al27)*Y(He4)*state.rho + screened_rates(k_Si30_p_to_He4_Al27_reaclib)*Y(Si30)*Y(H1)*state.rho) +
        (screened_rates(k_Si27_n_to_p_Al27_reaclib)*Y(Si27)*Y(N)*state.rho + -screened_rates(k_Al27_p_to_n_Si27_reaclib)*Y(Al27)*Y(H1)*state.rho) +
        (screened_rates(k_Si27_to_Al27_electron_capture_weaktab)*Y(Si27) + -screened_rates(k_Al27_to_Si27_beta_neg_weaktab)*Y(Al27));

    ydot_nuc(Si27) =
        (screened_rates(k_Mg23_He4_to_Si27_reaclib)*Y(He4)*Y(Mg23)*state.rho + -screened_rates(k_Si27_to_He4_Mg23_reaclib)*Y(Si27)) +
        (-screened_rates(k_Si27_n_to_Si28_reaclib)*Y(Si27)*Y(N)*state.rho + screened_rates(k_Si28_to_n_Si27_reaclib)*Y(Si28)) +
        (-screened_rates(k_Si27_n_to_p_Al27_reaclib)*Y(Si27)*Y(N)*state.rho + screened_rates(k_Al27_p_to_n_Si27_reaclib)*Y(Al27)*Y(H1)*state.rho) +
        (-screened_rates(k_Si27_n_to_He4_Mg24_reaclib)*Y(Si27)*Y(N)*state.rho + screened_rates(k_Mg24_He4_to_n_Si27_reaclib)*Y(He4)*Y(Mg24)*state.rho) +
        (-screened_rates(k_Si27_n_to_C12_O16_reaclib)*Y(Si27)*Y(N)*state.rho + screened_rates(k_O16_C12_to_n_Si27_reaclib)*Y(C12)*Y(O16)*state.rho) +
        (-screened_rates(k_Si27_to_Al27_electron_capture_weaktab)*Y(Si27) + screened_rates(k_Al27_to_Si27_beta_neg_weaktab)*Y(Al27));

    ydot_nuc(Si28) =
        (screened_rates(k_Mg24_He4_to_Si28_reaclib)*Y(He4)*Y(Mg24)*state.rho + -screened_rates(k_Si28_to_He4_Mg24_reaclib)*Y(Si28)) +
        (screened_rates(k_Al27_p_to_Si28_reaclib)*Y(Al27)*Y(H1)*state.rho + -screened_rates(k_Si28_to_p_Al27_reaclib)*Y(Si28)) +
        (screened_rates(k_Si27_n_to_Si28_reaclib)*Y(Si27)*Y(N)*state.rho + -screened_rates(k_Si28_to_n_Si27_reaclib)*Y(Si28)) +
        (-screened_rates(k_Si28_n_to_Si29_reaclib)*Y(Si28)*Y(N)*state.rho + screened_rates(k_Si29_to_n_Si28_reaclib)*Y(Si29)) +
        (-screened_rates(k_Si28_p_to_P29_reaclib)*Y(Si28)*Y(H1)*state.rho + screened_rates(k_P29_to_p_Si28_reaclib)*Y(P29)) +
        (-screened_rates(k_Si28_He4_to_S32_reaclib)*Y(He4)*Y(Si28)*state.rho + screened_rates(k_S32_to_He4_Si28_reaclib)*Y(S32)) +
        (0.5*screened_rates(k_O16_O16_to_He4_Si28_reaclib)*amrex::Math::powi<2>(Y(O16))*state.rho + -screened_rates(k_Si28_He4_to_O16_O16_reaclib)*Y(He4)*Y(Si28)*state.rho) +
        (screened_rates(k_Ne20_C12_to_He4_Si28_reaclib)*Y(C12)*Y(Ne20)*state.rho + -screened_rates(k_Si28_He4_to_C12_Ne20_reaclib)*Y(He4)*Y(Si28)*state.rho) +
        (screened_rates(k_P31_p_to_He4_Si28_reaclib)*Y(P31)*Y(H1)*state.rho + -screened_rates(k_Si28_He4_to_p_P31_reaclib)*Y(He4)*Y(Si28)*state.rho);

    ydot_nuc(Si29) =
        (screened_rates(k_Si28_n_to_Si29_reaclib)*Y(Si28)*Y(N)*state.rho + -screened_rates(k_Si29_to_n_Si28_reaclib)*Y(Si29)) +
        (-screened_rates(k_Si29_n_to_Si30_reaclib)*Y(Si29)*Y(N)*state.rho + screened_rates(k_Si30_to_n_Si29_reaclib)*Y(Si30)) +
        (-screened_rates(k_Si29_He4_to_S33_reaclib)*Y(He4)*Y(Si29)*state.rho + screened_rates(k_S33_to_He4_Si29_reaclib)*Y(S33)) +
        (screened_rates(k_P29_n_to_p_Si29_reaclib)*Y(P29)*Y(N)*state.rho + -screened_rates(k_Si29_p_to_n_P29_reaclib)*Y(Si29)*Y(H1)*state.rho) +
        (screened_rates(k_S32_n_to_He4_Si29_reaclib)*Y(S32)*Y(N)*state.rho + -screened_rates(k_Si29_He4_to_n_S32_reaclib)*Y(He4)*Y(Si29)*state.rho) +
        (screened_rates(k_P29_to_Si29_electron_capture_weaktab)*Y(P29) + -screened_rates(k_Si29_to_P29_beta_neg_weaktab)*Y(Si29));

    ydot_nuc(Si30) =
        (screened_rates(k_Si29_n_to_Si30_reaclib)*Y(Si29)*Y(N)*state.rho + -screened_rates(k_Si30_to_n_Si29_reaclib)*Y(Si30)) +
        (-screened_rates(k_Si30_n_to_Si31_reaclib)*Y(Si30)*Y(N)*state.rho + screened_rates(k_Si31_to_n_Si30_reaclib)*Y(Si31)) +
        (-screened_rates(k_Si30_p_to_P31_reaclib)*Y(Si30)*Y(H1)*state.rho + screened_rates(k_P31_to_p_Si30_reaclib)*Y(P31)) +
        (-screened_rates(k_Si30_He4_to_S34_reaclib)*Y(He4)*Y(Si30)*state.rho + screened_rates(k_S34_to_He4_Si30_reaclib)*Y(S34)) +
        (screened_rates(k_Al27_He4_to_p_Si30_reaclib)*Y(Al27)*Y(He4)*state.rho + -screened_rates(k_Si30_p_to_He4_Al27_reaclib)*Y(Si30)*Y(H1)*state.rho) +
        (screened_rates(k_P33_p_to_He4_Si30_reaclib)*Y(P33)*Y(H1)*state.rho + -screened_rates(k_Si30_He4_to_p_P33_reaclib)*Y(He4)*Y(Si30)*state.rho) +
        (screened_rates(k_S33_n_to_He4_Si30_reaclib)*Y(S33)*Y(N)*state.rho + -screened_rates(k_Si30_He4_to_n_S33_reaclib)*Y(He4)*Y(Si30)*state.rho);

    ydot_nuc(Si31) =
        (screened_rates(k_Si30_n_to_Si31_reaclib)*Y(Si30)*Y(N)*state.rho + -screened_rates(k_Si31_to_n_Si30_reaclib)*Y(Si31)) +
        (-screened_rates(k_Si31_He4_to_S35_reaclib)*Y(He4)*Y(Si31)*state.rho + screened_rates(k_S35_to_He4_Si31_reaclib)*Y(S35)) +
        (-screened_rates(k_Si31_p_to_n_P31_reaclib)*Y(Si31)*Y(H1)*state.rho + screened_rates(k_P31_n_to_p_Si31_reaclib)*Y(P31)*Y(N)*state.rho) +
        (-screened_rates(k_Si31_He4_to_n_S34_reaclib)*Y(He4)*Y(Si31)*state.rho + screened_rates(k_S34_n_to_He4_Si31_reaclib)*Y(S34)*Y(N)*state.rho) +
        (-screened_rates(k_Si31_to_P31_beta_neg_weaktab)*Y(Si31) + screened_rates(k_P31_to_Si31_electron_capture_weaktab)*Y(P31));

    ydot_nuc(P29) =
        (screened_rates(k_Si28_p_to_P29_reaclib)*Y(Si28)*Y(H1)*state.rho + -screened_rates(k_P29_to_p_Si28_reaclib)*Y(P29)) +
        (-screened_rates(k_P29_He4_to_Cl33_reaclib)*Y(He4)*Y(P29)*state.rho + screened_rates(k_Cl33_to_He4_P29_reaclib)*Y(Cl33)) +
        (-screened_rates(k_P29_n_to_p_Si29_reaclib)*Y(P29)*Y(N)*state.rho + screened_rates(k_Si29_p_to_n_P29_reaclib)*Y(Si29)*Y(H1)*state.rho) +
        (-screened_rates(k_P29_He4_to_p_S32_reaclib)*Y(He4)*Y(P29)*state.rho + screened_rates(k_S32_p_to_He4_P29_reaclib)*Y(S32)*Y(H1)*state.rho) +
        (-screened_rates(k_P29_to_Si29_electron_capture_weaktab)*Y(P29) + screened_rates(k_Si29_to_P29_beta_neg_weaktab)*Y(Si29));

    ydot_nuc(P31) =
        (screened_rates(k_Al27_He4_to_P31_reaclib)*Y(Al27)*Y(He4)*state.rho + -screened_rates(k_P31_to_He4_Al27_reaclib)*Y(P31)) +
        (screened_rates(k_Si30_p_to_P31_reaclib)*Y(Si30)*Y(H1)*state.rho + -screened_rates(k_P31_to_p_Si30_reaclib)*Y(P31)) +
        (-screened_rates(k_P31_p_to_S32_reaclib)*Y(P31)*Y(H1)*state.rho + screened_rates(k_S32_to_p_P31_reaclib)*Y(S32)) +
        (-screened_rates(k_P31_He4_to_Cl35_reaclib)*Y(He4)*Y(P31)*state.rho + screened_rates(k_Cl35_to_He4_P31_reaclib)*Y(Cl35)) +
        (0.5*screened_rates(k_O16_O16_to_p_P31_reaclib)*amrex::Math::powi<2>(Y(O16))*state.rho + -screened_rates(k_P31_p_to_O16_O16_reaclib)*Y(P31)*Y(H1)*state.rho) +
        (screened_rates(k_Ne20_C12_to_p_P31_reaclib)*Y(C12)*Y(Ne20)*state.rho + -screened_rates(k_P31_p_to_C12_Ne20_reaclib)*Y(P31)*Y(H1)*state.rho) +
        (screened_rates(k_Si31_p_to_n_P31_reaclib)*Y(Si31)*Y(H1)*state.rho + -screened_rates(k_P31_n_to_p_Si31_reaclib)*Y(P31)*Y(N)*state.rho) +
        (-screened_rates(k_P31_p_to_He4_Si28_reaclib)*Y(P31)*Y(H1)*state.rho + screened_rates(k_Si28_He4_to_p_P31_reaclib)*Y(He4)*Y(Si28)*state.rho) +
        (-screened_rates(k_P31_He4_to_p_S34_reaclib)*Y(He4)*Y(P31)*state.rho + screened_rates(k_S34_p_to_He4_P31_reaclib)*Y(S34)*Y(H1)*state.rho) +
        (screened_rates(k_Si31_to_P31_beta_neg_weaktab)*Y(Si31) + -screened_rates(k_P31_to_Si31_electron_capture_weaktab)*Y(P31));

    ydot_nuc(P33) =
        (-screened_rates(k_P33_p_to_S34_reaclib)*Y(P33)*Y(H1)*state.rho + screened_rates(k_S34_to_p_P33_reaclib)*Y(S34)) +
        (-screened_rates(k_P33_He4_to_Cl37_reaclib)*Y(He4)*Y(P33)*state.rho + screened_rates(k_Cl37_to_He4_P33_reaclib)*Y(Cl37)) +
        (-screened_rates(k_P33_p_to_He4_Si30_reaclib)*Y(P33)*Y(H1)*state.rho + screened_rates(k_Si30_He4_to_p_P33_reaclib)*Y(He4)*Y(Si30)*state.rho) +
        (screened_rates(k_S33_n_to_p_P33_reaclib)*Y(S33)*Y(N)*state.rho + -screened_rates(k_P33_p_to_n_S33_reaclib)*Y(P33)*Y(H1)*state.rho) +
        (-screened_rates(k_P33_to_S33_beta_neg_weaktab)*Y(P33) + screened_rates(k_S33_to_P33_electron_capture_weaktab)*Y(S33));

    ydot_nuc(S32) =
        (screened_rates(k_Si28_He4_to_S32_reaclib)*Y(He4)*Y(Si28)*state.rho + -screened_rates(k_S32_to_He4_Si28_reaclib)*Y(S32)) +
        (screened_rates(k_P31_p_to_S32_reaclib)*Y(P31)*Y(H1)*state.rho + -screened_rates(k_S32_to_p_P31_reaclib)*Y(S32)) +
        (-screened_rates(k_S32_n_to_S33_reaclib)*Y(S32)*Y(N)*state.rho + screened_rates(k_S33_to_n_S32_reaclib)*Y(S33)) +
        (-screened_rates(k_S32_p_to_Cl33_reaclib)*Y(S32)*Y(H1)*state.rho + screened_rates(k_Cl33_to_p_S32_reaclib)*Y(Cl33)) +
        (-screened_rates(k_S32_He4_to_Ar36_reaclib)*Y(He4)*Y(S32)*state.rho + screened_rates(k_Ar36_to_He4_S32_reaclib)*Y(Ar36)) +
        (screened_rates(k_P29_He4_to_p_S32_reaclib)*Y(He4)*Y(P29)*state.rho + -screened_rates(k_S32_p_to_He4_P29_reaclib)*Y(S32)*Y(H1)*state.rho) +
        (-screened_rates(k_S32_n_to_He4_Si29_reaclib)*Y(S32)*Y(N)*state.rho + screened_rates(k_Si29_He4_to_n_S32_reaclib)*Y(He4)*Y(Si29)*state.rho) +
        (screened_rates(k_Cl35_p_to_He4_S32_reaclib)*Y(Cl35)*Y(H1)*state.rho + -screened_rates(k_S32_He4_to_p_Cl35_reaclib)*Y(He4)*Y(S32)*state.rho);

    ydot_nuc(S33) =
        (screened_rates(k_Si29_He4_to_S33_reaclib)*Y(He4)*Y(Si29)*state.rho + -screened_rates(k_S33_to_He4_Si29_reaclib)*Y(S33)) +
        (screened_rates(k_S32_n_to_S33_reaclib)*Y(S32)*Y(N)*state.rho + -screened_rates(k_S33_to_n_S32_reaclib)*Y(S33)) +
        (-screened_rates(k_S33_n_to_S34_reaclib)*Y(S33)*Y(N)*state.rho + screened_rates(k_S34_to_n_S33_reaclib)*Y(S34)) +
        (-screened_rates(k_S33_He4_to_Ar37_reaclib)*Y(He4)*Y(S33)*state.rho + screened_rates(k_Ar37_to_He4_S33_reaclib)*Y(Ar37)) +
        (-screened_rates(k_S33_n_to_p_P33_reaclib)*Y(S33)*Y(N)*state.rho + screened_rates(k_P33_p_to_n_S33_reaclib)*Y(P33)*Y(H1)*state.rho) +
        (-screened_rates(k_S33_n_to_He4_Si30_reaclib)*Y(S33)*Y(N)*state.rho + screened_rates(k_Si30_He4_to_n_S33_reaclib)*Y(He4)*Y(Si30)*state.rho) +
        (screened_rates(k_Cl33_n_to_p_S33_reaclib)*Y(Cl33)*Y(N)*state.rho + -screened_rates(k_S33_p_to_n_Cl33_reaclib)*Y(S33)*Y(H1)*state.rho) +
        (screened_rates(k_Ar36_n_to_He4_S33_reaclib)*Y(Ar36)*Y(N)*state.rho + -screened_rates(k_S33_He4_to_n_Ar36_reaclib)*Y(He4)*Y(S33)*state.rho) +
        (screened_rates(k_Cl33_to_S33_electron_capture_weaktab)*Y(Cl33) + -screened_rates(k_S33_to_Cl33_beta_neg_weaktab)*Y(S33)) +
        (screened_rates(k_P33_to_S33_beta_neg_weaktab)*Y(P33) + -screened_rates(k_S33_to_P33_electron_capture_weaktab)*Y(S33));

    ydot_nuc(S34) =
        (screened_rates(k_Si30_He4_to_S34_reaclib)*Y(He4)*Y(Si30)*state.rho + -screened_rates(k_S34_to_He4_Si30_reaclib)*Y(S34)) +
        (screened_rates(k_P33_p_to_S34_reaclib)*Y(P33)*Y(H1)*state.rho + -screened_rates(k_S34_to_p_P33_reaclib)*Y(S34)) +
        (screened_rates(k_S33_n_to_S34_reaclib)*Y(S33)*Y(N)*state.rho + -screened_rates(k_S34_to_n_S33_reaclib)*Y(S34)) +
        (-screened_rates(k_S34_n_to_S35_reaclib)*Y(S34)*Y(N)*state.rho + screened_rates(k_S35_to_n_S34_reaclib)*Y(S35)) +
        (-screened_rates(k_S34_p_to_Cl35_reaclib)*Y(S34)*Y(H1)*state.rho + screened_rates(k_Cl35_to_p_S34_reaclib)*Y(Cl35)) +
        (-screened_rates(k_S34_He4_to_Ar38_reaclib)*Y(He4)*Y(S34)*state.rho + screened_rates(k_Ar38_to_He4_S34_reaclib)*Y(Ar38)) +
        (screened_rates(k_Si31_He4_to_n_S34_reaclib)*Y(He4)*Y(Si31)*state.rho + -screened_rates(k_S34_n_to_He4_Si31_reaclib)*Y(S34)*Y(N)*state.rho) +
        (screened_rates(k_P31_He4_to_p_S34_reaclib)*Y(He4)*Y(P31)*state.rho + -screened_rates(k_S34_p_to_He4_P31_reaclib)*Y(S34)*Y(H1)*state.rho) +
        (screened_rates(k_Cl37_p_to_He4_S34_reaclib)*Y(Cl37)*Y(H1)*state.rho + -screened_rates(k_S34_He4_to_p_Cl37_reaclib)*Y(He4)*Y(S34)*state.rho) +
        (screened_rates(k_Ar37_n_to_He4_S34_reaclib)*Y(Ar37)*Y(N)*state.rho + -screened_rates(k_S34_He4_to_n_Ar37_reaclib)*Y(He4)*Y(S34)*state.rho);

    ydot_nuc(S35) =
        (screened_rates(k_Si31_He4_to_S35_reaclib)*Y(He4)*Y(Si31)*state.rho + -screened_rates(k_S35_to_He4_Si31_reaclib)*Y(S35)) +
        (screened_rates(k_S34_n_to_S35_reaclib)*Y(S34)*Y(N)*state.rho + -screened_rates(k_S35_to_n_S34_reaclib)*Y(S35)) +
        (-screened_rates(k_S35_He4_to_n_Ar38_reaclib)*Y(He4)*Y(S35)*state.rho + screened_rates(k_Ar38_n_to_He4_S35_reaclib)*Y(Ar38)*Y(N)*state.rho) +
        (screened_rates(k_Cl35_n_to_p_S35_reaclib)*Y(Cl35)*Y(N)*state.rho + -screened_rates(k_S35_p_to_n_Cl35_reaclib)*Y(S35)*Y(H1)*state.rho) +
        (-screened_rates(k_S35_to_Cl35_beta_neg_weaktab)*Y(S35) + screened_rates(k_Cl35_to_S35_electron_capture_weaktab)*Y(Cl35));

    ydot_nuc(Cl33) =
        (screened_rates(k_P29_He4_to_Cl33_reaclib)*Y(He4)*Y(P29)*state.rho + -screened_rates(k_Cl33_to_He4_P29_reaclib)*Y(Cl33)) +
        (screened_rates(k_S32_p_to_Cl33_reaclib)*Y(S32)*Y(H1)*state.rho + -screened_rates(k_Cl33_to_p_S32_reaclib)*Y(Cl33)) +
        (-screened_rates(k_Cl33_He4_to_K37_reaclib)*Y(Cl33)*Y(He4)*state.rho + screened_rates(k_K37_to_He4_Cl33_reaclib)*Y(K37)) +
        (-screened_rates(k_Cl33_n_to_p_S33_reaclib)*Y(Cl33)*Y(N)*state.rho + screened_rates(k_S33_p_to_n_Cl33_reaclib)*Y(S33)*Y(H1)*state.rho) +
        (-screened_rates(k_Cl33_He4_to_p_Ar36_reaclib)*Y(Cl33)*Y(He4)*state.rho + screened_rates(k_Ar36_p_to_He4_Cl33_reaclib)*Y(Ar36)*Y(H1)*state.rho) +
        (-screened_rates(k_Cl33_to_S33_electron_capture_weaktab)*Y(Cl33) + screened_rates(k_S33_to_Cl33_beta_neg_weaktab)*Y(S33));

    ydot_nuc(Cl35) =
        (screened_rates(k_P31_He4_to_Cl35_reaclib)*Y(He4)*Y(P31)*state.rho + -screened_rates(k_Cl35_to_He4_P31_reaclib)*Y(Cl35)) +
        (screened_rates(k_S34_p_to_Cl35_reaclib)*Y(S34)*Y(H1)*state.rho + -screened_rates(k_Cl35_to_p_S34_reaclib)*Y(Cl35)) +
        (-screened_rates(k_Cl35_p_to_Ar36_reaclib)*Y(Cl35)*Y(H1)*state.rho + screened_rates(k_Ar36_to_p_Cl35_reaclib)*Y(Ar36)) +
        (-screened_rates(k_Cl35_He4_to_K39_reaclib)*Y(Cl35)*Y(He4)*state.rho + screened_rates(k_K39_to_He4_Cl35_reaclib)*Y(K39)) +
        (-screened_rates(k_Cl35_n_to_p_S35_reaclib)*Y(Cl35)*Y(N)*state.rho + screened_rates(k_S35_p_to_n_Cl35_reaclib)*Y(S35)*Y(H1)*state.rho) +
        (-screened_rates(k_Cl35_p_to_He4_S32_reaclib)*Y(Cl35)*Y(H1)*state.rho + screened_rates(k_S32_He4_to_p_Cl35_reaclib)*Y(He4)*Y(S32)*state.rho) +
        (-screened_rates(k_Cl35_He4_to_p_Ar38_reaclib)*Y(Cl35)*Y(He4)*state.rho + screened_rates(k_Ar38_p_to_He4_Cl35_reaclib)*Y(Ar38)*Y(H1)*state.rho) +
        (screened_rates(k_S35_to_Cl35_beta_neg_weaktab)*Y(S35) + -screened_rates(k_Cl35_to_S35_electron_capture_weaktab)*Y(Cl35));

    ydot_nuc(Cl37) =
        (screened_rates(k_P33_He4_to_Cl37_reaclib)*Y(He4)*Y(P33)*state.rho + -screened_rates(k_Cl37_to_He4_P33_reaclib)*Y(Cl37)) +
        (-screened_rates(k_Cl37_p_to_Ar38_reaclib)*Y(Cl37)*Y(H1)*state.rho + screened_rates(k_Ar38_to_p_Cl37_reaclib)*Y(Ar38)) +
        (-screened_rates(k_Cl37_p_to_He4_S34_reaclib)*Y(Cl37)*Y(H1)*state.rho + screened_rates(k_S34_He4_to_p_Cl37_reaclib)*Y(He4)*Y(S34)*state.rho) +
        (screened_rates(k_Ar37_n_to_p_Cl37_reaclib)*Y(Ar37)*Y(N)*state.rho + -screened_rates(k_Cl37_p_to_n_Ar37_reaclib)*Y(Cl37)*Y(H1)*state.rho) +
        (screened_rates(k_Ar37_to_Cl37_electron_capture_weaktab)*Y(Ar37) + -screened_rates(k_Cl37_to_Ar37_beta_neg_weaktab)*Y(Cl37));

    ydot_nuc(Ar36) =
        (screened_rates(k_S32_He4_to_Ar36_reaclib)*Y(He4)*Y(S32)*state.rho + -screened_rates(k_Ar36_to_He4_S32_reaclib)*Y(Ar36)) +
        (screened_rates(k_Cl35_p_to_Ar36_reaclib)*Y(Cl35)*Y(H1)*state.rho + -screened_rates(k_Ar36_to_p_Cl35_reaclib)*Y(Ar36)) +
        (-screened_rates(k_Ar36_n_to_Ar37_reaclib)*Y(Ar36)*Y(N)*state.rho + screened_rates(k_Ar37_to_n_Ar36_reaclib)*Y(Ar37)) +
        (-screened_rates(k_Ar36_p_to_K37_reaclib)*Y(Ar36)*Y(H1)*state.rho + screened_rates(k_K37_to_p_Ar36_reaclib)*Y(K37)) +
        (-screened_rates(k_Ar36_He4_to_Ca40_reaclib)*Y(Ar36)*Y(He4)*state.rho + screened_rates(k_Ca40_to_He4_Ar36_reaclib)*Y(Ca40)) +
        (screened_rates(k_Cl33_He4_to_p_Ar36_reaclib)*Y(Cl33)*Y(He4)*state.rho + -screened_rates(k_Ar36_p_to_He4_Cl33_reaclib)*Y(Ar36)*Y(H1)*state.rho) +
        (-screened_rates(k_Ar36_n_to_He4_S33_reaclib)*Y(Ar36)*Y(N)*state.rho + screened_rates(k_S33_He4_to_n_Ar36_reaclib)*Y(He4)*Y(S33)*state.rho) +
        (screened_rates(k_K39_p_to_He4_Ar36_reaclib)*Y(K39)*Y(H1)*state.rho + -screened_rates(k_Ar36_He4_to_p_K39_reaclib)*Y(Ar36)*Y(He4)*state.rho);

    ydot_nuc(Ar37) =
        (screened_rates(k_S33_He4_to_Ar37_reaclib)*Y(He4)*Y(S33)*state.rho + -screened_rates(k_Ar37_to_He4_S33_reaclib)*Y(Ar37)) +
        (screened_rates(k_Ar36_n_to_Ar37_reaclib)*Y(Ar36)*Y(N)*state.rho + -screened_rates(k_Ar37_to_n_Ar36_reaclib)*Y(Ar37)) +
        (-screened_rates(k_Ar37_n_to_Ar38_reaclib)*Y(Ar37)*Y(N)*state.rho + screened_rates(k_Ar38_to_n_Ar37_reaclib)*Y(Ar38)) +
        (-screened_rates(k_Ar37_n_to_p_Cl37_reaclib)*Y(Ar37)*Y(N)*state.rho + screened_rates(k_Cl37_p_to_n_Ar37_reaclib)*Y(Cl37)*Y(H1)*state.rho) +
        (-screened_rates(k_Ar37_n_to_He4_S34_reaclib)*Y(Ar37)*Y(N)*state.rho + screened_rates(k_S34_He4_to_n_Ar37_reaclib)*Y(He4)*Y(S34)*state.rho) +
        (screened_rates(k_K37_n_to_p_Ar37_reaclib)*Y(K37)*Y(N)*state.rho + -screened_rates(k_Ar37_p_to_n_K37_reaclib)*Y(Ar37)*Y(H1)*state.rho) +
        (screened_rates(k_Ca40_n_to_He4_Ar37_reaclib)*Y(Ca40)*Y(N)*state.rho + -screened_rates(k_Ar37_He4_to_n_Ca40_reaclib)*Y(Ar37)*Y(He4)*state.rho) +
        (-screened_rates(k_Ar37_to_Cl37_electron_capture_weaktab)*Y(Ar37) + screened_rates(k_Cl37_to_Ar37_beta_neg_weaktab)*Y(Cl37)) +
        (screened_rates(k_K37_to_Ar37_electron_capture_weaktab)*Y(K37) + -screened_rates(k_Ar37_to_K37_beta_neg_weaktab)*Y(Ar37));

    ydot_nuc(Ar38) =
        (screened_rates(k_S34_He4_to_Ar38_reaclib)*Y(He4)*Y(S34)*state.rho + -screened_rates(k_Ar38_to_He4_S34_reaclib)*Y(Ar38)) +
        (screened_rates(k_Cl37_p_to_Ar38_reaclib)*Y(Cl37)*Y(H1)*state.rho + -screened_rates(k_Ar38_to_p_Cl37_reaclib)*Y(Ar38)) +
        (screened_rates(k_Ar37_n_to_Ar38_reaclib)*Y(Ar37)*Y(N)*state.rho + -screened_rates(k_Ar38_to_n_Ar37_reaclib)*Y(Ar38)) +
        (-screened_rates(k_Ar38_p_to_K39_reaclib)*Y(Ar38)*Y(H1)*state.rho + screened_rates(k_K39_to_p_Ar38_reaclib)*Y(K39)) +
        (-screened_rates(k_Ar38_He4_to_Ca42_reaclib)*Y(Ar38)*Y(He4)*state.rho + screened_rates(k_Ca42_to_He4_Ar38_reaclib)*Y(Ca42)) +
        (screened_rates(k_S35_He4_to_n_Ar38_reaclib)*Y(He4)*Y(S35)*state.rho + -screened_rates(k_Ar38_n_to_He4_S35_reaclib)*Y(Ar38)*Y(N)*state.rho) +
        (screened_rates(k_Cl35_He4_to_p_Ar38_reaclib)*Y(Cl35)*Y(He4)*state.rho + -screened_rates(k_Ar38_p_to_He4_Cl35_reaclib)*Y(Ar38)*Y(H1)*state.rho);

    ydot_nuc(K37) =
        (screened_rates(k_Cl33_He4_to_K37_reaclib)*Y(Cl33)*Y(He4)*state.rho + -screened_rates(k_K37_to_He4_Cl33_reaclib)*Y(K37)) +
        (screened_rates(k_Ar36_p_to_K37_reaclib)*Y(Ar36)*Y(H1)*state.rho + -screened_rates(k_K37_to_p_Ar36_reaclib)*Y(K37)) +
        (-screened_rates(k_K37_He4_to_Sc41_reaclib)*Y(He4)*Y(K37)*state.rho + screened_rates(k_Sc41_to_He4_K37_reaclib)*Y(Sc41)) +
        (-screened_rates(k_K37_n_to_p_Ar37_reaclib)*Y(K37)*Y(N)*state.rho + screened_rates(k_Ar37_p_to_n_K37_reaclib)*Y(Ar37)*Y(H1)*state.rho) +
        (-screened_rates(k_K37_He4_to_p_Ca40_reaclib)*Y(He4)*Y(K37)*state.rho + screened_rates(k_Ca40_p_to_He4_K37_reaclib)*Y(Ca40)*Y(H1)*state.rho) +
        (-screened_rates(k_K37_to_Ar37_electron_capture_weaktab)*Y(K37) + screened_rates(k_Ar37_to_K37_beta_neg_weaktab)*Y(Ar37));

    ydot_nuc(K39) =
        (screened_rates(k_Cl35_He4_to_K39_reaclib)*Y(Cl35)*Y(He4)*state.rho + -screened_rates(k_K39_to_He4_Cl35_reaclib)*Y(K39)) +
        (screened_rates(k_Ar38_p_to_K39_reaclib)*Y(Ar38)*Y(H1)*state.rho + -screened_rates(k_K39_to_p_Ar38_reaclib)*Y(K39)) +
        (-screened_rates(k_K39_p_to_Ca40_reaclib)*Y(K39)*Y(H1)*state.rho + screened_rates(k_Ca40_to_p_K39_reaclib)*Y(Ca40)) +
        (-screened_rates(k_K39_He4_to_Sc43_reaclib)*Y(He4)*Y(K39)*state.rho + screened_rates(k_Sc43_to_He4_K39_reaclib)*Y(Sc43)) +
        (-screened_rates(k_K39_p_to_He4_Ar36_reaclib)*Y(K39)*Y(H1)*state.rho + screened_rates(k_Ar36_He4_to_p_K39_reaclib)*Y(Ar36)*Y(He4)*state.rho) +
        (screened_rates(k_Ca42_p_to_He4_K39_reaclib)*Y(Ca42)*Y(H1)*state.rho + -screened_rates(k_K39_He4_to_p_Ca42_reaclib)*Y(He4)*Y(K39)*state.rho);

    ydot_nuc(Ca40) =
        (screened_rates(k_Ar36_He4_to_Ca40_reaclib)*Y(Ar36)*Y(He4)*state.rho + -screened_rates(k_Ca40_to_He4_Ar36_reaclib)*Y(Ca40)) +
        (screened_rates(k_K39_p_to_Ca40_reaclib)*Y(K39)*Y(H1)*state.rho + -screened_rates(k_Ca40_to_p_K39_reaclib)*Y(Ca40)) +
        (-screened_rates(k_Ca40_p_to_Sc41_reaclib)*Y(Ca40)*Y(H1)*state.rho + screened_rates(k_Sc41_to_p_Ca40_reaclib)*Y(Sc41)) +
        (-screened_rates(k_Ca40_He4_to_Ti44_reaclib)*Y(Ca40)*Y(He4)*state.rho + screened_rates(k_Ti44_to_He4_Ca40_reaclib)*Y(Ti44)) +
        (screened_rates(k_K37_He4_to_p_Ca40_reaclib)*Y(He4)*Y(K37)*state.rho + -screened_rates(k_Ca40_p_to_He4_K37_reaclib)*Y(Ca40)*Y(H1)*state.rho) +
        (-screened_rates(k_Ca40_n_to_He4_Ar37_reaclib)*Y(Ca40)*Y(N)*state.rho + screened_rates(k_Ar37_He4_to_n_Ca40_reaclib)*Y(Ar37)*Y(He4)*state.rho) +
        (screened_rates(k_Sc43_p_to_He4_Ca40_reaclib)*Y(Sc43)*Y(H1)*state.rho + -screened_rates(k_Ca40_He4_to_p_Sc43_reaclib)*Y(Ca40)*Y(He4)*state.rho) +
        (-0.5*screened_rates(k_Ca40_n_n_to_Ca42_approx)*Y(Ca40)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + screened_rates(k_Ca42_to_Ca40_n_n_approx)*Y(Ca42));

    ydot_nuc(Ca42) =
        (screened_rates(k_Ar38_He4_to_Ca42_reaclib)*Y(Ar38)*Y(He4)*state.rho + -screened_rates(k_Ca42_to_He4_Ar38_reaclib)*Y(Ca42)) +
        (-screened_rates(k_Ca42_p_to_Sc43_reaclib)*Y(Ca42)*Y(H1)*state.rho + screened_rates(k_Sc43_to_p_Ca42_reaclib)*Y(Sc43)) +
        (-screened_rates(k_Ca42_He4_to_Ti46_reaclib)*Y(Ca42)*Y(He4)*state.rho + screened_rates(k_Ti46_to_He4_Ca42_reaclib)*Y(Ti46)) +
        (-screened_rates(k_Ca42_p_to_He4_K39_reaclib)*Y(Ca42)*Y(H1)*state.rho + screened_rates(k_K39_He4_to_p_Ca42_reaclib)*Y(He4)*Y(K39)*state.rho) +
        (screened_rates(k_Sc45_p_to_He4_Ca42_reaclib)*Y(Sc45)*Y(H1)*state.rho + -screened_rates(k_Ca42_He4_to_p_Sc45_reaclib)*Y(Ca42)*Y(He4)*state.rho) +
        (0.5*screened_rates(k_Ca40_n_n_to_Ca42_approx)*Y(Ca40)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_Ca42_to_Ca40_n_n_approx)*Y(Ca42));

    ydot_nuc(Sc41) =
        (screened_rates(k_K37_He4_to_Sc41_reaclib)*Y(He4)*Y(K37)*state.rho + -screened_rates(k_Sc41_to_He4_K37_reaclib)*Y(Sc41)) +
        (screened_rates(k_Ca40_p_to_Sc41_reaclib)*Y(Ca40)*Y(H1)*state.rho + -screened_rates(k_Sc41_to_p_Ca40_reaclib)*Y(Sc41)) +
        (-screened_rates(k_Sc41_He4_to_p_Ti44_reaclib)*Y(He4)*Y(Sc41)*state.rho + screened_rates(k_Ti44_p_to_He4_Sc41_reaclib)*Y(Ti44)*Y(H1)*state.rho) +
        (-0.5*screened_rates(k_Sc41_n_n_to_Sc43_approx)*Y(Sc41)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + screened_rates(k_Sc43_to_Sc41_n_n_approx)*Y(Sc43));

    ydot_nuc(Sc43) =
        (screened_rates(k_K39_He4_to_Sc43_reaclib)*Y(He4)*Y(K39)*state.rho + -screened_rates(k_Sc43_to_He4_K39_reaclib)*Y(Sc43)) +
        (screened_rates(k_Ca42_p_to_Sc43_reaclib)*Y(Ca42)*Y(H1)*state.rho + -screened_rates(k_Sc43_to_p_Ca42_reaclib)*Y(Sc43)) +
        (-screened_rates(k_Sc43_p_to_Ti44_reaclib)*Y(Sc43)*Y(H1)*state.rho + screened_rates(k_Ti44_to_p_Sc43_reaclib)*Y(Ti44)) +
        (-screened_rates(k_Sc43_He4_to_V47_reaclib)*Y(He4)*Y(Sc43)*state.rho + screened_rates(k_V47_to_He4_Sc43_reaclib)*Y(V47)) +
        (-screened_rates(k_Sc43_p_to_He4_Ca40_reaclib)*Y(Sc43)*Y(H1)*state.rho + screened_rates(k_Ca40_He4_to_p_Sc43_reaclib)*Y(Ca40)*Y(He4)*state.rho) +
        (-screened_rates(k_Sc43_He4_to_p_Ti46_reaclib)*Y(He4)*Y(Sc43)*state.rho + screened_rates(k_Ti46_p_to_He4_Sc43_reaclib)*Y(Ti46)*Y(H1)*state.rho) +
        (0.5*screened_rates(k_Sc41_n_n_to_Sc43_approx)*Y(Sc41)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_Sc43_to_Sc41_n_n_approx)*Y(Sc43)) +
        (-0.5*screened_rates(k_Sc43_n_n_to_Sc45_approx)*Y(Sc43)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + screened_rates(k_Sc45_to_Sc43_n_n_approx)*Y(Sc45));

    ydot_nuc(Sc45) =
        (-screened_rates(k_Sc45_p_to_Ti46_reaclib)*Y(Sc45)*Y(H1)*state.rho + screened_rates(k_Ti46_to_p_Sc45_reaclib)*Y(Ti46)) +
        (-screened_rates(k_Sc45_He4_to_V49_reaclib)*Y(He4)*Y(Sc45)*state.rho + screened_rates(k_V49_to_He4_Sc45_reaclib)*Y(V49)) +
        (-screened_rates(k_Sc45_p_to_He4_Ca42_reaclib)*Y(Sc45)*Y(H1)*state.rho + screened_rates(k_Ca42_He4_to_p_Sc45_reaclib)*Y(Ca42)*Y(He4)*state.rho) +
        (-screened_rates(k_Sc45_He4_to_p_Ti48_reaclib)*Y(He4)*Y(Sc45)*state.rho + screened_rates(k_Ti48_p_to_He4_Sc45_reaclib)*Y(Ti48)*Y(H1)*state.rho) +
        (0.5*screened_rates(k_Sc43_n_n_to_Sc45_approx)*Y(Sc43)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_Sc45_to_Sc43_n_n_approx)*Y(Sc45));

    ydot_nuc(Ti44) =
        (screened_rates(k_Ca40_He4_to_Ti44_reaclib)*Y(Ca40)*Y(He4)*state.rho + -screened_rates(k_Ti44_to_He4_Ca40_reaclib)*Y(Ti44)) +
        (screened_rates(k_Sc43_p_to_Ti44_reaclib)*Y(Sc43)*Y(H1)*state.rho + -screened_rates(k_Ti44_to_p_Sc43_reaclib)*Y(Ti44)) +
        (-screened_rates(k_Ti44_He4_to_Cr48_reaclib)*Y(He4)*Y(Ti44)*state.rho + screened_rates(k_Cr48_to_He4_Ti44_reaclib)*Y(Cr48)) +
        (screened_rates(k_Sc41_He4_to_p_Ti44_reaclib)*Y(He4)*Y(Sc41)*state.rho + -screened_rates(k_Ti44_p_to_He4_Sc41_reaclib)*Y(Ti44)*Y(H1)*state.rho) +
        (screened_rates(k_V47_p_to_He4_Ti44_reaclib)*Y(V47)*Y(H1)*state.rho + -screened_rates(k_Ti44_He4_to_p_V47_reaclib)*Y(He4)*Y(Ti44)*state.rho) +
        (-0.5*screened_rates(k_Ti44_n_n_to_Ti46_approx)*Y(Ti44)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + screened_rates(k_Ti46_to_Ti44_n_n_approx)*Y(Ti46));

    ydot_nuc(Ti46) =
        (screened_rates(k_Ca42_He4_to_Ti46_reaclib)*Y(Ca42)*Y(He4)*state.rho + -screened_rates(k_Ti46_to_He4_Ca42_reaclib)*Y(Ti46)) +
        (screened_rates(k_Sc45_p_to_Ti46_reaclib)*Y(Sc45)*Y(H1)*state.rho + -screened_rates(k_Ti46_to_p_Sc45_reaclib)*Y(Ti46)) +
        (-screened_rates(k_Ti46_p_to_V47_reaclib)*Y(Ti46)*Y(H1)*state.rho + screened_rates(k_V47_to_p_Ti46_reaclib)*Y(V47)) +
        (-screened_rates(k_Ti46_He4_to_Cr50_reaclib)*Y(He4)*Y(Ti46)*state.rho + screened_rates(k_Cr50_to_He4_Ti46_reaclib)*Y(Cr50)) +
        (screened_rates(k_Sc43_He4_to_p_Ti46_reaclib)*Y(He4)*Y(Sc43)*state.rho + -screened_rates(k_Ti46_p_to_He4_Sc43_reaclib)*Y(Ti46)*Y(H1)*state.rho) +
        (screened_rates(k_V49_p_to_He4_Ti46_reaclib)*Y(V49)*Y(H1)*state.rho + -screened_rates(k_Ti46_He4_to_p_V49_reaclib)*Y(He4)*Y(Ti46)*state.rho) +
        (0.5*screened_rates(k_Ti44_n_n_to_Ti46_approx)*Y(Ti44)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_Ti46_to_Ti44_n_n_approx)*Y(Ti46)) +
        (-0.5*screened_rates(k_Ti46_n_n_to_Ti48_approx)*Y(Ti46)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + screened_rates(k_Ti48_to_Ti46_n_n_approx)*Y(Ti48));

    ydot_nuc(Ti48) =
        (-screened_rates(k_Ti48_p_to_V49_reaclib)*Y(Ti48)*Y(H1)*state.rho + screened_rates(k_V49_to_p_Ti48_reaclib)*Y(V49)) +
        (-screened_rates(k_Ti48_He4_to_Cr52_reaclib)*Y(He4)*Y(Ti48)*state.rho + screened_rates(k_Cr52_to_He4_Ti48_reaclib)*Y(Cr52)) +
        (screened_rates(k_Sc45_He4_to_p_Ti48_reaclib)*Y(He4)*Y(Sc45)*state.rho + -screened_rates(k_Ti48_p_to_He4_Sc45_reaclib)*Y(Ti48)*Y(H1)*state.rho) +
        (screened_rates(k_V51_p_to_He4_Ti48_reaclib)*Y(V51)*Y(H1)*state.rho + -screened_rates(k_Ti48_He4_to_p_V51_reaclib)*Y(He4)*Y(Ti48)*state.rho) +
        (0.5*screened_rates(k_Ti46_n_n_to_Ti48_approx)*Y(Ti46)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_Ti48_to_Ti46_n_n_approx)*Y(Ti48));

    ydot_nuc(V47) =
        (screened_rates(k_Sc43_He4_to_V47_reaclib)*Y(He4)*Y(Sc43)*state.rho + -screened_rates(k_V47_to_He4_Sc43_reaclib)*Y(V47)) +
        (screened_rates(k_Ti46_p_to_V47_reaclib)*Y(Ti46)*Y(H1)*state.rho + -screened_rates(k_V47_to_p_Ti46_reaclib)*Y(V47)) +
        (-screened_rates(k_V47_p_to_Cr48_reaclib)*Y(V47)*Y(H1)*state.rho + screened_rates(k_Cr48_to_p_V47_reaclib)*Y(Cr48)) +
        (-screened_rates(k_V47_He4_to_Mn51_reaclib)*Y(He4)*Y(V47)*state.rho + screened_rates(k_Mn51_to_He4_V47_reaclib)*Y(Mn51)) +
        (-screened_rates(k_V47_p_to_He4_Ti44_reaclib)*Y(V47)*Y(H1)*state.rho + screened_rates(k_Ti44_He4_to_p_V47_reaclib)*Y(He4)*Y(Ti44)*state.rho) +
        (-screened_rates(k_V47_He4_to_p_Cr50_reaclib)*Y(He4)*Y(V47)*state.rho + screened_rates(k_Cr50_p_to_He4_V47_reaclib)*Y(Cr50)*Y(H1)*state.rho) +
        (-0.5*screened_rates(k_V47_n_n_to_V49_approx)*Y(V47)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + screened_rates(k_V49_to_V47_n_n_approx)*Y(V49));

    ydot_nuc(V49) =
        (screened_rates(k_Sc45_He4_to_V49_reaclib)*Y(He4)*Y(Sc45)*state.rho + -screened_rates(k_V49_to_He4_Sc45_reaclib)*Y(V49)) +
        (screened_rates(k_Ti48_p_to_V49_reaclib)*Y(Ti48)*Y(H1)*state.rho + -screened_rates(k_V49_to_p_Ti48_reaclib)*Y(V49)) +
        (-screened_rates(k_V49_p_to_Cr50_reaclib)*Y(V49)*Y(H1)*state.rho + screened_rates(k_Cr50_to_p_V49_reaclib)*Y(Cr50)) +
        (-screened_rates(k_V49_He4_to_Mn53_reaclib)*Y(He4)*Y(V49)*state.rho + screened_rates(k_Mn53_to_He4_V49_reaclib)*Y(Mn53)) +
        (-screened_rates(k_V49_p_to_He4_Ti46_reaclib)*Y(V49)*Y(H1)*state.rho + screened_rates(k_Ti46_He4_to_p_V49_reaclib)*Y(He4)*Y(Ti46)*state.rho) +
        (-screened_rates(k_V49_He4_to_p_Cr52_reaclib)*Y(He4)*Y(V49)*state.rho + screened_rates(k_Cr52_p_to_He4_V49_reaclib)*Y(Cr52)*Y(H1)*state.rho) +
        (0.5*screened_rates(k_V47_n_n_to_V49_approx)*Y(V47)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_V49_to_V47_n_n_approx)*Y(V49)) +
        (-0.5*screened_rates(k_V49_n_n_to_V51_approx)*Y(V49)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + screened_rates(k_V51_to_V49_n_n_approx)*Y(V51));

    ydot_nuc(V51) =
        (-screened_rates(k_V51_p_to_Cr52_reaclib)*Y(V51)*Y(H1)*state.rho + screened_rates(k_Cr52_to_p_V51_reaclib)*Y(Cr52)) +
        (-screened_rates(k_V51_He4_to_Mn55_reaclib)*Y(He4)*Y(V51)*state.rho + screened_rates(k_Mn55_to_He4_V51_reaclib)*Y(Mn55)) +
        (-screened_rates(k_V51_p_to_He4_Ti48_reaclib)*Y(V51)*Y(H1)*state.rho + screened_rates(k_Ti48_He4_to_p_V51_reaclib)*Y(He4)*Y(Ti48)*state.rho) +
        (screened_rates(k_Cr54_p_to_He4_V51_reaclib)*Y(Cr54)*Y(H1)*state.rho + -screened_rates(k_V51_He4_to_p_Cr54_reaclib)*Y(He4)*Y(V51)*state.rho) +
        (0.5*screened_rates(k_V49_n_n_to_V51_approx)*Y(V49)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_V51_to_V49_n_n_approx)*Y(V51));

    ydot_nuc(Cr48) =
        (screened_rates(k_Ti44_He4_to_Cr48_reaclib)*Y(He4)*Y(Ti44)*state.rho + -screened_rates(k_Cr48_to_He4_Ti44_reaclib)*Y(Cr48)) +
        (screened_rates(k_V47_p_to_Cr48_reaclib)*Y(V47)*Y(H1)*state.rho + -screened_rates(k_Cr48_to_p_V47_reaclib)*Y(Cr48)) +
        (-screened_rates(k_Cr48_He4_to_Fe52_reaclib)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_Fe52_to_He4_Cr48_reaclib)*Y(Fe52)) +
        (-screened_rates(k_Cr48_He4_to_p_Mn51_reaclib)*Y(Cr48)*Y(He4)*state.rho + screened_rates(k_Mn51_p_to_He4_Cr48_reaclib)*Y(Mn51)*Y(H1)*state.rho) +
        (-0.5*screened_rates(k_Cr48_n_n_to_Cr50_approx)*Y(Cr48)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + screened_rates(k_Cr50_to_Cr48_n_n_approx)*Y(Cr50));

    ydot_nuc(Cr50) =
        (screened_rates(k_Ti46_He4_to_Cr50_reaclib)*Y(He4)*Y(Ti46)*state.rho + -screened_rates(k_Cr50_to_He4_Ti46_reaclib)*Y(Cr50)) +
        (screened_rates(k_V49_p_to_Cr50_reaclib)*Y(V49)*Y(H1)*state.rho + -screened_rates(k_Cr50_to_p_V49_reaclib)*Y(Cr50)) +
        (-screened_rates(k_Cr50_p_to_Mn51_reaclib)*Y(Cr50)*Y(H1)*state.rho + screened_rates(k_Mn51_to_p_Cr50_reaclib)*Y(Mn51)) +
        (-screened_rates(k_Cr50_He4_to_Fe54_reaclib)*Y(Cr50)*Y(He4)*state.rho + screened_rates(k_Fe54_to_He4_Cr50_reaclib)*Y(Fe54)) +
        (screened_rates(k_V47_He4_to_p_Cr50_reaclib)*Y(He4)*Y(V47)*state.rho + -screened_rates(k_Cr50_p_to_He4_V47_reaclib)*Y(Cr50)*Y(H1)*state.rho) +
        (screened_rates(k_Mn53_p_to_He4_Cr50_reaclib)*Y(Mn53)*Y(H1)*state.rho + -screened_rates(k_Cr50_He4_to_p_Mn53_reaclib)*Y(Cr50)*Y(He4)*state.rho) +
        (0.5*screened_rates(k_Cr48_n_n_to_Cr50_approx)*Y(Cr48)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_Cr50_to_Cr48_n_n_approx)*Y(Cr50)) +
        (-0.5*screened_rates(k_Cr50_n_n_to_Cr52_approx)*Y(Cr50)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + screened_rates(k_Cr52_to_Cr50_n_n_approx)*Y(Cr52));

    ydot_nuc(Cr52) =
        (screened_rates(k_Ti48_He4_to_Cr52_reaclib)*Y(He4)*Y(Ti48)*state.rho + -screened_rates(k_Cr52_to_He4_Ti48_reaclib)*Y(Cr52)) +
        (screened_rates(k_V51_p_to_Cr52_reaclib)*Y(V51)*Y(H1)*state.rho + -screened_rates(k_Cr52_to_p_V51_reaclib)*Y(Cr52)) +
        (-screened_rates(k_Cr52_p_to_Mn53_reaclib)*Y(Cr52)*Y(H1)*state.rho + screened_rates(k_Mn53_to_p_Cr52_reaclib)*Y(Mn53)) +
        (-screened_rates(k_Cr52_He4_to_Fe56_reaclib)*Y(Cr52)*Y(He4)*state.rho + screened_rates(k_Fe56_to_He4_Cr52_reaclib)*Y(Fe56)) +
        (screened_rates(k_V49_He4_to_p_Cr52_reaclib)*Y(He4)*Y(V49)*state.rho + -screened_rates(k_Cr52_p_to_He4_V49_reaclib)*Y(Cr52)*Y(H1)*state.rho) +
        (screened_rates(k_Mn55_p_to_He4_Cr52_reaclib)*Y(Mn55)*Y(H1)*state.rho + -screened_rates(k_Cr52_He4_to_p_Mn55_reaclib)*Y(Cr52)*Y(He4)*state.rho) +
        (0.5*screened_rates(k_Cr50_n_n_to_Cr52_approx)*Y(Cr50)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_Cr52_to_Cr50_n_n_approx)*Y(Cr52)) +
        (-0.5*screened_rates(k_Cr52_n_n_to_Cr54_approx)*Y(Cr52)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + screened_rates(k_Cr54_to_Cr52_n_n_approx)*Y(Cr54));

    ydot_nuc(Cr54) =
        (-screened_rates(k_Cr54_p_to_Mn55_reaclib)*Y(Cr54)*Y(H1)*state.rho + screened_rates(k_Mn55_to_p_Cr54_reaclib)*Y(Mn55)) +
        (-screened_rates(k_Cr54_He4_to_Fe58_reaclib)*Y(Cr54)*Y(He4)*state.rho + screened_rates(k_Fe58_to_He4_Cr54_reaclib)*Y(Fe58)) +
        (-screened_rates(k_Cr54_p_to_He4_V51_reaclib)*Y(Cr54)*Y(H1)*state.rho + screened_rates(k_V51_He4_to_p_Cr54_reaclib)*Y(He4)*Y(V51)*state.rho) +
        (0.5*screened_rates(k_Cr52_n_n_to_Cr54_approx)*Y(Cr52)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_Cr54_to_Cr52_n_n_approx)*Y(Cr54));

    ydot_nuc(Mn51) =
        (screened_rates(k_V47_He4_to_Mn51_reaclib)*Y(He4)*Y(V47)*state.rho + -screened_rates(k_Mn51_to_He4_V47_reaclib)*Y(Mn51)) +
        (screened_rates(k_Cr50_p_to_Mn51_reaclib)*Y(Cr50)*Y(H1)*state.rho + -screened_rates(k_Mn51_to_p_Cr50_reaclib)*Y(Mn51)) +
        (-screened_rates(k_Mn51_p_to_Fe52_reaclib)*Y(Mn51)*Y(H1)*state.rho + screened_rates(k_Fe52_to_p_Mn51_reaclib)*Y(Fe52)) +
        (-screened_rates(k_Mn51_He4_to_Co55_reaclib)*Y(He4)*Y(Mn51)*state.rho + screened_rates(k_Co55_to_He4_Mn51_reaclib)*Y(Co55)) +
        (screened_rates(k_Cr48_He4_to_p_Mn51_reaclib)*Y(Cr48)*Y(He4)*state.rho + -screened_rates(k_Mn51_p_to_He4_Cr48_reaclib)*Y(Mn51)*Y(H1)*state.rho) +
        (-screened_rates(k_Mn51_He4_to_p_Fe54_reaclib)*Y(He4)*Y(Mn51)*state.rho + screened_rates(k_Fe54_p_to_He4_Mn51_reaclib)*Y(Fe54)*Y(H1)*state.rho) +
        (-0.5*screened_rates(k_Mn51_n_n_to_Mn53_approx)*Y(Mn51)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + screened_rates(k_Mn53_to_Mn51_n_n_approx)*Y(Mn53));

    ydot_nuc(Mn53) =
        (screened_rates(k_V49_He4_to_Mn53_reaclib)*Y(He4)*Y(V49)*state.rho + -screened_rates(k_Mn53_to_He4_V49_reaclib)*Y(Mn53)) +
        (screened_rates(k_Cr52_p_to_Mn53_reaclib)*Y(Cr52)*Y(H1)*state.rho + -screened_rates(k_Mn53_to_p_Cr52_reaclib)*Y(Mn53)) +
        (-screened_rates(k_Mn53_p_to_Fe54_reaclib)*Y(Mn53)*Y(H1)*state.rho + screened_rates(k_Fe54_to_p_Mn53_reaclib)*Y(Fe54)) +
        (-screened_rates(k_Mn53_He4_to_Co57_reaclib)*Y(He4)*Y(Mn53)*state.rho + screened_rates(k_Co57_to_He4_Mn53_reaclib)*Y(Co57)) +
        (-screened_rates(k_Mn53_p_to_He4_Cr50_reaclib)*Y(Mn53)*Y(H1)*state.rho + screened_rates(k_Cr50_He4_to_p_Mn53_reaclib)*Y(Cr50)*Y(He4)*state.rho) +
        (-screened_rates(k_Mn53_He4_to_p_Fe56_reaclib)*Y(He4)*Y(Mn53)*state.rho + screened_rates(k_Fe56_p_to_He4_Mn53_reaclib)*Y(Fe56)*Y(H1)*state.rho) +
        (screened_rates(k_Co56_n_to_He4_Mn53_reaclib)*Y(Co56)*Y(N)*state.rho + -screened_rates(k_Mn53_He4_to_n_Co56_reaclib)*Y(He4)*Y(Mn53)*state.rho) +
        (0.5*screened_rates(k_Mn51_n_n_to_Mn53_approx)*Y(Mn51)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_Mn53_to_Mn51_n_n_approx)*Y(Mn53)) +
        (-0.5*screened_rates(k_Mn53_n_n_to_Mn55_approx)*Y(Mn53)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + screened_rates(k_Mn55_to_Mn53_n_n_approx)*Y(Mn55));

    ydot_nuc(Mn55) =
        (screened_rates(k_V51_He4_to_Mn55_reaclib)*Y(He4)*Y(V51)*state.rho + -screened_rates(k_Mn55_to_He4_V51_reaclib)*Y(Mn55)) +
        (screened_rates(k_Cr54_p_to_Mn55_reaclib)*Y(Cr54)*Y(H1)*state.rho + -screened_rates(k_Mn55_to_p_Cr54_reaclib)*Y(Mn55)) +
        (-screened_rates(k_Mn55_p_to_Fe56_reaclib)*Y(Mn55)*Y(H1)*state.rho + screened_rates(k_Fe56_to_p_Mn55_reaclib)*Y(Fe56)) +
        (-screened_rates(k_Mn55_He4_to_Co59_reaclib)*Y(He4)*Y(Mn55)*state.rho + screened_rates(k_Co59_to_He4_Mn55_reaclib)*Y(Co59)) +
        (-screened_rates(k_Mn55_p_to_He4_Cr52_reaclib)*Y(Mn55)*Y(H1)*state.rho + screened_rates(k_Cr52_He4_to_p_Mn55_reaclib)*Y(Cr52)*Y(He4)*state.rho) +
        (screened_rates(k_Fe58_p_to_He4_Mn55_reaclib)*Y(Fe58)*Y(H1)*state.rho + -screened_rates(k_Mn55_He4_to_p_Fe58_reaclib)*Y(He4)*Y(Mn55)*state.rho) +
        (0.5*screened_rates(k_Mn53_n_n_to_Mn55_approx)*Y(Mn53)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_Mn55_to_Mn53_n_n_approx)*Y(Mn55));

    ydot_nuc(Fe52) =
        (screened_rates(k_Cr48_He4_to_Fe52_reaclib)*Y(Cr48)*Y(He4)*state.rho + -screened_rates(k_Fe52_to_He4_Cr48_reaclib)*Y(Fe52)) +
        (screened_rates(k_Mn51_p_to_Fe52_reaclib)*Y(Mn51)*Y(H1)*state.rho + -screened_rates(k_Fe52_to_p_Mn51_reaclib)*Y(Fe52)) +
        (-screened_rates(k_Fe52_He4_to_Ni56_reaclib)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_Ni56_to_He4_Fe52_reaclib)*Y(Ni56)) +
        (-screened_rates(k_Fe52_He4_to_p_Co55_reaclib)*Y(Fe52)*Y(He4)*state.rho + screened_rates(k_Co55_p_to_He4_Fe52_reaclib)*Y(Co55)*Y(H1)*state.rho) +
        (-0.5*screened_rates(k_Fe52_n_n_to_Fe54_approx)*Y(Fe52)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + screened_rates(k_Fe54_to_Fe52_n_n_approx)*Y(Fe54));

    ydot_nuc(Fe54) =
        (screened_rates(k_Cr50_He4_to_Fe54_reaclib)*Y(Cr50)*Y(He4)*state.rho + -screened_rates(k_Fe54_to_He4_Cr50_reaclib)*Y(Fe54)) +
        (screened_rates(k_Mn53_p_to_Fe54_reaclib)*Y(Mn53)*Y(H1)*state.rho + -screened_rates(k_Fe54_to_p_Mn53_reaclib)*Y(Fe54)) +
        (-screened_rates(k_Fe54_p_to_Co55_reaclib)*Y(Fe54)*Y(H1)*state.rho + screened_rates(k_Co55_to_p_Fe54_reaclib)*Y(Co55)) +
        (-screened_rates(k_Fe54_He4_to_Ni58_reaclib)*Y(Fe54)*Y(He4)*state.rho + screened_rates(k_Ni58_to_He4_Fe54_reaclib)*Y(Ni58)) +
        (screened_rates(k_Mn51_He4_to_p_Fe54_reaclib)*Y(He4)*Y(Mn51)*state.rho + -screened_rates(k_Fe54_p_to_He4_Mn51_reaclib)*Y(Fe54)*Y(H1)*state.rho) +
        (screened_rates(k_Co57_p_to_He4_Fe54_reaclib)*Y(Co57)*Y(H1)*state.rho + -screened_rates(k_Fe54_He4_to_p_Co57_reaclib)*Y(Fe54)*Y(He4)*state.rho) +
        (0.5*screened_rates(k_Fe52_n_n_to_Fe54_approx)*Y(Fe52)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_Fe54_to_Fe52_n_n_approx)*Y(Fe54)) +
        (-0.5*screened_rates(k_Fe54_n_n_to_Fe56_approx)*Y(Fe54)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + screened_rates(k_Fe56_to_Fe54_n_n_approx)*Y(Fe56));

    ydot_nuc(Fe56) =
        (screened_rates(k_Cr52_He4_to_Fe56_reaclib)*Y(Cr52)*Y(He4)*state.rho + -screened_rates(k_Fe56_to_He4_Cr52_reaclib)*Y(Fe56)) +
        (screened_rates(k_Mn55_p_to_Fe56_reaclib)*Y(Mn55)*Y(H1)*state.rho + -screened_rates(k_Fe56_to_p_Mn55_reaclib)*Y(Fe56)) +
        (-screened_rates(k_Fe56_p_to_Co57_reaclib)*Y(Fe56)*Y(H1)*state.rho + screened_rates(k_Co57_to_p_Fe56_reaclib)*Y(Co57)) +
        (-screened_rates(k_Fe56_He4_to_Ni60_reaclib)*Y(Fe56)*Y(He4)*state.rho + screened_rates(k_Ni60_to_He4_Fe56_reaclib)*Y(Ni60)) +
        (screened_rates(k_Mn53_He4_to_p_Fe56_reaclib)*Y(He4)*Y(Mn53)*state.rho + -screened_rates(k_Fe56_p_to_He4_Mn53_reaclib)*Y(Fe56)*Y(H1)*state.rho) +
        (screened_rates(k_Co56_n_to_p_Fe56_reaclib)*Y(Co56)*Y(N)*state.rho + -screened_rates(k_Fe56_p_to_n_Co56_reaclib)*Y(Fe56)*Y(H1)*state.rho) +
        (screened_rates(k_Co59_p_to_He4_Fe56_reaclib)*Y(Co59)*Y(H1)*state.rho + -screened_rates(k_Fe56_He4_to_p_Co59_reaclib)*Y(Fe56)*Y(He4)*state.rho) +
        (0.5*screened_rates(k_Fe54_n_n_to_Fe56_approx)*Y(Fe54)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_Fe56_to_Fe54_n_n_approx)*Y(Fe56)) +
        (-0.5*screened_rates(k_Fe56_n_n_to_Fe58_approx)*Y(Fe56)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + screened_rates(k_Fe58_to_Fe56_n_n_approx)*Y(Fe58)) +
        (screened_rates(k_Co56_to_Fe56_electron_capture_weaktab)*Y(Co56) + -screened_rates(k_Fe56_to_Co56_beta_neg_weaktab)*Y(Fe56));

    ydot_nuc(Fe58) =
        (screened_rates(k_Cr54_He4_to_Fe58_reaclib)*Y(Cr54)*Y(He4)*state.rho + -screened_rates(k_Fe58_to_He4_Cr54_reaclib)*Y(Fe58)) +
        (-screened_rates(k_Fe58_p_to_Co59_reaclib)*Y(Fe58)*Y(H1)*state.rho + screened_rates(k_Co59_to_p_Fe58_reaclib)*Y(Co59)) +
        (-screened_rates(k_Fe58_p_to_He4_Mn55_reaclib)*Y(Fe58)*Y(H1)*state.rho + screened_rates(k_Mn55_He4_to_p_Fe58_reaclib)*Y(He4)*Y(Mn55)*state.rho) +
        (0.5*screened_rates(k_Fe56_n_n_to_Fe58_approx)*Y(Fe56)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_Fe58_to_Fe56_n_n_approx)*Y(Fe58));

    ydot_nuc(Co55) =
        (screened_rates(k_Mn51_He4_to_Co55_reaclib)*Y(He4)*Y(Mn51)*state.rho + -screened_rates(k_Co55_to_He4_Mn51_reaclib)*Y(Co55)) +
        (screened_rates(k_Fe54_p_to_Co55_reaclib)*Y(Fe54)*Y(H1)*state.rho + -screened_rates(k_Co55_to_p_Fe54_reaclib)*Y(Co55)) +
        (-screened_rates(k_Co55_n_to_Co56_reaclib)*Y(Co55)*Y(N)*state.rho + screened_rates(k_Co56_to_n_Co55_reaclib)*Y(Co56)) +
        (-screened_rates(k_Co55_p_to_Ni56_reaclib)*Y(Co55)*Y(H1)*state.rho + screened_rates(k_Ni56_to_p_Co55_reaclib)*Y(Ni56)) +
        (-screened_rates(k_Co55_He4_to_Cu59_reaclib)*Y(Co55)*Y(He4)*state.rho + screened_rates(k_Cu59_to_He4_Co55_reaclib)*Y(Cu59)) +
        (screened_rates(k_Fe52_He4_to_p_Co55_reaclib)*Y(Fe52)*Y(He4)*state.rho + -screened_rates(k_Co55_p_to_He4_Fe52_reaclib)*Y(Co55)*Y(H1)*state.rho) +
        (-screened_rates(k_Co55_He4_to_p_Ni58_reaclib)*Y(Co55)*Y(He4)*state.rho + screened_rates(k_Ni58_p_to_He4_Co55_reaclib)*Y(Ni58)*Y(H1)*state.rho);

    ydot_nuc(Co56) =
        (screened_rates(k_Co55_n_to_Co56_reaclib)*Y(Co55)*Y(N)*state.rho + -screened_rates(k_Co56_to_n_Co55_reaclib)*Y(Co56)) +
        (-screened_rates(k_Co56_n_to_Co57_reaclib)*Y(Co56)*Y(N)*state.rho + screened_rates(k_Co57_to_n_Co56_reaclib)*Y(Co57)) +
        (-screened_rates(k_Co56_n_to_p_Fe56_reaclib)*Y(Co56)*Y(N)*state.rho + screened_rates(k_Fe56_p_to_n_Co56_reaclib)*Y(Fe56)*Y(H1)*state.rho) +
        (-screened_rates(k_Co56_n_to_He4_Mn53_reaclib)*Y(Co56)*Y(N)*state.rho + screened_rates(k_Mn53_He4_to_n_Co56_reaclib)*Y(He4)*Y(Mn53)*state.rho) +
        (screened_rates(k_Ni56_n_to_p_Co56_reaclib)*Y(Ni56)*Y(N)*state.rho + -screened_rates(k_Co56_p_to_n_Ni56_reaclib)*Y(Co56)*Y(H1)*state.rho) +
        (screened_rates(k_Cu59_n_to_He4_Co56_reaclib)*Y(Cu59)*Y(N)*state.rho + -screened_rates(k_Co56_He4_to_n_Cu59_reaclib)*Y(Co56)*Y(He4)*state.rho) +
        (-screened_rates(k_Co56_to_Fe56_electron_capture_weaktab)*Y(Co56) + screened_rates(k_Fe56_to_Co56_beta_neg_weaktab)*Y(Fe56)) +
        (screened_rates(k_Ni56_to_Co56_electron_capture_weaktab)*Y(Ni56) + -screened_rates(k_Co56_to_Ni56_beta_neg_weaktab)*Y(Co56));

    ydot_nuc(Co57) =
        (screened_rates(k_Mn53_He4_to_Co57_reaclib)*Y(He4)*Y(Mn53)*state.rho + -screened_rates(k_Co57_to_He4_Mn53_reaclib)*Y(Co57)) +
        (screened_rates(k_Fe56_p_to_Co57_reaclib)*Y(Fe56)*Y(H1)*state.rho + -screened_rates(k_Co57_to_p_Fe56_reaclib)*Y(Co57)) +
        (screened_rates(k_Co56_n_to_Co57_reaclib)*Y(Co56)*Y(N)*state.rho + -screened_rates(k_Co57_to_n_Co56_reaclib)*Y(Co57)) +
        (-screened_rates(k_Co57_p_to_Ni58_reaclib)*Y(Co57)*Y(H1)*state.rho + screened_rates(k_Ni58_to_p_Co57_reaclib)*Y(Ni58)) +
        (-screened_rates(k_Co57_p_to_He4_Fe54_reaclib)*Y(Co57)*Y(H1)*state.rho + screened_rates(k_Fe54_He4_to_p_Co57_reaclib)*Y(Fe54)*Y(He4)*state.rho) +
        (-screened_rates(k_Co57_He4_to_p_Ni60_reaclib)*Y(Co57)*Y(He4)*state.rho + screened_rates(k_Ni60_p_to_He4_Co57_reaclib)*Y(Ni60)*Y(H1)*state.rho) +
        (-0.5*screened_rates(k_Co57_n_n_to_Co59_approx)*Y(Co57)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + screened_rates(k_Co59_to_Co57_n_n_approx)*Y(Co59));

    ydot_nuc(Co59) =
        (screened_rates(k_Mn55_He4_to_Co59_reaclib)*Y(He4)*Y(Mn55)*state.rho + -screened_rates(k_Co59_to_He4_Mn55_reaclib)*Y(Co59)) +
        (screened_rates(k_Fe58_p_to_Co59_reaclib)*Y(Fe58)*Y(H1)*state.rho + -screened_rates(k_Co59_to_p_Fe58_reaclib)*Y(Co59)) +
        (-screened_rates(k_Co59_p_to_Ni60_reaclib)*Y(Co59)*Y(H1)*state.rho + screened_rates(k_Ni60_to_p_Co59_reaclib)*Y(Ni60)) +
        (-screened_rates(k_Co59_p_to_He4_Fe56_reaclib)*Y(Co59)*Y(H1)*state.rho + screened_rates(k_Fe56_He4_to_p_Co59_reaclib)*Y(Fe56)*Y(He4)*state.rho) +
        (0.5*screened_rates(k_Co57_n_n_to_Co59_approx)*Y(Co57)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_Co59_to_Co57_n_n_approx)*Y(Co59));

    ydot_nuc(Ni56) =
        (screened_rates(k_Fe52_He4_to_Ni56_reaclib)*Y(Fe52)*Y(He4)*state.rho + -screened_rates(k_Ni56_to_He4_Fe52_reaclib)*Y(Ni56)) +
        (screened_rates(k_Co55_p_to_Ni56_reaclib)*Y(Co55)*Y(H1)*state.rho + -screened_rates(k_Ni56_to_p_Co55_reaclib)*Y(Ni56)) +
        (-screened_rates(k_Ni56_He4_to_Zn60_reaclib)*Y(He4)*Y(Ni56)*state.rho + screened_rates(k_Zn60_to_He4_Ni56_reaclib)*Y(Zn60)) +
        (-screened_rates(k_Ni56_n_to_p_Co56_reaclib)*Y(Ni56)*Y(N)*state.rho + screened_rates(k_Co56_p_to_n_Ni56_reaclib)*Y(Co56)*Y(H1)*state.rho) +
        (screened_rates(k_Cu59_p_to_He4_Ni56_reaclib)*Y(Cu59)*Y(H1)*state.rho + -screened_rates(k_Ni56_He4_to_p_Cu59_reaclib)*Y(He4)*Y(Ni56)*state.rho) +
        (-0.5*screened_rates(k_Ni56_n_n_to_Ni58_approx)*Y(Ni56)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + screened_rates(k_Ni58_to_Ni56_n_n_approx)*Y(Ni58)) +
        (-screened_rates(k_Ni56_to_Co56_electron_capture_weaktab)*Y(Ni56) + screened_rates(k_Co56_to_Ni56_beta_neg_weaktab)*Y(Co56));

    ydot_nuc(Ni58) =
        (screened_rates(k_Fe54_He4_to_Ni58_reaclib)*Y(Fe54)*Y(He4)*state.rho + -screened_rates(k_Ni58_to_He4_Fe54_reaclib)*Y(Ni58)) +
        (screened_rates(k_Co57_p_to_Ni58_reaclib)*Y(Co57)*Y(H1)*state.rho + -screened_rates(k_Ni58_to_p_Co57_reaclib)*Y(Ni58)) +
        (-screened_rates(k_Ni58_p_to_Cu59_reaclib)*Y(Ni58)*Y(H1)*state.rho + screened_rates(k_Cu59_to_p_Ni58_reaclib)*Y(Cu59)) +
        (screened_rates(k_Co55_He4_to_p_Ni58_reaclib)*Y(Co55)*Y(He4)*state.rho + -screened_rates(k_Ni58_p_to_He4_Co55_reaclib)*Y(Ni58)*Y(H1)*state.rho) +
        (0.5*screened_rates(k_Ni56_n_n_to_Ni58_approx)*Y(Ni56)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_Ni58_to_Ni56_n_n_approx)*Y(Ni58)) +
        (-0.5*screened_rates(k_Ni58_n_n_to_Ni60_approx)*Y(Ni58)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + screened_rates(k_Ni60_to_Ni58_n_n_approx)*Y(Ni60));

    ydot_nuc(Ni60) =
        (screened_rates(k_Fe56_He4_to_Ni60_reaclib)*Y(Fe56)*Y(He4)*state.rho + -screened_rates(k_Ni60_to_He4_Fe56_reaclib)*Y(Ni60)) +
        (screened_rates(k_Co59_p_to_Ni60_reaclib)*Y(Co59)*Y(H1)*state.rho + -screened_rates(k_Ni60_to_p_Co59_reaclib)*Y(Ni60)) +
        (screened_rates(k_Co57_He4_to_p_Ni60_reaclib)*Y(Co57)*Y(He4)*state.rho + -screened_rates(k_Ni60_p_to_He4_Co57_reaclib)*Y(Ni60)*Y(H1)*state.rho) +
        (0.5*screened_rates(k_Ni58_n_n_to_Ni60_approx)*Y(Ni58)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho) + -screened_rates(k_Ni60_to_Ni58_n_n_approx)*Y(Ni60));

    ydot_nuc(Cu59) =
        (screened_rates(k_Co55_He4_to_Cu59_reaclib)*Y(Co55)*Y(He4)*state.rho + -screened_rates(k_Cu59_to_He4_Co55_reaclib)*Y(Cu59)) +
        (screened_rates(k_Ni58_p_to_Cu59_reaclib)*Y(Ni58)*Y(H1)*state.rho + -screened_rates(k_Cu59_to_p_Ni58_reaclib)*Y(Cu59)) +
        (-screened_rates(k_Cu59_p_to_Zn60_reaclib)*Y(Cu59)*Y(H1)*state.rho + screened_rates(k_Zn60_to_p_Cu59_reaclib)*Y(Zn60)) +
        (-screened_rates(k_Cu59_n_to_He4_Co56_reaclib)*Y(Cu59)*Y(N)*state.rho + screened_rates(k_Co56_He4_to_n_Cu59_reaclib)*Y(Co56)*Y(He4)*state.rho) +
        (-screened_rates(k_Cu59_p_to_He4_Ni56_reaclib)*Y(Cu59)*Y(H1)*state.rho + screened_rates(k_Ni56_He4_to_p_Cu59_reaclib)*Y(He4)*Y(Ni56)*state.rho);

    ydot_nuc(Zn60) =
        (screened_rates(k_Ni56_He4_to_Zn60_reaclib)*Y(He4)*Y(Ni56)*state.rho + -screened_rates(k_Zn60_to_He4_Ni56_reaclib)*Y(Zn60)) +
        (screened_rates(k_Cu59_p_to_Zn60_reaclib)*Y(Cu59)*Y(H1)*state.rho + -screened_rates(k_Zn60_to_p_Cu59_reaclib)*Y(Zn60));

}


AMREX_GPU_DEVICE AMREX_INLINE
void actual_rhs (const burn_t& state, amrex::Array1D<amrex::Real, 1, neqs>& ydot)
{

    // Set molar abundances
    amrex::Array1D<amrex::Real, 1, NumSpec> Y;
    for (int i = 1; i <= NumSpec; ++i) {
        Y(i) = state.xn[i-1] * aion_inv[i-1];
    }

    // build the rates

    rate_t rate_eval;

    constexpr int do_T_derivatives = 0;

    evaluate_rates<do_T_derivatives, rate_t>(state, Y, rate_eval);

    rhs_nuc(state, ydot, Y, rate_eval.screened_rates);

    // ion binding energy contributions

    amrex::Real enuc;
    ener_gener_rate(ydot, enuc);

    // include any weak rate neutrino losses
    enuc += rate_eval.enuc_weak;

    // Get the thermal neutrino losses

    amrex::Real sneut{};
#ifdef NEUTRINOS
    constexpr int do_derivatives{0};
    amrex::Real dsneutdt{}, dsneutdd{}, dsnuda{}, dsnudz{};

    neutrino_cooling<do_derivatives>(state.T, state.rho, state.abar, state.zbar,
                                     sneut, dsneutdt, dsneutdd, dsnuda, dsnudz);
#endif

    // Append the energy equation (this is erg/g/s)

    ydot(net_ienuc) = enuc - sneut;

}


template<class MatrixType>
AMREX_GPU_DEVICE AMREX_INLINE
void jac_nuc(const burn_t& state,
             MatrixType& jac,
             const amrex::Array1D<amrex::Real, 1, NumSpec>& Y,
             const rate_derivs_t& rate_eval)
{

    using namespace Rates;

    const auto& screened_rates = rate_eval.screened_rates;

    amrex::Real scratch;

    // n row

    scratch = -screened_rates(k_Ar36_n_to_Ar37_reaclib)*Y(Ar36)*state.rho
              - screened_rates(k_Ar36_n_to_He4_S33_reaclib)*Y(Ar36)*state.rho
              - screened_rates(k_Ar37_n_to_Ar38_reaclib)*Y(Ar37)*state.rho
              - screened_rates(k_Ar37_n_to_He4_S34_reaclib)*Y(Ar37)*state.rho
              - screened_rates(k_Ar37_n_to_p_Cl37_reaclib)*Y(Ar37)*state.rho
              - screened_rates(k_Ar38_n_to_He4_S35_reaclib)*Y(Ar38)*state.rho
              - 2.0*screened_rates(k_Ca40_n_n_to_Ca42_approx)*Y(Ca40)*Y(N)*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_Ca40_n_to_He4_Ar37_reaclib)*Y(Ca40)*state.rho
              - screened_rates(k_Cl33_n_to_p_S33_reaclib)*Y(Cl33)*state.rho
              - screened_rates(k_Cl35_n_to_p_S35_reaclib)*Y(Cl35)*state.rho
              - screened_rates(k_Co55_n_to_Co56_reaclib)*Y(Co55)*state.rho
              - screened_rates(k_Co56_n_to_Co57_reaclib)*Y(Co56)*state.rho
              - screened_rates(k_Co56_n_to_He4_Mn53_reaclib)*Y(Co56)*state.rho
              - screened_rates(k_Co56_n_to_p_Fe56_reaclib)*Y(Co56)*state.rho
              - 2.0*screened_rates(k_Co57_n_n_to_Co59_approx)*Y(Co57)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 2.0*screened_rates(k_Cr48_n_n_to_Cr50_approx)*Y(Cr48)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 2.0*screened_rates(k_Cr50_n_n_to_Cr52_approx)*Y(Cr50)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 2.0*screened_rates(k_Cr52_n_n_to_Cr54_approx)*Y(Cr52)*Y(N)*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_Cu59_n_to_He4_Co56_reaclib)*Y(Cu59)*state.rho
              - 2.0*screened_rates(k_Fe52_n_n_to_Fe54_approx)*Y(Fe52)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 2.0*screened_rates(k_Fe54_n_n_to_Fe56_approx)*Y(Fe54)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 2.0*screened_rates(k_Fe56_n_n_to_Fe58_approx)*Y(Fe56)*Y(N)*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_K37_n_to_p_Ar37_reaclib)*Y(K37)*state.rho
              - screened_rates(k_Mg23_n_to_C12_C12_reaclib)*Y(Mg23)*state.rho
              - screened_rates(k_Mg23_n_to_He4_Ne20_reaclib)*Y(Mg23)*state.rho
              - screened_rates(k_Mg23_n_to_Mg24_reaclib)*Y(Mg23)*state.rho
              - 2.0*screened_rates(k_Mn51_n_n_to_Mn53_approx)*Y(Mn51)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 2.0*screened_rates(k_Mn53_n_n_to_Mn55_approx)*Y(Mn53)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 2.0*screened_rates(k_Ni56_n_n_to_Ni58_approx)*Y(Ni56)*Y(N)*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_Ni56_n_to_p_Co56_reaclib)*Y(Ni56)*state.rho
              - 2.0*screened_rates(k_Ni58_n_n_to_Ni60_approx)*Y(Ni58)*Y(N)*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_P29_n_to_p_Si29_reaclib)*Y(P29)*state.rho
              - screened_rates(k_P31_n_to_p_Si31_reaclib)*Y(P31)*state.rho
              - screened_rates(k_S32_n_to_He4_Si29_reaclib)*Y(S32)*state.rho
              - screened_rates(k_S32_n_to_S33_reaclib)*Y(S32)*state.rho
              - screened_rates(k_S33_n_to_He4_Si30_reaclib)*Y(S33)*state.rho
              - screened_rates(k_S33_n_to_S34_reaclib)*Y(S33)*state.rho
              - screened_rates(k_S33_n_to_p_P33_reaclib)*Y(S33)*state.rho
              - screened_rates(k_S34_n_to_He4_Si31_reaclib)*Y(S34)*state.rho
              - screened_rates(k_S34_n_to_S35_reaclib)*Y(S34)*state.rho
              - 2.0*screened_rates(k_Sc41_n_n_to_Sc43_approx)*Y(Sc41)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 2.0*screened_rates(k_Sc43_n_n_to_Sc45_approx)*Y(Sc43)*Y(N)*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_Si27_n_to_C12_O16_reaclib)*Y(Si27)*state.rho
              - screened_rates(k_Si27_n_to_He4_Mg24_reaclib)*Y(Si27)*state.rho
              - screened_rates(k_Si27_n_to_Si28_reaclib)*Y(Si27)*state.rho
              - screened_rates(k_Si27_n_to_p_Al27_reaclib)*Y(Si27)*state.rho
              - screened_rates(k_Si28_n_to_Si29_reaclib)*Y(Si28)*state.rho
              - screened_rates(k_Si29_n_to_Si30_reaclib)*Y(Si29)*state.rho
              - screened_rates(k_Si30_n_to_Si31_reaclib)*Y(Si30)*state.rho
              - 2.0*screened_rates(k_Ti44_n_n_to_Ti46_approx)*Y(Ti44)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 2.0*screened_rates(k_Ti46_n_n_to_Ti48_approx)*Y(Ti46)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 2.0*screened_rates(k_V47_n_n_to_V49_approx)*Y(V47)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 2.0*screened_rates(k_V49_n_n_to_V51_approx)*Y(V49)*Y(N)*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_n_to_p_beta_neg_weaktab)
              - Y(Ca40)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Ca40_n_n_to_Ca42_approx_dYN
              + 2.0*Y(Ca42)*rate_eval.drate_Ca42_to_Ca40_n_n_approx_dYN
              - Y(Co57)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Co57_n_n_to_Co59_approx_dYN
              + 2.0*Y(Co59)*rate_eval.drate_Co59_to_Co57_n_n_approx_dYN
              - Y(Cr48)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Cr48_n_n_to_Cr50_approx_dYN
              - Y(Cr50)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Cr50_n_n_to_Cr52_approx_dYN
              + 2.0*Y(Cr50)*rate_eval.drate_Cr50_to_Cr48_n_n_approx_dYN
              - Y(Cr52)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Cr52_n_n_to_Cr54_approx_dYN
              + 2.0*Y(Cr52)*rate_eval.drate_Cr52_to_Cr50_n_n_approx_dYN
              + 2.0*Y(Cr54)*rate_eval.drate_Cr54_to_Cr52_n_n_approx_dYN
              - Y(Fe52)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Fe52_n_n_to_Fe54_approx_dYN
              - Y(Fe54)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Fe54_n_n_to_Fe56_approx_dYN
              + 2.0*Y(Fe54)*rate_eval.drate_Fe54_to_Fe52_n_n_approx_dYN
              - Y(Fe56)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Fe56_n_n_to_Fe58_approx_dYN
              + 2.0*Y(Fe56)*rate_eval.drate_Fe56_to_Fe54_n_n_approx_dYN
              + 2.0*Y(Fe58)*rate_eval.drate_Fe58_to_Fe56_n_n_approx_dYN
              - Y(Mn51)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Mn51_n_n_to_Mn53_approx_dYN
              - Y(Mn53)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Mn53_n_n_to_Mn55_approx_dYN
              + 2.0*Y(Mn53)*rate_eval.drate_Mn53_to_Mn51_n_n_approx_dYN
              + 2.0*Y(Mn55)*rate_eval.drate_Mn55_to_Mn53_n_n_approx_dYN
              - Y(Ni56)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Ni56_n_n_to_Ni58_approx_dYN
              - Y(Ni58)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Ni58_n_n_to_Ni60_approx_dYN
              + 2.0*Y(Ni58)*rate_eval.drate_Ni58_to_Ni56_n_n_approx_dYN
              + 2.0*Y(Ni60)*rate_eval.drate_Ni60_to_Ni58_n_n_approx_dYN
              - Y(Sc41)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Sc41_n_n_to_Sc43_approx_dYN
              - Y(Sc43)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Sc43_n_n_to_Sc45_approx_dYN
              + 2.0*Y(Sc43)*rate_eval.drate_Sc43_to_Sc41_n_n_approx_dYN
              + 2.0*Y(Sc45)*rate_eval.drate_Sc45_to_Sc43_n_n_approx_dYN
              - Y(Ti44)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Ti44_n_n_to_Ti46_approx_dYN
              - Y(Ti46)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Ti46_n_n_to_Ti48_approx_dYN
              + 2.0*Y(Ti46)*rate_eval.drate_Ti46_to_Ti44_n_n_approx_dYN
              + 2.0*Y(Ti48)*rate_eval.drate_Ti48_to_Ti46_n_n_approx_dYN
              - Y(V47)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_V47_n_n_to_V49_approx_dYN
              - Y(V49)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_V49_n_n_to_V51_approx_dYN
              + 2.0*Y(V49)*rate_eval.drate_V49_to_V47_n_n_approx_dYN
              + 2.0*Y(V51)*rate_eval.drate_V51_to_V49_n_n_approx_dYN;
    jac.set(N, N, scratch);

    scratch = screened_rates(k_Al27_p_to_n_Si27_reaclib)*Y(Al27)*state.rho
              + screened_rates(k_Ar37_p_to_n_K37_reaclib)*Y(Ar37)*state.rho
              + screened_rates(k_Cl37_p_to_n_Ar37_reaclib)*Y(Cl37)*state.rho
              + screened_rates(k_Co56_p_to_n_Ni56_reaclib)*Y(Co56)*state.rho
              + screened_rates(k_Fe56_p_to_n_Co56_reaclib)*Y(Fe56)*state.rho
              + screened_rates(k_P33_p_to_n_S33_reaclib)*Y(P33)*state.rho
              + screened_rates(k_S33_p_to_n_Cl33_reaclib)*Y(S33)*state.rho
              + screened_rates(k_S35_p_to_n_Cl35_reaclib)*Y(S35)*state.rho
              + screened_rates(k_Si29_p_to_n_P29_reaclib)*Y(Si29)*state.rho
              + screened_rates(k_Si31_p_to_n_P31_reaclib)*Y(Si31)*state.rho
              + screened_rates(k_p_to_n_electron_capture_weaktab);
    jac.set(N, H1, scratch);

    scratch = screened_rates(k_Ar37_He4_to_n_Ca40_reaclib)*Y(Ar37)*state.rho
              + screened_rates(k_Co56_He4_to_n_Cu59_reaclib)*Y(Co56)*state.rho
              + screened_rates(k_Mg24_He4_to_n_Si27_reaclib)*Y(Mg24)*state.rho
              + screened_rates(k_Mn53_He4_to_n_Co56_reaclib)*Y(Mn53)*state.rho
              + screened_rates(k_Ne20_He4_to_n_Mg23_reaclib)*Y(Ne20)*state.rho
              + screened_rates(k_S33_He4_to_n_Ar36_reaclib)*Y(S33)*state.rho
              + screened_rates(k_S34_He4_to_n_Ar37_reaclib)*Y(S34)*state.rho
              + screened_rates(k_S35_He4_to_n_Ar38_reaclib)*Y(S35)*state.rho
              + screened_rates(k_Si29_He4_to_n_S32_reaclib)*Y(Si29)*state.rho
              + screened_rates(k_Si30_He4_to_n_S33_reaclib)*Y(Si30)*state.rho
              + screened_rates(k_Si31_He4_to_n_S34_reaclib)*Y(Si31)*state.rho;
    jac.set(N, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_n_Mg23_reaclib)*Y(C12)*state.rho
              + screened_rates(k_O16_C12_to_n_Si27_reaclib)*Y(O16)*state.rho;
    jac.set(N, C12, scratch);

    scratch = screened_rates(k_O16_C12_to_n_Si27_reaclib)*Y(C12)*state.rho;
    jac.set(N, O16, scratch);

    scratch = screened_rates(k_Ne20_He4_to_n_Mg23_reaclib)*Y(He4)*state.rho;
    jac.set(N, Ne20, scratch);

    scratch = -screened_rates(k_Mg23_n_to_C12_C12_reaclib)*Y(N)*state.rho
              - screened_rates(k_Mg23_n_to_He4_Ne20_reaclib)*Y(N)*state.rho
              - screened_rates(k_Mg23_n_to_Mg24_reaclib)*Y(N)*state.rho;
    jac.set(N, Mg23, scratch);

    scratch = screened_rates(k_Mg24_He4_to_n_Si27_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Mg24_to_n_Mg23_reaclib);
    jac.set(N, Mg24, scratch);

    scratch = screened_rates(k_Al27_p_to_n_Si27_reaclib)*Y(H1)*state.rho;
    jac.set(N, Al27, scratch);

    scratch = -screened_rates(k_Si27_n_to_C12_O16_reaclib)*Y(N)*state.rho
              - screened_rates(k_Si27_n_to_He4_Mg24_reaclib)*Y(N)*state.rho
              - screened_rates(k_Si27_n_to_Si28_reaclib)*Y(N)*state.rho
              - screened_rates(k_Si27_n_to_p_Al27_reaclib)*Y(N)*state.rho;
    jac.set(N, Si27, scratch);

    scratch = -screened_rates(k_Si28_n_to_Si29_reaclib)*Y(N)*state.rho
              + screened_rates(k_Si28_to_n_Si27_reaclib);
    jac.set(N, Si28, scratch);

    scratch = screened_rates(k_Si29_He4_to_n_S32_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Si29_n_to_Si30_reaclib)*Y(N)*state.rho
              + screened_rates(k_Si29_p_to_n_P29_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Si29_to_n_Si28_reaclib);
    jac.set(N, Si29, scratch);

    scratch = screened_rates(k_Si30_He4_to_n_S33_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Si30_n_to_Si31_reaclib)*Y(N)*state.rho
              + screened_rates(k_Si30_to_n_Si29_reaclib);
    jac.set(N, Si30, scratch);

    scratch = screened_rates(k_Si31_He4_to_n_S34_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Si31_p_to_n_P31_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Si31_to_n_Si30_reaclib);
    jac.set(N, Si31, scratch);

    scratch = -screened_rates(k_P29_n_to_p_Si29_reaclib)*Y(N)*state.rho;
    jac.set(N, P29, scratch);

    scratch = -screened_rates(k_P31_n_to_p_Si31_reaclib)*Y(N)*state.rho;
    jac.set(N, P31, scratch);

    scratch = screened_rates(k_P33_p_to_n_S33_reaclib)*Y(H1)*state.rho;
    jac.set(N, P33, scratch);

    scratch = -screened_rates(k_S32_n_to_He4_Si29_reaclib)*Y(N)*state.rho
              - screened_rates(k_S32_n_to_S33_reaclib)*Y(N)*state.rho;
    jac.set(N, S32, scratch);

    scratch = screened_rates(k_S33_He4_to_n_Ar36_reaclib)*Y(He4)*state.rho
              - screened_rates(k_S33_n_to_He4_Si30_reaclib)*Y(N)*state.rho
              - screened_rates(k_S33_n_to_S34_reaclib)*Y(N)*state.rho
              - screened_rates(k_S33_n_to_p_P33_reaclib)*Y(N)*state.rho
              + screened_rates(k_S33_p_to_n_Cl33_reaclib)*Y(H1)*state.rho
              + screened_rates(k_S33_to_n_S32_reaclib);
    jac.set(N, S33, scratch);

    scratch = screened_rates(k_S34_He4_to_n_Ar37_reaclib)*Y(He4)*state.rho
              - screened_rates(k_S34_n_to_He4_Si31_reaclib)*Y(N)*state.rho
              - screened_rates(k_S34_n_to_S35_reaclib)*Y(N)*state.rho
              + screened_rates(k_S34_to_n_S33_reaclib);
    jac.set(N, S34, scratch);

    scratch = screened_rates(k_S35_He4_to_n_Ar38_reaclib)*Y(He4)*state.rho
              + screened_rates(k_S35_p_to_n_Cl35_reaclib)*Y(H1)*state.rho
              + screened_rates(k_S35_to_n_S34_reaclib);
    jac.set(N, S35, scratch);

    scratch = -screened_rates(k_Cl33_n_to_p_S33_reaclib)*Y(N)*state.rho;
    jac.set(N, Cl33, scratch);

    scratch = -screened_rates(k_Cl35_n_to_p_S35_reaclib)*Y(N)*state.rho;
    jac.set(N, Cl35, scratch);

    scratch = screened_rates(k_Cl37_p_to_n_Ar37_reaclib)*Y(H1)*state.rho;
    jac.set(N, Cl37, scratch);

    scratch = -screened_rates(k_Ar36_n_to_Ar37_reaclib)*Y(N)*state.rho
              - screened_rates(k_Ar36_n_to_He4_S33_reaclib)*Y(N)*state.rho;
    jac.set(N, Ar36, scratch);

    scratch = screened_rates(k_Ar37_He4_to_n_Ca40_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ar37_n_to_Ar38_reaclib)*Y(N)*state.rho
              - screened_rates(k_Ar37_n_to_He4_S34_reaclib)*Y(N)*state.rho
              - screened_rates(k_Ar37_n_to_p_Cl37_reaclib)*Y(N)*state.rho
              + screened_rates(k_Ar37_p_to_n_K37_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Ar37_to_n_Ar36_reaclib);
    jac.set(N, Ar37, scratch);

    scratch = -screened_rates(k_Ar38_n_to_He4_S35_reaclib)*Y(N)*state.rho
              + screened_rates(k_Ar38_to_n_Ar37_reaclib);
    jac.set(N, Ar38, scratch);

    scratch = -screened_rates(k_K37_n_to_p_Ar37_reaclib)*Y(N)*state.rho;
    jac.set(N, K37, scratch);

    jac.set(N, K39, 0.0);

    scratch = -screened_rates(k_Ca40_n_n_to_Ca42_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_Ca40_n_to_He4_Ar37_reaclib)*Y(N)*state.rho;
    jac.set(N, Ca40, scratch);

    scratch = 2.0*screened_rates(k_Ca42_to_Ca40_n_n_approx);
    jac.set(N, Ca42, scratch);

    scratch = -screened_rates(k_Sc41_n_n_to_Sc43_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(N, Sc41, scratch);

    scratch = -screened_rates(k_Sc43_n_n_to_Sc45_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              + 2.0*screened_rates(k_Sc43_to_Sc41_n_n_approx);
    jac.set(N, Sc43, scratch);

    scratch = 2.0*screened_rates(k_Sc45_to_Sc43_n_n_approx);
    jac.set(N, Sc45, scratch);

    scratch = -screened_rates(k_Ti44_n_n_to_Ti46_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(N, Ti44, scratch);

    scratch = -screened_rates(k_Ti46_n_n_to_Ti48_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              + 2.0*screened_rates(k_Ti46_to_Ti44_n_n_approx);
    jac.set(N, Ti46, scratch);

    scratch = 2.0*screened_rates(k_Ti48_to_Ti46_n_n_approx);
    jac.set(N, Ti48, scratch);

    scratch = -screened_rates(k_V47_n_n_to_V49_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(N, V47, scratch);

    scratch = -screened_rates(k_V49_n_n_to_V51_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              + 2.0*screened_rates(k_V49_to_V47_n_n_approx);
    jac.set(N, V49, scratch);

    scratch = 2.0*screened_rates(k_V51_to_V49_n_n_approx);
    jac.set(N, V51, scratch);

    scratch = -screened_rates(k_Cr48_n_n_to_Cr50_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(N, Cr48, scratch);

    scratch = -screened_rates(k_Cr50_n_n_to_Cr52_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              + 2.0*screened_rates(k_Cr50_to_Cr48_n_n_approx);
    jac.set(N, Cr50, scratch);

    scratch = -screened_rates(k_Cr52_n_n_to_Cr54_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              + 2.0*screened_rates(k_Cr52_to_Cr50_n_n_approx);
    jac.set(N, Cr52, scratch);

    scratch = 2.0*screened_rates(k_Cr54_to_Cr52_n_n_approx);
    jac.set(N, Cr54, scratch);

    scratch = -screened_rates(k_Mn51_n_n_to_Mn53_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(N, Mn51, scratch);

    scratch = screened_rates(k_Mn53_He4_to_n_Co56_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Mn53_n_n_to_Mn55_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              + 2.0*screened_rates(k_Mn53_to_Mn51_n_n_approx);
    jac.set(N, Mn53, scratch);

    scratch = 2.0*screened_rates(k_Mn55_to_Mn53_n_n_approx);
    jac.set(N, Mn55, scratch);

    scratch = -screened_rates(k_Fe52_n_n_to_Fe54_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(N, Fe52, scratch);

    scratch = -screened_rates(k_Fe54_n_n_to_Fe56_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              + 2.0*screened_rates(k_Fe54_to_Fe52_n_n_approx);
    jac.set(N, Fe54, scratch);

    scratch = -screened_rates(k_Fe56_n_n_to_Fe58_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              + screened_rates(k_Fe56_p_to_n_Co56_reaclib)*Y(H1)*state.rho
              + 2.0*screened_rates(k_Fe56_to_Fe54_n_n_approx);
    jac.set(N, Fe56, scratch);

    scratch = 2.0*screened_rates(k_Fe58_to_Fe56_n_n_approx);
    jac.set(N, Fe58, scratch);

    scratch = -screened_rates(k_Co55_n_to_Co56_reaclib)*Y(N)*state.rho;
    jac.set(N, Co55, scratch);

    scratch = screened_rates(k_Co56_He4_to_n_Cu59_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Co56_n_to_Co57_reaclib)*Y(N)*state.rho
              - screened_rates(k_Co56_n_to_He4_Mn53_reaclib)*Y(N)*state.rho
              - screened_rates(k_Co56_n_to_p_Fe56_reaclib)*Y(N)*state.rho
              + screened_rates(k_Co56_p_to_n_Ni56_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Co56_to_n_Co55_reaclib);
    jac.set(N, Co56, scratch);

    scratch = -screened_rates(k_Co57_n_n_to_Co59_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              + screened_rates(k_Co57_to_n_Co56_reaclib);
    jac.set(N, Co57, scratch);

    scratch = 2.0*screened_rates(k_Co59_to_Co57_n_n_approx);
    jac.set(N, Co59, scratch);

    scratch = -screened_rates(k_Ni56_n_n_to_Ni58_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_Ni56_n_to_p_Co56_reaclib)*Y(N)*state.rho;
    jac.set(N, Ni56, scratch);

    scratch = -screened_rates(k_Ni58_n_n_to_Ni60_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              + 2.0*screened_rates(k_Ni58_to_Ni56_n_n_approx);
    jac.set(N, Ni58, scratch);

    scratch = 2.0*screened_rates(k_Ni60_to_Ni58_n_n_approx);
    jac.set(N, Ni60, scratch);

    scratch = -screened_rates(k_Cu59_n_to_He4_Co56_reaclib)*Y(N)*state.rho;
    jac.set(N, Cu59, scratch);

    jac.set(N, Zn60, 0.0);


    // p row

    scratch = screened_rates(k_Ar37_n_to_p_Cl37_reaclib)*Y(Ar37)*state.rho
              + screened_rates(k_Cl33_n_to_p_S33_reaclib)*Y(Cl33)*state.rho
              + screened_rates(k_Cl35_n_to_p_S35_reaclib)*Y(Cl35)*state.rho
              + screened_rates(k_Co56_n_to_p_Fe56_reaclib)*Y(Co56)*state.rho
              + screened_rates(k_K37_n_to_p_Ar37_reaclib)*Y(K37)*state.rho
              + screened_rates(k_Ni56_n_to_p_Co56_reaclib)*Y(Ni56)*state.rho
              + screened_rates(k_P29_n_to_p_Si29_reaclib)*Y(P29)*state.rho
              + screened_rates(k_P31_n_to_p_Si31_reaclib)*Y(P31)*state.rho
              + screened_rates(k_S33_n_to_p_P33_reaclib)*Y(S33)*state.rho
              + screened_rates(k_Si27_n_to_p_Al27_reaclib)*Y(Si27)*state.rho
              + screened_rates(k_n_to_p_beta_neg_weaktab);
    jac.set(H1, N, scratch);

    scratch = -screened_rates(k_Al27_p_to_C12_O16_reaclib)*Y(Al27)*state.rho
              - screened_rates(k_Al27_p_to_He4_Mg24_reaclib)*Y(Al27)*state.rho
              - screened_rates(k_Al27_p_to_Si28_reaclib)*Y(Al27)*state.rho
              - screened_rates(k_Al27_p_to_n_Si27_reaclib)*Y(Al27)*state.rho
              - screened_rates(k_Ar36_p_to_He4_Cl33_reaclib)*Y(Ar36)*state.rho
              - screened_rates(k_Ar36_p_to_K37_reaclib)*Y(Ar36)*state.rho
              - screened_rates(k_Ar37_p_to_n_K37_reaclib)*Y(Ar37)*state.rho
              - screened_rates(k_Ar38_p_to_He4_Cl35_reaclib)*Y(Ar38)*state.rho
              - screened_rates(k_Ar38_p_to_K39_reaclib)*Y(Ar38)*state.rho
              - screened_rates(k_Ca40_p_to_He4_K37_reaclib)*Y(Ca40)*state.rho
              - screened_rates(k_Ca40_p_to_Sc41_reaclib)*Y(Ca40)*state.rho
              - screened_rates(k_Ca42_p_to_He4_K39_reaclib)*Y(Ca42)*state.rho
              - screened_rates(k_Ca42_p_to_Sc43_reaclib)*Y(Ca42)*state.rho
              - screened_rates(k_Cl35_p_to_Ar36_reaclib)*Y(Cl35)*state.rho
              - screened_rates(k_Cl35_p_to_He4_S32_reaclib)*Y(Cl35)*state.rho
              - screened_rates(k_Cl37_p_to_Ar38_reaclib)*Y(Cl37)*state.rho
              - screened_rates(k_Cl37_p_to_He4_S34_reaclib)*Y(Cl37)*state.rho
              - screened_rates(k_Cl37_p_to_n_Ar37_reaclib)*Y(Cl37)*state.rho
              - screened_rates(k_Co55_p_to_He4_Fe52_reaclib)*Y(Co55)*state.rho
              - screened_rates(k_Co55_p_to_Ni56_reaclib)*Y(Co55)*state.rho
              - screened_rates(k_Co56_p_to_n_Ni56_reaclib)*Y(Co56)*state.rho
              - screened_rates(k_Co57_p_to_He4_Fe54_reaclib)*Y(Co57)*state.rho
              - screened_rates(k_Co57_p_to_Ni58_reaclib)*Y(Co57)*state.rho
              - screened_rates(k_Co59_p_to_He4_Fe56_reaclib)*Y(Co59)*state.rho
              - screened_rates(k_Co59_p_to_Ni60_reaclib)*Y(Co59)*state.rho
              - screened_rates(k_Cr50_p_to_He4_V47_reaclib)*Y(Cr50)*state.rho
              - screened_rates(k_Cr50_p_to_Mn51_reaclib)*Y(Cr50)*state.rho
              - screened_rates(k_Cr52_p_to_He4_V49_reaclib)*Y(Cr52)*state.rho
              - screened_rates(k_Cr52_p_to_Mn53_reaclib)*Y(Cr52)*state.rho
              - screened_rates(k_Cr54_p_to_He4_V51_reaclib)*Y(Cr54)*state.rho
              - screened_rates(k_Cr54_p_to_Mn55_reaclib)*Y(Cr54)*state.rho
              - screened_rates(k_Cu59_p_to_He4_Ni56_reaclib)*Y(Cu59)*state.rho
              - screened_rates(k_Cu59_p_to_Zn60_reaclib)*Y(Cu59)*state.rho
              - screened_rates(k_Fe54_p_to_Co55_reaclib)*Y(Fe54)*state.rho
              - screened_rates(k_Fe54_p_to_He4_Mn51_reaclib)*Y(Fe54)*state.rho
              - screened_rates(k_Fe56_p_to_Co57_reaclib)*Y(Fe56)*state.rho
              - screened_rates(k_Fe56_p_to_He4_Mn53_reaclib)*Y(Fe56)*state.rho
              - screened_rates(k_Fe56_p_to_n_Co56_reaclib)*Y(Fe56)*state.rho
              - screened_rates(k_Fe58_p_to_Co59_reaclib)*Y(Fe58)*state.rho
              - screened_rates(k_Fe58_p_to_He4_Mn55_reaclib)*Y(Fe58)*state.rho
              - screened_rates(k_K39_p_to_Ca40_reaclib)*Y(K39)*state.rho
              - screened_rates(k_K39_p_to_He4_Ar36_reaclib)*Y(K39)*state.rho
              - screened_rates(k_Mn51_p_to_Fe52_reaclib)*Y(Mn51)*state.rho
              - screened_rates(k_Mn51_p_to_He4_Cr48_reaclib)*Y(Mn51)*state.rho
              - screened_rates(k_Mn53_p_to_Fe54_reaclib)*Y(Mn53)*state.rho
              - screened_rates(k_Mn53_p_to_He4_Cr50_reaclib)*Y(Mn53)*state.rho
              - screened_rates(k_Mn55_p_to_Fe56_reaclib)*Y(Mn55)*state.rho
              - screened_rates(k_Mn55_p_to_He4_Cr52_reaclib)*Y(Mn55)*state.rho
              - screened_rates(k_Ni58_p_to_Cu59_reaclib)*Y(Ni58)*state.rho
              - screened_rates(k_Ni58_p_to_He4_Co55_reaclib)*Y(Ni58)*state.rho
              - screened_rates(k_Ni60_p_to_He4_Co57_reaclib)*Y(Ni60)*state.rho
              - screened_rates(k_P31_p_to_C12_Ne20_reaclib)*Y(P31)*state.rho
              - screened_rates(k_P31_p_to_He4_Si28_reaclib)*Y(P31)*state.rho
              - screened_rates(k_P31_p_to_O16_O16_reaclib)*Y(P31)*state.rho
              - screened_rates(k_P31_p_to_S32_reaclib)*Y(P31)*state.rho
              - screened_rates(k_P33_p_to_He4_Si30_reaclib)*Y(P33)*state.rho
              - screened_rates(k_P33_p_to_S34_reaclib)*Y(P33)*state.rho
              - screened_rates(k_P33_p_to_n_S33_reaclib)*Y(P33)*state.rho
              - screened_rates(k_S32_p_to_Cl33_reaclib)*Y(S32)*state.rho
              - screened_rates(k_S32_p_to_He4_P29_reaclib)*Y(S32)*state.rho
              - screened_rates(k_S33_p_to_n_Cl33_reaclib)*Y(S33)*state.rho
              - screened_rates(k_S34_p_to_Cl35_reaclib)*Y(S34)*state.rho
              - screened_rates(k_S34_p_to_He4_P31_reaclib)*Y(S34)*state.rho
              - screened_rates(k_S35_p_to_n_Cl35_reaclib)*Y(S35)*state.rho
              - screened_rates(k_Sc43_p_to_He4_Ca40_reaclib)*Y(Sc43)*state.rho
              - screened_rates(k_Sc43_p_to_Ti44_reaclib)*Y(Sc43)*state.rho
              - screened_rates(k_Sc45_p_to_He4_Ca42_reaclib)*Y(Sc45)*state.rho
              - screened_rates(k_Sc45_p_to_Ti46_reaclib)*Y(Sc45)*state.rho
              - screened_rates(k_Si28_p_to_P29_reaclib)*Y(Si28)*state.rho
              - screened_rates(k_Si29_p_to_n_P29_reaclib)*Y(Si29)*state.rho
              - screened_rates(k_Si30_p_to_He4_Al27_reaclib)*Y(Si30)*state.rho
              - screened_rates(k_Si30_p_to_P31_reaclib)*Y(Si30)*state.rho
              - screened_rates(k_Si31_p_to_n_P31_reaclib)*Y(Si31)*state.rho
              - screened_rates(k_Ti44_p_to_He4_Sc41_reaclib)*Y(Ti44)*state.rho
              - screened_rates(k_Ti46_p_to_He4_Sc43_reaclib)*Y(Ti46)*state.rho
              - screened_rates(k_Ti46_p_to_V47_reaclib)*Y(Ti46)*state.rho
              - screened_rates(k_Ti48_p_to_He4_Sc45_reaclib)*Y(Ti48)*state.rho
              - screened_rates(k_Ti48_p_to_V49_reaclib)*Y(Ti48)*state.rho
              - screened_rates(k_V47_p_to_Cr48_reaclib)*Y(V47)*state.rho
              - screened_rates(k_V47_p_to_He4_Ti44_reaclib)*Y(V47)*state.rho
              - screened_rates(k_V49_p_to_Cr50_reaclib)*Y(V49)*state.rho
              - screened_rates(k_V49_p_to_He4_Ti46_reaclib)*Y(V49)*state.rho
              - screened_rates(k_V51_p_to_Cr52_reaclib)*Y(V51)*state.rho
              - screened_rates(k_V51_p_to_He4_Ti48_reaclib)*Y(V51)*state.rho
              - screened_rates(k_p_to_n_electron_capture_weaktab);
    jac.set(H1, H1, scratch);

    scratch = screened_rates(k_Al27_He4_to_p_Si30_reaclib)*Y(Al27)*state.rho
              + screened_rates(k_Ar36_He4_to_p_K39_reaclib)*Y(Ar36)*state.rho
              + screened_rates(k_Ca40_He4_to_p_Sc43_reaclib)*Y(Ca40)*state.rho
              + screened_rates(k_Ca42_He4_to_p_Sc45_reaclib)*Y(Ca42)*state.rho
              + screened_rates(k_Cl33_He4_to_p_Ar36_reaclib)*Y(Cl33)*state.rho
              + screened_rates(k_Cl35_He4_to_p_Ar38_reaclib)*Y(Cl35)*state.rho
              + screened_rates(k_Co55_He4_to_p_Ni58_reaclib)*Y(Co55)*state.rho
              + screened_rates(k_Co57_He4_to_p_Ni60_reaclib)*Y(Co57)*state.rho
              + screened_rates(k_Cr48_He4_to_p_Mn51_reaclib)*Y(Cr48)*state.rho
              + screened_rates(k_Cr50_He4_to_p_Mn53_reaclib)*Y(Cr50)*state.rho
              + screened_rates(k_Cr52_He4_to_p_Mn55_reaclib)*Y(Cr52)*state.rho
              + screened_rates(k_Fe52_He4_to_p_Co55_reaclib)*Y(Fe52)*state.rho
              + screened_rates(k_Fe54_He4_to_p_Co57_reaclib)*Y(Fe54)*state.rho
              + screened_rates(k_Fe56_He4_to_p_Co59_reaclib)*Y(Fe56)*state.rho
              + screened_rates(k_K37_He4_to_p_Ca40_reaclib)*Y(K37)*state.rho
              + screened_rates(k_K39_He4_to_p_Ca42_reaclib)*Y(K39)*state.rho
              + screened_rates(k_Mg24_He4_to_p_Al27_reaclib)*Y(Mg24)*state.rho
              + screened_rates(k_Mn51_He4_to_p_Fe54_reaclib)*Y(Mn51)*state.rho
              + screened_rates(k_Mn53_He4_to_p_Fe56_reaclib)*Y(Mn53)*state.rho
              + screened_rates(k_Mn55_He4_to_p_Fe58_reaclib)*Y(Mn55)*state.rho
              + screened_rates(k_Ni56_He4_to_p_Cu59_reaclib)*Y(Ni56)*state.rho
              + screened_rates(k_P29_He4_to_p_S32_reaclib)*Y(P29)*state.rho
              + screened_rates(k_P31_He4_to_p_S34_reaclib)*Y(P31)*state.rho
              + screened_rates(k_S32_He4_to_p_Cl35_reaclib)*Y(S32)*state.rho
              + screened_rates(k_S34_He4_to_p_Cl37_reaclib)*Y(S34)*state.rho
              + screened_rates(k_Sc41_He4_to_p_Ti44_reaclib)*Y(Sc41)*state.rho
              + screened_rates(k_Sc43_He4_to_p_Ti46_reaclib)*Y(Sc43)*state.rho
              + screened_rates(k_Sc45_He4_to_p_Ti48_reaclib)*Y(Sc45)*state.rho
              + screened_rates(k_Si28_He4_to_p_P31_reaclib)*Y(Si28)*state.rho
              + screened_rates(k_Si30_He4_to_p_P33_reaclib)*Y(Si30)*state.rho
              + screened_rates(k_Ti44_He4_to_p_V47_reaclib)*Y(Ti44)*state.rho
              + screened_rates(k_Ti46_He4_to_p_V49_reaclib)*Y(Ti46)*state.rho
              + screened_rates(k_Ti48_He4_to_p_V51_reaclib)*Y(Ti48)*state.rho
              + screened_rates(k_V47_He4_to_p_Cr50_reaclib)*Y(V47)*state.rho
              + screened_rates(k_V49_He4_to_p_Cr52_reaclib)*Y(V49)*state.rho
              + screened_rates(k_V51_He4_to_p_Cr54_reaclib)*Y(V51)*state.rho;
    jac.set(H1, He4, scratch);

    scratch = screened_rates(k_Ne20_C12_to_p_P31_reaclib)*Y(Ne20)*state.rho
              + screened_rates(k_O16_C12_to_p_Al27_reaclib)*Y(O16)*state.rho;
    jac.set(H1, C12, scratch);

    scratch = screened_rates(k_O16_C12_to_p_Al27_reaclib)*Y(C12)*state.rho
              + 1.0*screened_rates(k_O16_O16_to_p_P31_reaclib)*Y(O16)*state.rho;
    jac.set(H1, O16, scratch);

    scratch = screened_rates(k_Ne20_C12_to_p_P31_reaclib)*Y(C12)*state.rho;
    jac.set(H1, Ne20, scratch);

    jac.set(H1, Mg23, 0.0);

    scratch = screened_rates(k_Mg24_He4_to_p_Al27_reaclib)*Y(He4)*state.rho;
    jac.set(H1, Mg24, scratch);

    scratch = screened_rates(k_Al27_He4_to_p_Si30_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Al27_p_to_C12_O16_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Al27_p_to_He4_Mg24_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Al27_p_to_Si28_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Al27_p_to_n_Si27_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Al27, scratch);

    scratch = screened_rates(k_Si27_n_to_p_Al27_reaclib)*Y(N)*state.rho;
    jac.set(H1, Si27, scratch);

    scratch = screened_rates(k_Si28_He4_to_p_P31_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Si28_p_to_P29_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Si28_to_p_Al27_reaclib);
    jac.set(H1, Si28, scratch);

    scratch = -screened_rates(k_Si29_p_to_n_P29_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Si29, scratch);

    scratch = screened_rates(k_Si30_He4_to_p_P33_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Si30_p_to_He4_Al27_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Si30_p_to_P31_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Si30, scratch);

    scratch = -screened_rates(k_Si31_p_to_n_P31_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Si31, scratch);

    scratch = screened_rates(k_P29_He4_to_p_S32_reaclib)*Y(He4)*state.rho
              + screened_rates(k_P29_n_to_p_Si29_reaclib)*Y(N)*state.rho
              + screened_rates(k_P29_to_p_Si28_reaclib);
    jac.set(H1, P29, scratch);

    scratch = screened_rates(k_P31_He4_to_p_S34_reaclib)*Y(He4)*state.rho
              + screened_rates(k_P31_n_to_p_Si31_reaclib)*Y(N)*state.rho
              - screened_rates(k_P31_p_to_C12_Ne20_reaclib)*Y(H1)*state.rho
              - screened_rates(k_P31_p_to_He4_Si28_reaclib)*Y(H1)*state.rho
              - screened_rates(k_P31_p_to_O16_O16_reaclib)*Y(H1)*state.rho
              - screened_rates(k_P31_p_to_S32_reaclib)*Y(H1)*state.rho
              + screened_rates(k_P31_to_p_Si30_reaclib);
    jac.set(H1, P31, scratch);

    scratch = -screened_rates(k_P33_p_to_He4_Si30_reaclib)*Y(H1)*state.rho
              - screened_rates(k_P33_p_to_S34_reaclib)*Y(H1)*state.rho
              - screened_rates(k_P33_p_to_n_S33_reaclib)*Y(H1)*state.rho;
    jac.set(H1, P33, scratch);

    scratch = screened_rates(k_S32_He4_to_p_Cl35_reaclib)*Y(He4)*state.rho
              - screened_rates(k_S32_p_to_Cl33_reaclib)*Y(H1)*state.rho
              - screened_rates(k_S32_p_to_He4_P29_reaclib)*Y(H1)*state.rho
              + screened_rates(k_S32_to_p_P31_reaclib);
    jac.set(H1, S32, scratch);

    scratch = screened_rates(k_S33_n_to_p_P33_reaclib)*Y(N)*state.rho
              - screened_rates(k_S33_p_to_n_Cl33_reaclib)*Y(H1)*state.rho;
    jac.set(H1, S33, scratch);

    scratch = screened_rates(k_S34_He4_to_p_Cl37_reaclib)*Y(He4)*state.rho
              - screened_rates(k_S34_p_to_Cl35_reaclib)*Y(H1)*state.rho
              - screened_rates(k_S34_p_to_He4_P31_reaclib)*Y(H1)*state.rho
              + screened_rates(k_S34_to_p_P33_reaclib);
    jac.set(H1, S34, scratch);

    scratch = -screened_rates(k_S35_p_to_n_Cl35_reaclib)*Y(H1)*state.rho;
    jac.set(H1, S35, scratch);

    scratch = screened_rates(k_Cl33_He4_to_p_Ar36_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Cl33_n_to_p_S33_reaclib)*Y(N)*state.rho
              + screened_rates(k_Cl33_to_p_S32_reaclib);
    jac.set(H1, Cl33, scratch);

    scratch = screened_rates(k_Cl35_He4_to_p_Ar38_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Cl35_n_to_p_S35_reaclib)*Y(N)*state.rho
              - screened_rates(k_Cl35_p_to_Ar36_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Cl35_p_to_He4_S32_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Cl35_to_p_S34_reaclib);
    jac.set(H1, Cl35, scratch);

    scratch = -screened_rates(k_Cl37_p_to_Ar38_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Cl37_p_to_He4_S34_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Cl37_p_to_n_Ar37_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Cl37, scratch);

    scratch = screened_rates(k_Ar36_He4_to_p_K39_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ar36_p_to_He4_Cl33_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Ar36_p_to_K37_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Ar36_to_p_Cl35_reaclib);
    jac.set(H1, Ar36, scratch);

    scratch = screened_rates(k_Ar37_n_to_p_Cl37_reaclib)*Y(N)*state.rho
              - screened_rates(k_Ar37_p_to_n_K37_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Ar37, scratch);

    scratch = -screened_rates(k_Ar38_p_to_He4_Cl35_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Ar38_p_to_K39_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Ar38_to_p_Cl37_reaclib);
    jac.set(H1, Ar38, scratch);

    scratch = screened_rates(k_K37_He4_to_p_Ca40_reaclib)*Y(He4)*state.rho
              + screened_rates(k_K37_n_to_p_Ar37_reaclib)*Y(N)*state.rho
              + screened_rates(k_K37_to_p_Ar36_reaclib);
    jac.set(H1, K37, scratch);

    scratch = screened_rates(k_K39_He4_to_p_Ca42_reaclib)*Y(He4)*state.rho
              - screened_rates(k_K39_p_to_Ca40_reaclib)*Y(H1)*state.rho
              - screened_rates(k_K39_p_to_He4_Ar36_reaclib)*Y(H1)*state.rho
              + screened_rates(k_K39_to_p_Ar38_reaclib);
    jac.set(H1, K39, scratch);

    scratch = screened_rates(k_Ca40_He4_to_p_Sc43_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ca40_p_to_He4_K37_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Ca40_p_to_Sc41_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Ca40_to_p_K39_reaclib);
    jac.set(H1, Ca40, scratch);

    scratch = screened_rates(k_Ca42_He4_to_p_Sc45_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ca42_p_to_He4_K39_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Ca42_p_to_Sc43_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Ca42, scratch);

    scratch = screened_rates(k_Sc41_He4_to_p_Ti44_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Sc41_to_p_Ca40_reaclib);
    jac.set(H1, Sc41, scratch);

    scratch = screened_rates(k_Sc43_He4_to_p_Ti46_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Sc43_p_to_He4_Ca40_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Sc43_p_to_Ti44_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Sc43_to_p_Ca42_reaclib);
    jac.set(H1, Sc43, scratch);

    scratch = screened_rates(k_Sc45_He4_to_p_Ti48_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Sc45_p_to_He4_Ca42_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Sc45_p_to_Ti46_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Sc45, scratch);

    scratch = screened_rates(k_Ti44_He4_to_p_V47_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ti44_p_to_He4_Sc41_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Ti44_to_p_Sc43_reaclib);
    jac.set(H1, Ti44, scratch);

    scratch = screened_rates(k_Ti46_He4_to_p_V49_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ti46_p_to_He4_Sc43_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Ti46_p_to_V47_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Ti46_to_p_Sc45_reaclib);
    jac.set(H1, Ti46, scratch);

    scratch = screened_rates(k_Ti48_He4_to_p_V51_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ti48_p_to_He4_Sc45_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Ti48_p_to_V49_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Ti48, scratch);

    scratch = screened_rates(k_V47_He4_to_p_Cr50_reaclib)*Y(He4)*state.rho
              - screened_rates(k_V47_p_to_Cr48_reaclib)*Y(H1)*state.rho
              - screened_rates(k_V47_p_to_He4_Ti44_reaclib)*Y(H1)*state.rho
              + screened_rates(k_V47_to_p_Ti46_reaclib);
    jac.set(H1, V47, scratch);

    scratch = screened_rates(k_V49_He4_to_p_Cr52_reaclib)*Y(He4)*state.rho
              - screened_rates(k_V49_p_to_Cr50_reaclib)*Y(H1)*state.rho
              - screened_rates(k_V49_p_to_He4_Ti46_reaclib)*Y(H1)*state.rho
              + screened_rates(k_V49_to_p_Ti48_reaclib);
    jac.set(H1, V49, scratch);

    scratch = screened_rates(k_V51_He4_to_p_Cr54_reaclib)*Y(He4)*state.rho
              - screened_rates(k_V51_p_to_Cr52_reaclib)*Y(H1)*state.rho
              - screened_rates(k_V51_p_to_He4_Ti48_reaclib)*Y(H1)*state.rho;
    jac.set(H1, V51, scratch);

    scratch = screened_rates(k_Cr48_He4_to_p_Mn51_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Cr48_to_p_V47_reaclib);
    jac.set(H1, Cr48, scratch);

    scratch = screened_rates(k_Cr50_He4_to_p_Mn53_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Cr50_p_to_He4_V47_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Cr50_p_to_Mn51_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Cr50_to_p_V49_reaclib);
    jac.set(H1, Cr50, scratch);

    scratch = screened_rates(k_Cr52_He4_to_p_Mn55_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Cr52_p_to_He4_V49_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Cr52_p_to_Mn53_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Cr52_to_p_V51_reaclib);
    jac.set(H1, Cr52, scratch);

    scratch = -screened_rates(k_Cr54_p_to_He4_V51_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Cr54_p_to_Mn55_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Cr54, scratch);

    scratch = screened_rates(k_Mn51_He4_to_p_Fe54_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Mn51_p_to_Fe52_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Mn51_p_to_He4_Cr48_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Mn51_to_p_Cr50_reaclib);
    jac.set(H1, Mn51, scratch);

    scratch = screened_rates(k_Mn53_He4_to_p_Fe56_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Mn53_p_to_Fe54_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Mn53_p_to_He4_Cr50_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Mn53_to_p_Cr52_reaclib);
    jac.set(H1, Mn53, scratch);

    scratch = screened_rates(k_Mn55_He4_to_p_Fe58_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Mn55_p_to_Fe56_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Mn55_p_to_He4_Cr52_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Mn55_to_p_Cr54_reaclib);
    jac.set(H1, Mn55, scratch);

    scratch = screened_rates(k_Fe52_He4_to_p_Co55_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Fe52_to_p_Mn51_reaclib);
    jac.set(H1, Fe52, scratch);

    scratch = screened_rates(k_Fe54_He4_to_p_Co57_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Fe54_p_to_Co55_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Fe54_p_to_He4_Mn51_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Fe54_to_p_Mn53_reaclib);
    jac.set(H1, Fe54, scratch);

    scratch = screened_rates(k_Fe56_He4_to_p_Co59_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Fe56_p_to_Co57_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Fe56_p_to_He4_Mn53_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Fe56_p_to_n_Co56_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Fe56_to_p_Mn55_reaclib);
    jac.set(H1, Fe56, scratch);

    scratch = -screened_rates(k_Fe58_p_to_Co59_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Fe58_p_to_He4_Mn55_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Fe58, scratch);

    scratch = screened_rates(k_Co55_He4_to_p_Ni58_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Co55_p_to_He4_Fe52_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Co55_p_to_Ni56_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Co55_to_p_Fe54_reaclib);
    jac.set(H1, Co55, scratch);

    scratch = screened_rates(k_Co56_n_to_p_Fe56_reaclib)*Y(N)*state.rho
              - screened_rates(k_Co56_p_to_n_Ni56_reaclib)*Y(H1)*state.rho;
    jac.set(H1, Co56, scratch);

    scratch = screened_rates(k_Co57_He4_to_p_Ni60_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Co57_p_to_He4_Fe54_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Co57_p_to_Ni58_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Co57_to_p_Fe56_reaclib);
    jac.set(H1, Co57, scratch);

    scratch = -screened_rates(k_Co59_p_to_He4_Fe56_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Co59_p_to_Ni60_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Co59_to_p_Fe58_reaclib);
    jac.set(H1, Co59, scratch);

    scratch = screened_rates(k_Ni56_He4_to_p_Cu59_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Ni56_n_to_p_Co56_reaclib)*Y(N)*state.rho
              + screened_rates(k_Ni56_to_p_Co55_reaclib);
    jac.set(H1, Ni56, scratch);

    scratch = -screened_rates(k_Ni58_p_to_Cu59_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Ni58_p_to_He4_Co55_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Ni58_to_p_Co57_reaclib);
    jac.set(H1, Ni58, scratch);

    scratch = -screened_rates(k_Ni60_p_to_He4_Co57_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Ni60_to_p_Co59_reaclib);
    jac.set(H1, Ni60, scratch);

    scratch = -screened_rates(k_Cu59_p_to_He4_Ni56_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Cu59_p_to_Zn60_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Cu59_to_p_Ni58_reaclib);
    jac.set(H1, Cu59, scratch);

    scratch = screened_rates(k_Zn60_to_p_Cu59_reaclib);
    jac.set(H1, Zn60, scratch);


    // He4 row

    scratch = screened_rates(k_Ar36_n_to_He4_S33_reaclib)*Y(Ar36)*state.rho
              + screened_rates(k_Ar37_n_to_He4_S34_reaclib)*Y(Ar37)*state.rho
              + screened_rates(k_Ar38_n_to_He4_S35_reaclib)*Y(Ar38)*state.rho
              + screened_rates(k_Ca40_n_to_He4_Ar37_reaclib)*Y(Ca40)*state.rho
              + screened_rates(k_Co56_n_to_He4_Mn53_reaclib)*Y(Co56)*state.rho
              + screened_rates(k_Cu59_n_to_He4_Co56_reaclib)*Y(Cu59)*state.rho
              + screened_rates(k_Mg23_n_to_He4_Ne20_reaclib)*Y(Mg23)*state.rho
              + screened_rates(k_S32_n_to_He4_Si29_reaclib)*Y(S32)*state.rho
              + screened_rates(k_S33_n_to_He4_Si30_reaclib)*Y(S33)*state.rho
              + screened_rates(k_S34_n_to_He4_Si31_reaclib)*Y(S34)*state.rho
              + screened_rates(k_Si27_n_to_He4_Mg24_reaclib)*Y(Si27)*state.rho;
    jac.set(He4, N, scratch);

    scratch = screened_rates(k_Al27_p_to_He4_Mg24_reaclib)*Y(Al27)*state.rho
              + screened_rates(k_Ar36_p_to_He4_Cl33_reaclib)*Y(Ar36)*state.rho
              + screened_rates(k_Ar38_p_to_He4_Cl35_reaclib)*Y(Ar38)*state.rho
              + screened_rates(k_Ca40_p_to_He4_K37_reaclib)*Y(Ca40)*state.rho
              + screened_rates(k_Ca42_p_to_He4_K39_reaclib)*Y(Ca42)*state.rho
              + screened_rates(k_Cl35_p_to_He4_S32_reaclib)*Y(Cl35)*state.rho
              + screened_rates(k_Cl37_p_to_He4_S34_reaclib)*Y(Cl37)*state.rho
              + screened_rates(k_Co55_p_to_He4_Fe52_reaclib)*Y(Co55)*state.rho
              + screened_rates(k_Co57_p_to_He4_Fe54_reaclib)*Y(Co57)*state.rho
              + screened_rates(k_Co59_p_to_He4_Fe56_reaclib)*Y(Co59)*state.rho
              + screened_rates(k_Cr50_p_to_He4_V47_reaclib)*Y(Cr50)*state.rho
              + screened_rates(k_Cr52_p_to_He4_V49_reaclib)*Y(Cr52)*state.rho
              + screened_rates(k_Cr54_p_to_He4_V51_reaclib)*Y(Cr54)*state.rho
              + screened_rates(k_Cu59_p_to_He4_Ni56_reaclib)*Y(Cu59)*state.rho
              + screened_rates(k_Fe54_p_to_He4_Mn51_reaclib)*Y(Fe54)*state.rho
              + screened_rates(k_Fe56_p_to_He4_Mn53_reaclib)*Y(Fe56)*state.rho
              + screened_rates(k_Fe58_p_to_He4_Mn55_reaclib)*Y(Fe58)*state.rho
              + screened_rates(k_K39_p_to_He4_Ar36_reaclib)*Y(K39)*state.rho
              + screened_rates(k_Mn51_p_to_He4_Cr48_reaclib)*Y(Mn51)*state.rho
              + screened_rates(k_Mn53_p_to_He4_Cr50_reaclib)*Y(Mn53)*state.rho
              + screened_rates(k_Mn55_p_to_He4_Cr52_reaclib)*Y(Mn55)*state.rho
              + screened_rates(k_Ni58_p_to_He4_Co55_reaclib)*Y(Ni58)*state.rho
              + screened_rates(k_Ni60_p_to_He4_Co57_reaclib)*Y(Ni60)*state.rho
              + screened_rates(k_P31_p_to_He4_Si28_reaclib)*Y(P31)*state.rho
              + screened_rates(k_P33_p_to_He4_Si30_reaclib)*Y(P33)*state.rho
              + screened_rates(k_S32_p_to_He4_P29_reaclib)*Y(S32)*state.rho
              + screened_rates(k_S34_p_to_He4_P31_reaclib)*Y(S34)*state.rho
              + screened_rates(k_Sc43_p_to_He4_Ca40_reaclib)*Y(Sc43)*state.rho
              + screened_rates(k_Sc45_p_to_He4_Ca42_reaclib)*Y(Sc45)*state.rho
              + screened_rates(k_Si30_p_to_He4_Al27_reaclib)*Y(Si30)*state.rho
              + screened_rates(k_Ti44_p_to_He4_Sc41_reaclib)*Y(Ti44)*state.rho
              + screened_rates(k_Ti46_p_to_He4_Sc43_reaclib)*Y(Ti46)*state.rho
              + screened_rates(k_Ti48_p_to_He4_Sc45_reaclib)*Y(Ti48)*state.rho
              + screened_rates(k_V47_p_to_He4_Ti44_reaclib)*Y(V47)*state.rho
              + screened_rates(k_V49_p_to_He4_Ti46_reaclib)*Y(V49)*state.rho
              + screened_rates(k_V51_p_to_He4_Ti48_reaclib)*Y(V51)*state.rho;
    jac.set(He4, H1, scratch);

    scratch = -screened_rates(k_Al27_He4_to_P31_reaclib)*Y(Al27)*state.rho
              - screened_rates(k_Al27_He4_to_p_Si30_reaclib)*Y(Al27)*state.rho
              - screened_rates(k_Ar36_He4_to_Ca40_reaclib)*Y(Ar36)*state.rho
              - screened_rates(k_Ar36_He4_to_p_K39_reaclib)*Y(Ar36)*state.rho
              - screened_rates(k_Ar37_He4_to_n_Ca40_reaclib)*Y(Ar37)*state.rho
              - screened_rates(k_Ar38_He4_to_Ca42_reaclib)*Y(Ar38)*state.rho
              - screened_rates(k_C12_He4_to_O16_reaclib)*Y(C12)*state.rho
              - screened_rates(k_Ca40_He4_to_Ti44_reaclib)*Y(Ca40)*state.rho
              - screened_rates(k_Ca40_He4_to_p_Sc43_reaclib)*Y(Ca40)*state.rho
              - screened_rates(k_Ca42_He4_to_Ti46_reaclib)*Y(Ca42)*state.rho
              - screened_rates(k_Ca42_He4_to_p_Sc45_reaclib)*Y(Ca42)*state.rho
              - screened_rates(k_Cl33_He4_to_K37_reaclib)*Y(Cl33)*state.rho
              - screened_rates(k_Cl33_He4_to_p_Ar36_reaclib)*Y(Cl33)*state.rho
              - screened_rates(k_Cl35_He4_to_K39_reaclib)*Y(Cl35)*state.rho
              - screened_rates(k_Cl35_He4_to_p_Ar38_reaclib)*Y(Cl35)*state.rho
              - screened_rates(k_Co55_He4_to_Cu59_reaclib)*Y(Co55)*state.rho
              - screened_rates(k_Co55_He4_to_p_Ni58_reaclib)*Y(Co55)*state.rho
              - screened_rates(k_Co56_He4_to_n_Cu59_reaclib)*Y(Co56)*state.rho
              - screened_rates(k_Co57_He4_to_p_Ni60_reaclib)*Y(Co57)*state.rho
              - screened_rates(k_Cr48_He4_to_Fe52_reaclib)*Y(Cr48)*state.rho
              - screened_rates(k_Cr48_He4_to_p_Mn51_reaclib)*Y(Cr48)*state.rho
              - screened_rates(k_Cr50_He4_to_Fe54_reaclib)*Y(Cr50)*state.rho
              - screened_rates(k_Cr50_He4_to_p_Mn53_reaclib)*Y(Cr50)*state.rho
              - screened_rates(k_Cr52_He4_to_Fe56_reaclib)*Y(Cr52)*state.rho
              - screened_rates(k_Cr52_He4_to_p_Mn55_reaclib)*Y(Cr52)*state.rho
              - screened_rates(k_Cr54_He4_to_Fe58_reaclib)*Y(Cr54)*state.rho
              - screened_rates(k_Fe52_He4_to_Ni56_reaclib)*Y(Fe52)*state.rho
              - screened_rates(k_Fe52_He4_to_p_Co55_reaclib)*Y(Fe52)*state.rho
              - screened_rates(k_Fe54_He4_to_Ni58_reaclib)*Y(Fe54)*state.rho
              - screened_rates(k_Fe54_He4_to_p_Co57_reaclib)*Y(Fe54)*state.rho
              - screened_rates(k_Fe56_He4_to_Ni60_reaclib)*Y(Fe56)*state.rho
              - screened_rates(k_Fe56_He4_to_p_Co59_reaclib)*Y(Fe56)*state.rho
              - 1.5*screened_rates(k_He4_He4_He4_to_C12_reaclib)*amrex::Math::powi<2>(Y(He4))*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_K37_He4_to_Sc41_reaclib)*Y(K37)*state.rho
              - screened_rates(k_K37_He4_to_p_Ca40_reaclib)*Y(K37)*state.rho
              - screened_rates(k_K39_He4_to_Sc43_reaclib)*Y(K39)*state.rho
              - screened_rates(k_K39_He4_to_p_Ca42_reaclib)*Y(K39)*state.rho
              - screened_rates(k_Mg23_He4_to_Si27_reaclib)*Y(Mg23)*state.rho
              - screened_rates(k_Mg24_He4_to_C12_O16_reaclib)*Y(Mg24)*state.rho
              - screened_rates(k_Mg24_He4_to_Si28_reaclib)*Y(Mg24)*state.rho
              - screened_rates(k_Mg24_He4_to_n_Si27_reaclib)*Y(Mg24)*state.rho
              - screened_rates(k_Mg24_He4_to_p_Al27_reaclib)*Y(Mg24)*state.rho
              - screened_rates(k_Mn51_He4_to_Co55_reaclib)*Y(Mn51)*state.rho
              - screened_rates(k_Mn51_He4_to_p_Fe54_reaclib)*Y(Mn51)*state.rho
              - screened_rates(k_Mn53_He4_to_Co57_reaclib)*Y(Mn53)*state.rho
              - screened_rates(k_Mn53_He4_to_n_Co56_reaclib)*Y(Mn53)*state.rho
              - screened_rates(k_Mn53_He4_to_p_Fe56_reaclib)*Y(Mn53)*state.rho
              - screened_rates(k_Mn55_He4_to_Co59_reaclib)*Y(Mn55)*state.rho
              - screened_rates(k_Mn55_He4_to_p_Fe58_reaclib)*Y(Mn55)*state.rho
              - screened_rates(k_Ne20_He4_to_C12_C12_approx)*Y(Ne20)*state.rho
              - screened_rates(k_Ne20_He4_to_Mg24_approx)*Y(Ne20)*state.rho
              - screened_rates(k_Ne20_He4_to_n_Mg23_reaclib)*Y(Ne20)*state.rho
              - screened_rates(k_Ni56_He4_to_Zn60_reaclib)*Y(Ni56)*state.rho
              - screened_rates(k_Ni56_He4_to_p_Cu59_reaclib)*Y(Ni56)*state.rho
              - screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(O16)*state.rho
              - screened_rates(k_P29_He4_to_Cl33_reaclib)*Y(P29)*state.rho
              - screened_rates(k_P29_He4_to_p_S32_reaclib)*Y(P29)*state.rho
              - screened_rates(k_P31_He4_to_Cl35_reaclib)*Y(P31)*state.rho
              - screened_rates(k_P31_He4_to_p_S34_reaclib)*Y(P31)*state.rho
              - screened_rates(k_P33_He4_to_Cl37_reaclib)*Y(P33)*state.rho
              - screened_rates(k_S32_He4_to_Ar36_reaclib)*Y(S32)*state.rho
              - screened_rates(k_S32_He4_to_p_Cl35_reaclib)*Y(S32)*state.rho
              - screened_rates(k_S33_He4_to_Ar37_reaclib)*Y(S33)*state.rho
              - screened_rates(k_S33_He4_to_n_Ar36_reaclib)*Y(S33)*state.rho
              - screened_rates(k_S34_He4_to_Ar38_reaclib)*Y(S34)*state.rho
              - screened_rates(k_S34_He4_to_n_Ar37_reaclib)*Y(S34)*state.rho
              - screened_rates(k_S34_He4_to_p_Cl37_reaclib)*Y(S34)*state.rho
              - screened_rates(k_S35_He4_to_n_Ar38_reaclib)*Y(S35)*state.rho
              - screened_rates(k_Sc41_He4_to_p_Ti44_reaclib)*Y(Sc41)*state.rho
              - screened_rates(k_Sc43_He4_to_V47_reaclib)*Y(Sc43)*state.rho
              - screened_rates(k_Sc43_He4_to_p_Ti46_reaclib)*Y(Sc43)*state.rho
              - screened_rates(k_Sc45_He4_to_V49_reaclib)*Y(Sc45)*state.rho
              - screened_rates(k_Sc45_He4_to_p_Ti48_reaclib)*Y(Sc45)*state.rho
              - screened_rates(k_Si28_He4_to_C12_Ne20_reaclib)*Y(Si28)*state.rho
              - screened_rates(k_Si28_He4_to_O16_O16_reaclib)*Y(Si28)*state.rho
              - screened_rates(k_Si28_He4_to_S32_reaclib)*Y(Si28)*state.rho
              - screened_rates(k_Si28_He4_to_p_P31_reaclib)*Y(Si28)*state.rho
              - screened_rates(k_Si29_He4_to_S33_reaclib)*Y(Si29)*state.rho
              - screened_rates(k_Si29_He4_to_n_S32_reaclib)*Y(Si29)*state.rho
              - screened_rates(k_Si30_He4_to_S34_reaclib)*Y(Si30)*state.rho
              - screened_rates(k_Si30_He4_to_n_S33_reaclib)*Y(Si30)*state.rho
              - screened_rates(k_Si30_He4_to_p_P33_reaclib)*Y(Si30)*state.rho
              - screened_rates(k_Si31_He4_to_S35_reaclib)*Y(Si31)*state.rho
              - screened_rates(k_Si31_He4_to_n_S34_reaclib)*Y(Si31)*state.rho
              - screened_rates(k_Ti44_He4_to_Cr48_reaclib)*Y(Ti44)*state.rho
              - screened_rates(k_Ti44_He4_to_p_V47_reaclib)*Y(Ti44)*state.rho
              - screened_rates(k_Ti46_He4_to_Cr50_reaclib)*Y(Ti46)*state.rho
              - screened_rates(k_Ti46_He4_to_p_V49_reaclib)*Y(Ti46)*state.rho
              - screened_rates(k_Ti48_He4_to_Cr52_reaclib)*Y(Ti48)*state.rho
              - screened_rates(k_Ti48_He4_to_p_V51_reaclib)*Y(Ti48)*state.rho
              - screened_rates(k_V47_He4_to_Mn51_reaclib)*Y(V47)*state.rho
              - screened_rates(k_V47_He4_to_p_Cr50_reaclib)*Y(V47)*state.rho
              - screened_rates(k_V49_He4_to_Mn53_reaclib)*Y(V49)*state.rho
              - screened_rates(k_V49_He4_to_p_Cr52_reaclib)*Y(V49)*state.rho
              - screened_rates(k_V51_He4_to_Mn55_reaclib)*Y(V51)*state.rho
              - screened_rates(k_V51_He4_to_p_Cr54_reaclib)*Y(V51)*state.rho;
    jac.set(He4, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_Ne20_He4_approx)*Y(C12)*state.rho
              - screened_rates(k_C12_He4_to_O16_reaclib)*Y(He4)*state.rho
              + 3.0*screened_rates(k_C12_to_3He4_reaclib)
              + screened_rates(k_Ne20_C12_to_He4_Si28_reaclib)*Y(Ne20)*state.rho
              + screened_rates(k_O16_C12_to_He4_Mg24_reaclib)*Y(O16)*state.rho;
    jac.set(He4, C12, scratch);

    scratch = screened_rates(k_O16_C12_to_He4_Mg24_reaclib)*Y(C12)*state.rho
              - screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(He4)*state.rho
              + 1.0*screened_rates(k_O16_O16_to_He4_Si28_reaclib)*Y(O16)*state.rho
              + screened_rates(k_O16_to_He4_C12_reaclib);
    jac.set(He4, O16, scratch);

    scratch = screened_rates(k_Ne20_C12_to_He4_Si28_reaclib)*Y(C12)*state.rho
              - screened_rates(k_Ne20_He4_to_C12_C12_approx)*Y(He4)*state.rho
              - screened_rates(k_Ne20_He4_to_Mg24_approx)*Y(He4)*state.rho
              - screened_rates(k_Ne20_He4_to_n_Mg23_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Ne20_to_He4_O16_reaclib);
    jac.set(He4, Ne20, scratch);

    scratch = -screened_rates(k_Mg23_He4_to_Si27_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Mg23_n_to_He4_Ne20_reaclib)*Y(N)*state.rho;
    jac.set(He4, Mg23, scratch);

    scratch = -screened_rates(k_Mg24_He4_to_C12_O16_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Mg24_He4_to_Si28_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Mg24_He4_to_n_Si27_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Mg24_He4_to_p_Al27_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Mg24_to_Ne20_He4_approx);
    jac.set(He4, Mg24, scratch);

    scratch = -screened_rates(k_Al27_He4_to_P31_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Al27_He4_to_p_Si30_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Al27_p_to_He4_Mg24_reaclib)*Y(H1)*state.rho;
    jac.set(He4, Al27, scratch);

    scratch = screened_rates(k_Si27_n_to_He4_Mg24_reaclib)*Y(N)*state.rho
              + screened_rates(k_Si27_to_He4_Mg23_reaclib);
    jac.set(He4, Si27, scratch);

    scratch = -screened_rates(k_Si28_He4_to_C12_Ne20_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Si28_He4_to_O16_O16_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Si28_He4_to_S32_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Si28_He4_to_p_P31_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Si28_to_He4_Mg24_reaclib);
    jac.set(He4, Si28, scratch);

    scratch = -screened_rates(k_Si29_He4_to_S33_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Si29_He4_to_n_S32_reaclib)*Y(He4)*state.rho;
    jac.set(He4, Si29, scratch);

    scratch = -screened_rates(k_Si30_He4_to_S34_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Si30_He4_to_n_S33_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Si30_He4_to_p_P33_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Si30_p_to_He4_Al27_reaclib)*Y(H1)*state.rho;
    jac.set(He4, Si30, scratch);

    scratch = -screened_rates(k_Si31_He4_to_S35_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Si31_He4_to_n_S34_reaclib)*Y(He4)*state.rho;
    jac.set(He4, Si31, scratch);

    scratch = -screened_rates(k_P29_He4_to_Cl33_reaclib)*Y(He4)*state.rho
              - screened_rates(k_P29_He4_to_p_S32_reaclib)*Y(He4)*state.rho;
    jac.set(He4, P29, scratch);

    scratch = -screened_rates(k_P31_He4_to_Cl35_reaclib)*Y(He4)*state.rho
              - screened_rates(k_P31_He4_to_p_S34_reaclib)*Y(He4)*state.rho
              + screened_rates(k_P31_p_to_He4_Si28_reaclib)*Y(H1)*state.rho
              + screened_rates(k_P31_to_He4_Al27_reaclib);
    jac.set(He4, P31, scratch);

    scratch = -screened_rates(k_P33_He4_to_Cl37_reaclib)*Y(He4)*state.rho
              + screened_rates(k_P33_p_to_He4_Si30_reaclib)*Y(H1)*state.rho;
    jac.set(He4, P33, scratch);

    scratch = -screened_rates(k_S32_He4_to_Ar36_reaclib)*Y(He4)*state.rho
              - screened_rates(k_S32_He4_to_p_Cl35_reaclib)*Y(He4)*state.rho
              + screened_rates(k_S32_n_to_He4_Si29_reaclib)*Y(N)*state.rho
              + screened_rates(k_S32_p_to_He4_P29_reaclib)*Y(H1)*state.rho
              + screened_rates(k_S32_to_He4_Si28_reaclib);
    jac.set(He4, S32, scratch);

    scratch = -screened_rates(k_S33_He4_to_Ar37_reaclib)*Y(He4)*state.rho
              - screened_rates(k_S33_He4_to_n_Ar36_reaclib)*Y(He4)*state.rho
              + screened_rates(k_S33_n_to_He4_Si30_reaclib)*Y(N)*state.rho
              + screened_rates(k_S33_to_He4_Si29_reaclib);
    jac.set(He4, S33, scratch);

    scratch = -screened_rates(k_S34_He4_to_Ar38_reaclib)*Y(He4)*state.rho
              - screened_rates(k_S34_He4_to_n_Ar37_reaclib)*Y(He4)*state.rho
              - screened_rates(k_S34_He4_to_p_Cl37_reaclib)*Y(He4)*state.rho
              + screened_rates(k_S34_n_to_He4_Si31_reaclib)*Y(N)*state.rho
              + screened_rates(k_S34_p_to_He4_P31_reaclib)*Y(H1)*state.rho
              + screened_rates(k_S34_to_He4_Si30_reaclib);
    jac.set(He4, S34, scratch);

    scratch = -screened_rates(k_S35_He4_to_n_Ar38_reaclib)*Y(He4)*state.rho
              + screened_rates(k_S35_to_He4_Si31_reaclib);
    jac.set(He4, S35, scratch);

    scratch = -screened_rates(k_Cl33_He4_to_K37_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Cl33_He4_to_p_Ar36_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Cl33_to_He4_P29_reaclib);
    jac.set(He4, Cl33, scratch);

    scratch = -screened_rates(k_Cl35_He4_to_K39_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Cl35_He4_to_p_Ar38_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Cl35_p_to_He4_S32_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Cl35_to_He4_P31_reaclib);
    jac.set(He4, Cl35, scratch);

    scratch = screened_rates(k_Cl37_p_to_He4_S34_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Cl37_to_He4_P33_reaclib);
    jac.set(He4, Cl37, scratch);

    scratch = -screened_rates(k_Ar36_He4_to_Ca40_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ar36_He4_to_p_K39_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Ar36_n_to_He4_S33_reaclib)*Y(N)*state.rho
              + screened_rates(k_Ar36_p_to_He4_Cl33_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Ar36_to_He4_S32_reaclib);
    jac.set(He4, Ar36, scratch);

    scratch = -screened_rates(k_Ar37_He4_to_n_Ca40_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Ar37_n_to_He4_S34_reaclib)*Y(N)*state.rho
              + screened_rates(k_Ar37_to_He4_S33_reaclib);
    jac.set(He4, Ar37, scratch);

    scratch = -screened_rates(k_Ar38_He4_to_Ca42_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Ar38_n_to_He4_S35_reaclib)*Y(N)*state.rho
              + screened_rates(k_Ar38_p_to_He4_Cl35_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Ar38_to_He4_S34_reaclib);
    jac.set(He4, Ar38, scratch);

    scratch = -screened_rates(k_K37_He4_to_Sc41_reaclib)*Y(He4)*state.rho
              - screened_rates(k_K37_He4_to_p_Ca40_reaclib)*Y(He4)*state.rho
              + screened_rates(k_K37_to_He4_Cl33_reaclib);
    jac.set(He4, K37, scratch);

    scratch = -screened_rates(k_K39_He4_to_Sc43_reaclib)*Y(He4)*state.rho
              - screened_rates(k_K39_He4_to_p_Ca42_reaclib)*Y(He4)*state.rho
              + screened_rates(k_K39_p_to_He4_Ar36_reaclib)*Y(H1)*state.rho
              + screened_rates(k_K39_to_He4_Cl35_reaclib);
    jac.set(He4, K39, scratch);

    scratch = -screened_rates(k_Ca40_He4_to_Ti44_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ca40_He4_to_p_Sc43_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Ca40_n_to_He4_Ar37_reaclib)*Y(N)*state.rho
              + screened_rates(k_Ca40_p_to_He4_K37_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Ca40_to_He4_Ar36_reaclib);
    jac.set(He4, Ca40, scratch);

    scratch = -screened_rates(k_Ca42_He4_to_Ti46_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ca42_He4_to_p_Sc45_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Ca42_p_to_He4_K39_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Ca42_to_He4_Ar38_reaclib);
    jac.set(He4, Ca42, scratch);

    scratch = -screened_rates(k_Sc41_He4_to_p_Ti44_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Sc41_to_He4_K37_reaclib);
    jac.set(He4, Sc41, scratch);

    scratch = -screened_rates(k_Sc43_He4_to_V47_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Sc43_He4_to_p_Ti46_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Sc43_p_to_He4_Ca40_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Sc43_to_He4_K39_reaclib);
    jac.set(He4, Sc43, scratch);

    scratch = -screened_rates(k_Sc45_He4_to_V49_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Sc45_He4_to_p_Ti48_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Sc45_p_to_He4_Ca42_reaclib)*Y(H1)*state.rho;
    jac.set(He4, Sc45, scratch);

    scratch = -screened_rates(k_Ti44_He4_to_Cr48_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ti44_He4_to_p_V47_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Ti44_p_to_He4_Sc41_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Ti44_to_He4_Ca40_reaclib);
    jac.set(He4, Ti44, scratch);

    scratch = -screened_rates(k_Ti46_He4_to_Cr50_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ti46_He4_to_p_V49_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Ti46_p_to_He4_Sc43_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Ti46_to_He4_Ca42_reaclib);
    jac.set(He4, Ti46, scratch);

    scratch = -screened_rates(k_Ti48_He4_to_Cr52_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ti48_He4_to_p_V51_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Ti48_p_to_He4_Sc45_reaclib)*Y(H1)*state.rho;
    jac.set(He4, Ti48, scratch);

    scratch = -screened_rates(k_V47_He4_to_Mn51_reaclib)*Y(He4)*state.rho
              - screened_rates(k_V47_He4_to_p_Cr50_reaclib)*Y(He4)*state.rho
              + screened_rates(k_V47_p_to_He4_Ti44_reaclib)*Y(H1)*state.rho
              + screened_rates(k_V47_to_He4_Sc43_reaclib);
    jac.set(He4, V47, scratch);

    scratch = -screened_rates(k_V49_He4_to_Mn53_reaclib)*Y(He4)*state.rho
              - screened_rates(k_V49_He4_to_p_Cr52_reaclib)*Y(He4)*state.rho
              + screened_rates(k_V49_p_to_He4_Ti46_reaclib)*Y(H1)*state.rho
              + screened_rates(k_V49_to_He4_Sc45_reaclib);
    jac.set(He4, V49, scratch);

    scratch = -screened_rates(k_V51_He4_to_Mn55_reaclib)*Y(He4)*state.rho
              - screened_rates(k_V51_He4_to_p_Cr54_reaclib)*Y(He4)*state.rho
              + screened_rates(k_V51_p_to_He4_Ti48_reaclib)*Y(H1)*state.rho;
    jac.set(He4, V51, scratch);

    scratch = -screened_rates(k_Cr48_He4_to_Fe52_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Cr48_He4_to_p_Mn51_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Cr48_to_He4_Ti44_reaclib);
    jac.set(He4, Cr48, scratch);

    scratch = -screened_rates(k_Cr50_He4_to_Fe54_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Cr50_He4_to_p_Mn53_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Cr50_p_to_He4_V47_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Cr50_to_He4_Ti46_reaclib);
    jac.set(He4, Cr50, scratch);

    scratch = -screened_rates(k_Cr52_He4_to_Fe56_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Cr52_He4_to_p_Mn55_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Cr52_p_to_He4_V49_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Cr52_to_He4_Ti48_reaclib);
    jac.set(He4, Cr52, scratch);

    scratch = -screened_rates(k_Cr54_He4_to_Fe58_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Cr54_p_to_He4_V51_reaclib)*Y(H1)*state.rho;
    jac.set(He4, Cr54, scratch);

    scratch = -screened_rates(k_Mn51_He4_to_Co55_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Mn51_He4_to_p_Fe54_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Mn51_p_to_He4_Cr48_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Mn51_to_He4_V47_reaclib);
    jac.set(He4, Mn51, scratch);

    scratch = -screened_rates(k_Mn53_He4_to_Co57_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Mn53_He4_to_n_Co56_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Mn53_He4_to_p_Fe56_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Mn53_p_to_He4_Cr50_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Mn53_to_He4_V49_reaclib);
    jac.set(He4, Mn53, scratch);

    scratch = -screened_rates(k_Mn55_He4_to_Co59_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Mn55_He4_to_p_Fe58_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Mn55_p_to_He4_Cr52_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Mn55_to_He4_V51_reaclib);
    jac.set(He4, Mn55, scratch);

    scratch = -screened_rates(k_Fe52_He4_to_Ni56_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Fe52_He4_to_p_Co55_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Fe52_to_He4_Cr48_reaclib);
    jac.set(He4, Fe52, scratch);

    scratch = -screened_rates(k_Fe54_He4_to_Ni58_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Fe54_He4_to_p_Co57_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Fe54_p_to_He4_Mn51_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Fe54_to_He4_Cr50_reaclib);
    jac.set(He4, Fe54, scratch);

    scratch = -screened_rates(k_Fe56_He4_to_Ni60_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Fe56_He4_to_p_Co59_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Fe56_p_to_He4_Mn53_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Fe56_to_He4_Cr52_reaclib);
    jac.set(He4, Fe56, scratch);

    scratch = screened_rates(k_Fe58_p_to_He4_Mn55_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Fe58_to_He4_Cr54_reaclib);
    jac.set(He4, Fe58, scratch);

    scratch = -screened_rates(k_Co55_He4_to_Cu59_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Co55_He4_to_p_Ni58_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Co55_p_to_He4_Fe52_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Co55_to_He4_Mn51_reaclib);
    jac.set(He4, Co55, scratch);

    scratch = -screened_rates(k_Co56_He4_to_n_Cu59_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Co56_n_to_He4_Mn53_reaclib)*Y(N)*state.rho;
    jac.set(He4, Co56, scratch);

    scratch = -screened_rates(k_Co57_He4_to_p_Ni60_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Co57_p_to_He4_Fe54_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Co57_to_He4_Mn53_reaclib);
    jac.set(He4, Co57, scratch);

    scratch = screened_rates(k_Co59_p_to_He4_Fe56_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Co59_to_He4_Mn55_reaclib);
    jac.set(He4, Co59, scratch);

    scratch = -screened_rates(k_Ni56_He4_to_Zn60_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ni56_He4_to_p_Cu59_reaclib)*Y(He4)*state.rho
              + screened_rates(k_Ni56_to_He4_Fe52_reaclib);
    jac.set(He4, Ni56, scratch);

    scratch = screened_rates(k_Ni58_p_to_He4_Co55_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Ni58_to_He4_Fe54_reaclib);
    jac.set(He4, Ni58, scratch);

    scratch = screened_rates(k_Ni60_p_to_He4_Co57_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Ni60_to_He4_Fe56_reaclib);
    jac.set(He4, Ni60, scratch);

    scratch = screened_rates(k_Cu59_n_to_He4_Co56_reaclib)*Y(N)*state.rho
              + screened_rates(k_Cu59_p_to_He4_Ni56_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Cu59_to_He4_Co55_reaclib);
    jac.set(He4, Cu59, scratch);

    scratch = screened_rates(k_Zn60_to_He4_Ni56_reaclib);
    jac.set(He4, Zn60, scratch);


    // C12 row

    scratch = 2.0*screened_rates(k_Mg23_n_to_C12_C12_reaclib)*Y(Mg23)*state.rho
              + screened_rates(k_Si27_n_to_C12_O16_reaclib)*Y(Si27)*state.rho;
    jac.set(C12, N, scratch);

    scratch = screened_rates(k_Al27_p_to_C12_O16_reaclib)*Y(Al27)*state.rho
              + screened_rates(k_P31_p_to_C12_Ne20_reaclib)*Y(P31)*state.rho;
    jac.set(C12, H1, scratch);

    scratch = -screened_rates(k_C12_He4_to_O16_reaclib)*Y(C12)*state.rho
              + 0.5*screened_rates(k_He4_He4_He4_to_C12_reaclib)*amrex::Math::powi<2>(Y(He4))*amrex::Math::powi<2>(state.rho)
              + screened_rates(k_Mg24_He4_to_C12_O16_reaclib)*Y(Mg24)*state.rho
              + 2.0*screened_rates(k_Ne20_He4_to_C12_C12_approx)*Y(Ne20)*state.rho
              + screened_rates(k_Si28_He4_to_C12_Ne20_reaclib)*Y(Si28)*state.rho;
    jac.set(C12, He4, scratch);

    scratch = -2.0*screened_rates(k_C12_C12_to_Mg24_approx)*Y(C12)*state.rho
              - 2.0*screened_rates(k_C12_C12_to_Ne20_He4_approx)*Y(C12)*state.rho
              - 2.0*screened_rates(k_C12_C12_to_n_Mg23_reaclib)*Y(C12)*state.rho
              - screened_rates(k_C12_He4_to_O16_reaclib)*Y(He4)*state.rho
              - screened_rates(k_C12_to_3He4_reaclib)
              - screened_rates(k_Ne20_C12_to_He4_Si28_reaclib)*Y(Ne20)*state.rho
              - screened_rates(k_Ne20_C12_to_p_P31_reaclib)*Y(Ne20)*state.rho
              - screened_rates(k_O16_C12_to_He4_Mg24_reaclib)*Y(O16)*state.rho
              - screened_rates(k_O16_C12_to_n_Si27_reaclib)*Y(O16)*state.rho
              - screened_rates(k_O16_C12_to_p_Al27_reaclib)*Y(O16)*state.rho;
    jac.set(C12, C12, scratch);

    scratch = -screened_rates(k_O16_C12_to_He4_Mg24_reaclib)*Y(C12)*state.rho
              - screened_rates(k_O16_C12_to_n_Si27_reaclib)*Y(C12)*state.rho
              - screened_rates(k_O16_C12_to_p_Al27_reaclib)*Y(C12)*state.rho
              + screened_rates(k_O16_to_He4_C12_reaclib);
    jac.set(C12, O16, scratch);

    scratch = -screened_rates(k_Ne20_C12_to_He4_Si28_reaclib)*Y(C12)*state.rho
              - screened_rates(k_Ne20_C12_to_p_P31_reaclib)*Y(C12)*state.rho
              + 2.0*screened_rates(k_Ne20_He4_to_C12_C12_approx)*Y(He4)*state.rho;
    jac.set(C12, Ne20, scratch);

    scratch = 2.0*screened_rates(k_Mg23_n_to_C12_C12_reaclib)*Y(N)*state.rho;
    jac.set(C12, Mg23, scratch);

    scratch = screened_rates(k_Mg24_He4_to_C12_O16_reaclib)*Y(He4)*state.rho
              + 2.0*screened_rates(k_Mg24_to_C12_C12_approx);
    jac.set(C12, Mg24, scratch);

    scratch = screened_rates(k_Al27_p_to_C12_O16_reaclib)*Y(H1)*state.rho;
    jac.set(C12, Al27, scratch);

    scratch = screened_rates(k_Si27_n_to_C12_O16_reaclib)*Y(N)*state.rho;
    jac.set(C12, Si27, scratch);

    scratch = screened_rates(k_Si28_He4_to_C12_Ne20_reaclib)*Y(He4)*state.rho;
    jac.set(C12, Si28, scratch);

    jac.set(C12, Si29, 0.0);

    jac.set(C12, Si30, 0.0);

    jac.set(C12, Si31, 0.0);

    jac.set(C12, P29, 0.0);

    scratch = screened_rates(k_P31_p_to_C12_Ne20_reaclib)*Y(H1)*state.rho;
    jac.set(C12, P31, scratch);

    jac.set(C12, P33, 0.0);

    jac.set(C12, S32, 0.0);

    jac.set(C12, S33, 0.0);

    jac.set(C12, S34, 0.0);

    jac.set(C12, S35, 0.0);

    jac.set(C12, Cl33, 0.0);

    jac.set(C12, Cl35, 0.0);

    jac.set(C12, Cl37, 0.0);

    jac.set(C12, Ar36, 0.0);

    jac.set(C12, Ar37, 0.0);

    jac.set(C12, Ar38, 0.0);

    jac.set(C12, K37, 0.0);

    jac.set(C12, K39, 0.0);

    jac.set(C12, Ca40, 0.0);

    jac.set(C12, Ca42, 0.0);

    jac.set(C12, Sc41, 0.0);

    jac.set(C12, Sc43, 0.0);

    jac.set(C12, Sc45, 0.0);

    jac.set(C12, Ti44, 0.0);

    jac.set(C12, Ti46, 0.0);

    jac.set(C12, Ti48, 0.0);

    jac.set(C12, V47, 0.0);

    jac.set(C12, V49, 0.0);

    jac.set(C12, V51, 0.0);

    jac.set(C12, Cr48, 0.0);

    jac.set(C12, Cr50, 0.0);

    jac.set(C12, Cr52, 0.0);

    jac.set(C12, Cr54, 0.0);

    jac.set(C12, Mn51, 0.0);

    jac.set(C12, Mn53, 0.0);

    jac.set(C12, Mn55, 0.0);

    jac.set(C12, Fe52, 0.0);

    jac.set(C12, Fe54, 0.0);

    jac.set(C12, Fe56, 0.0);

    jac.set(C12, Fe58, 0.0);

    jac.set(C12, Co55, 0.0);

    jac.set(C12, Co56, 0.0);

    jac.set(C12, Co57, 0.0);

    jac.set(C12, Co59, 0.0);

    jac.set(C12, Ni56, 0.0);

    jac.set(C12, Ni58, 0.0);

    jac.set(C12, Ni60, 0.0);

    jac.set(C12, Cu59, 0.0);

    jac.set(C12, Zn60, 0.0);


    // O16 row

    scratch = screened_rates(k_Si27_n_to_C12_O16_reaclib)*Y(Si27)*state.rho;
    jac.set(O16, N, scratch);

    scratch = screened_rates(k_Al27_p_to_C12_O16_reaclib)*Y(Al27)*state.rho
              + 2.0*screened_rates(k_P31_p_to_O16_O16_reaclib)*Y(P31)*state.rho;
    jac.set(O16, H1, scratch);

    scratch = screened_rates(k_C12_He4_to_O16_reaclib)*Y(C12)*state.rho
              + screened_rates(k_Mg24_He4_to_C12_O16_reaclib)*Y(Mg24)*state.rho
              - screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(O16)*state.rho
              + 2.0*screened_rates(k_Si28_He4_to_O16_O16_reaclib)*Y(Si28)*state.rho;
    jac.set(O16, He4, scratch);

    scratch = screened_rates(k_C12_He4_to_O16_reaclib)*Y(He4)*state.rho
              - screened_rates(k_O16_C12_to_He4_Mg24_reaclib)*Y(O16)*state.rho
              - screened_rates(k_O16_C12_to_n_Si27_reaclib)*Y(O16)*state.rho
              - screened_rates(k_O16_C12_to_p_Al27_reaclib)*Y(O16)*state.rho;
    jac.set(O16, C12, scratch);

    scratch = -screened_rates(k_O16_C12_to_He4_Mg24_reaclib)*Y(C12)*state.rho
              - screened_rates(k_O16_C12_to_n_Si27_reaclib)*Y(C12)*state.rho
              - screened_rates(k_O16_C12_to_p_Al27_reaclib)*Y(C12)*state.rho
              - screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(He4)*state.rho
              - 2.0*screened_rates(k_O16_O16_to_He4_Si28_reaclib)*Y(O16)*state.rho
              - 2.0*screened_rates(k_O16_O16_to_p_P31_reaclib)*Y(O16)*state.rho
              - screened_rates(k_O16_to_He4_C12_reaclib);
    jac.set(O16, O16, scratch);

    scratch = screened_rates(k_Ne20_to_He4_O16_reaclib);
    jac.set(O16, Ne20, scratch);

    jac.set(O16, Mg23, 0.0);

    scratch = screened_rates(k_Mg24_He4_to_C12_O16_reaclib)*Y(He4)*state.rho;
    jac.set(O16, Mg24, scratch);

    scratch = screened_rates(k_Al27_p_to_C12_O16_reaclib)*Y(H1)*state.rho;
    jac.set(O16, Al27, scratch);

    scratch = screened_rates(k_Si27_n_to_C12_O16_reaclib)*Y(N)*state.rho;
    jac.set(O16, Si27, scratch);

    scratch = 2.0*screened_rates(k_Si28_He4_to_O16_O16_reaclib)*Y(He4)*state.rho;
    jac.set(O16, Si28, scratch);

    jac.set(O16, Si29, 0.0);

    jac.set(O16, Si30, 0.0);

    jac.set(O16, Si31, 0.0);

    jac.set(O16, P29, 0.0);

    scratch = 2.0*screened_rates(k_P31_p_to_O16_O16_reaclib)*Y(H1)*state.rho;
    jac.set(O16, P31, scratch);

    jac.set(O16, P33, 0.0);

    jac.set(O16, S32, 0.0);

    jac.set(O16, S33, 0.0);

    jac.set(O16, S34, 0.0);

    jac.set(O16, S35, 0.0);

    jac.set(O16, Cl33, 0.0);

    jac.set(O16, Cl35, 0.0);

    jac.set(O16, Cl37, 0.0);

    jac.set(O16, Ar36, 0.0);

    jac.set(O16, Ar37, 0.0);

    jac.set(O16, Ar38, 0.0);

    jac.set(O16, K37, 0.0);

    jac.set(O16, K39, 0.0);

    jac.set(O16, Ca40, 0.0);

    jac.set(O16, Ca42, 0.0);

    jac.set(O16, Sc41, 0.0);

    jac.set(O16, Sc43, 0.0);

    jac.set(O16, Sc45, 0.0);

    jac.set(O16, Ti44, 0.0);

    jac.set(O16, Ti46, 0.0);

    jac.set(O16, Ti48, 0.0);

    jac.set(O16, V47, 0.0);

    jac.set(O16, V49, 0.0);

    jac.set(O16, V51, 0.0);

    jac.set(O16, Cr48, 0.0);

    jac.set(O16, Cr50, 0.0);

    jac.set(O16, Cr52, 0.0);

    jac.set(O16, Cr54, 0.0);

    jac.set(O16, Mn51, 0.0);

    jac.set(O16, Mn53, 0.0);

    jac.set(O16, Mn55, 0.0);

    jac.set(O16, Fe52, 0.0);

    jac.set(O16, Fe54, 0.0);

    jac.set(O16, Fe56, 0.0);

    jac.set(O16, Fe58, 0.0);

    jac.set(O16, Co55, 0.0);

    jac.set(O16, Co56, 0.0);

    jac.set(O16, Co57, 0.0);

    jac.set(O16, Co59, 0.0);

    jac.set(O16, Ni56, 0.0);

    jac.set(O16, Ni58, 0.0);

    jac.set(O16, Ni60, 0.0);

    jac.set(O16, Cu59, 0.0);

    jac.set(O16, Zn60, 0.0);


    // Ne20 row

    scratch = screened_rates(k_Mg23_n_to_He4_Ne20_reaclib)*Y(Mg23)*state.rho;
    jac.set(Ne20, N, scratch);

    scratch = screened_rates(k_P31_p_to_C12_Ne20_reaclib)*Y(P31)*state.rho;
    jac.set(Ne20, H1, scratch);

    scratch = -screened_rates(k_Ne20_He4_to_C12_C12_approx)*Y(Ne20)*state.rho
              - screened_rates(k_Ne20_He4_to_Mg24_approx)*Y(Ne20)*state.rho
              - screened_rates(k_Ne20_He4_to_n_Mg23_reaclib)*Y(Ne20)*state.rho
              + screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(O16)*state.rho
              + screened_rates(k_Si28_He4_to_C12_Ne20_reaclib)*Y(Si28)*state.rho;
    jac.set(Ne20, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_Ne20_He4_approx)*Y(C12)*state.rho
              - screened_rates(k_Ne20_C12_to_He4_Si28_reaclib)*Y(Ne20)*state.rho
              - screened_rates(k_Ne20_C12_to_p_P31_reaclib)*Y(Ne20)*state.rho;
    jac.set(Ne20, C12, scratch);

    scratch = screened_rates(k_O16_He4_to_Ne20_reaclib)*Y(He4)*state.rho;
    jac.set(Ne20, O16, scratch);

    scratch = -screened_rates(k_Ne20_C12_to_He4_Si28_reaclib)*Y(C12)*state.rho
              - screened_rates(k_Ne20_C12_to_p_P31_reaclib)*Y(C12)*state.rho
              - screened_rates(k_Ne20_He4_to_C12_C12_approx)*Y(He4)*state.rho
              - screened_rates(k_Ne20_He4_to_Mg24_approx)*Y(He4)*state.rho
              - screened_rates(k_Ne20_He4_to_n_Mg23_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ne20_to_He4_O16_reaclib);
    jac.set(Ne20, Ne20, scratch);

    scratch = screened_rates(k_Mg23_n_to_He4_Ne20_reaclib)*Y(N)*state.rho;
    jac.set(Ne20, Mg23, scratch);

    scratch = screened_rates(k_Mg24_to_Ne20_He4_approx);
    jac.set(Ne20, Mg24, scratch);

    jac.set(Ne20, Al27, 0.0);

    jac.set(Ne20, Si27, 0.0);

    scratch = screened_rates(k_Si28_He4_to_C12_Ne20_reaclib)*Y(He4)*state.rho;
    jac.set(Ne20, Si28, scratch);

    jac.set(Ne20, Si29, 0.0);

    jac.set(Ne20, Si30, 0.0);

    jac.set(Ne20, Si31, 0.0);

    jac.set(Ne20, P29, 0.0);

    scratch = screened_rates(k_P31_p_to_C12_Ne20_reaclib)*Y(H1)*state.rho;
    jac.set(Ne20, P31, scratch);

    jac.set(Ne20, P33, 0.0);

    jac.set(Ne20, S32, 0.0);

    jac.set(Ne20, S33, 0.0);

    jac.set(Ne20, S34, 0.0);

    jac.set(Ne20, S35, 0.0);

    jac.set(Ne20, Cl33, 0.0);

    jac.set(Ne20, Cl35, 0.0);

    jac.set(Ne20, Cl37, 0.0);

    jac.set(Ne20, Ar36, 0.0);

    jac.set(Ne20, Ar37, 0.0);

    jac.set(Ne20, Ar38, 0.0);

    jac.set(Ne20, K37, 0.0);

    jac.set(Ne20, K39, 0.0);

    jac.set(Ne20, Ca40, 0.0);

    jac.set(Ne20, Ca42, 0.0);

    jac.set(Ne20, Sc41, 0.0);

    jac.set(Ne20, Sc43, 0.0);

    jac.set(Ne20, Sc45, 0.0);

    jac.set(Ne20, Ti44, 0.0);

    jac.set(Ne20, Ti46, 0.0);

    jac.set(Ne20, Ti48, 0.0);

    jac.set(Ne20, V47, 0.0);

    jac.set(Ne20, V49, 0.0);

    jac.set(Ne20, V51, 0.0);

    jac.set(Ne20, Cr48, 0.0);

    jac.set(Ne20, Cr50, 0.0);

    jac.set(Ne20, Cr52, 0.0);

    jac.set(Ne20, Cr54, 0.0);

    jac.set(Ne20, Mn51, 0.0);

    jac.set(Ne20, Mn53, 0.0);

    jac.set(Ne20, Mn55, 0.0);

    jac.set(Ne20, Fe52, 0.0);

    jac.set(Ne20, Fe54, 0.0);

    jac.set(Ne20, Fe56, 0.0);

    jac.set(Ne20, Fe58, 0.0);

    jac.set(Ne20, Co55, 0.0);

    jac.set(Ne20, Co56, 0.0);

    jac.set(Ne20, Co57, 0.0);

    jac.set(Ne20, Co59, 0.0);

    jac.set(Ne20, Ni56, 0.0);

    jac.set(Ne20, Ni58, 0.0);

    jac.set(Ne20, Ni60, 0.0);

    jac.set(Ne20, Cu59, 0.0);

    jac.set(Ne20, Zn60, 0.0);


    // Mg23 row

    scratch = -screened_rates(k_Mg23_n_to_C12_C12_reaclib)*Y(Mg23)*state.rho
              - screened_rates(k_Mg23_n_to_He4_Ne20_reaclib)*Y(Mg23)*state.rho
              - screened_rates(k_Mg23_n_to_Mg24_reaclib)*Y(Mg23)*state.rho;
    jac.set(Mg23, N, scratch);

    jac.set(Mg23, H1, 0.0);

    scratch = -screened_rates(k_Mg23_He4_to_Si27_reaclib)*Y(Mg23)*state.rho
              + screened_rates(k_Ne20_He4_to_n_Mg23_reaclib)*Y(Ne20)*state.rho;
    jac.set(Mg23, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_n_Mg23_reaclib)*Y(C12)*state.rho;
    jac.set(Mg23, C12, scratch);

    jac.set(Mg23, O16, 0.0);

    scratch = screened_rates(k_Ne20_He4_to_n_Mg23_reaclib)*Y(He4)*state.rho;
    jac.set(Mg23, Ne20, scratch);

    scratch = -screened_rates(k_Mg23_He4_to_Si27_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Mg23_n_to_C12_C12_reaclib)*Y(N)*state.rho
              - screened_rates(k_Mg23_n_to_He4_Ne20_reaclib)*Y(N)*state.rho
              - screened_rates(k_Mg23_n_to_Mg24_reaclib)*Y(N)*state.rho;
    jac.set(Mg23, Mg23, scratch);

    scratch = screened_rates(k_Mg24_to_n_Mg23_reaclib);
    jac.set(Mg23, Mg24, scratch);

    jac.set(Mg23, Al27, 0.0);

    scratch = screened_rates(k_Si27_to_He4_Mg23_reaclib);
    jac.set(Mg23, Si27, scratch);

    jac.set(Mg23, Si28, 0.0);

    jac.set(Mg23, Si29, 0.0);

    jac.set(Mg23, Si30, 0.0);

    jac.set(Mg23, Si31, 0.0);

    jac.set(Mg23, P29, 0.0);

    jac.set(Mg23, P31, 0.0);

    jac.set(Mg23, P33, 0.0);

    jac.set(Mg23, S32, 0.0);

    jac.set(Mg23, S33, 0.0);

    jac.set(Mg23, S34, 0.0);

    jac.set(Mg23, S35, 0.0);

    jac.set(Mg23, Cl33, 0.0);

    jac.set(Mg23, Cl35, 0.0);

    jac.set(Mg23, Cl37, 0.0);

    jac.set(Mg23, Ar36, 0.0);

    jac.set(Mg23, Ar37, 0.0);

    jac.set(Mg23, Ar38, 0.0);

    jac.set(Mg23, K37, 0.0);

    jac.set(Mg23, K39, 0.0);

    jac.set(Mg23, Ca40, 0.0);

    jac.set(Mg23, Ca42, 0.0);

    jac.set(Mg23, Sc41, 0.0);

    jac.set(Mg23, Sc43, 0.0);

    jac.set(Mg23, Sc45, 0.0);

    jac.set(Mg23, Ti44, 0.0);

    jac.set(Mg23, Ti46, 0.0);

    jac.set(Mg23, Ti48, 0.0);

    jac.set(Mg23, V47, 0.0);

    jac.set(Mg23, V49, 0.0);

    jac.set(Mg23, V51, 0.0);

    jac.set(Mg23, Cr48, 0.0);

    jac.set(Mg23, Cr50, 0.0);

    jac.set(Mg23, Cr52, 0.0);

    jac.set(Mg23, Cr54, 0.0);

    jac.set(Mg23, Mn51, 0.0);

    jac.set(Mg23, Mn53, 0.0);

    jac.set(Mg23, Mn55, 0.0);

    jac.set(Mg23, Fe52, 0.0);

    jac.set(Mg23, Fe54, 0.0);

    jac.set(Mg23, Fe56, 0.0);

    jac.set(Mg23, Fe58, 0.0);

    jac.set(Mg23, Co55, 0.0);

    jac.set(Mg23, Co56, 0.0);

    jac.set(Mg23, Co57, 0.0);

    jac.set(Mg23, Co59, 0.0);

    jac.set(Mg23, Ni56, 0.0);

    jac.set(Mg23, Ni58, 0.0);

    jac.set(Mg23, Ni60, 0.0);

    jac.set(Mg23, Cu59, 0.0);

    jac.set(Mg23, Zn60, 0.0);


    // Mg24 row

    scratch = screened_rates(k_Mg23_n_to_Mg24_reaclib)*Y(Mg23)*state.rho
              + screened_rates(k_Si27_n_to_He4_Mg24_reaclib)*Y(Si27)*state.rho;
    jac.set(Mg24, N, scratch);

    scratch = screened_rates(k_Al27_p_to_He4_Mg24_reaclib)*Y(Al27)*state.rho;
    jac.set(Mg24, H1, scratch);

    scratch = -screened_rates(k_Mg24_He4_to_C12_O16_reaclib)*Y(Mg24)*state.rho
              - screened_rates(k_Mg24_He4_to_Si28_reaclib)*Y(Mg24)*state.rho
              - screened_rates(k_Mg24_He4_to_n_Si27_reaclib)*Y(Mg24)*state.rho
              - screened_rates(k_Mg24_He4_to_p_Al27_reaclib)*Y(Mg24)*state.rho
              + screened_rates(k_Ne20_He4_to_Mg24_approx)*Y(Ne20)*state.rho;
    jac.set(Mg24, He4, scratch);

    scratch = 1.0*screened_rates(k_C12_C12_to_Mg24_approx)*Y(C12)*state.rho
              + screened_rates(k_O16_C12_to_He4_Mg24_reaclib)*Y(O16)*state.rho;
    jac.set(Mg24, C12, scratch);

    scratch = screened_rates(k_O16_C12_to_He4_Mg24_reaclib)*Y(C12)*state.rho;
    jac.set(Mg24, O16, scratch);

    scratch = screened_rates(k_Ne20_He4_to_Mg24_approx)*Y(He4)*state.rho;
    jac.set(Mg24, Ne20, scratch);

    scratch = screened_rates(k_Mg23_n_to_Mg24_reaclib)*Y(N)*state.rho;
    jac.set(Mg24, Mg23, scratch);

    scratch = -screened_rates(k_Mg24_He4_to_C12_O16_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Mg24_He4_to_Si28_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Mg24_He4_to_n_Si27_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Mg24_He4_to_p_Al27_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Mg24_to_C12_C12_approx)
              - screened_rates(k_Mg24_to_Ne20_He4_approx)
              - screened_rates(k_Mg24_to_n_Mg23_reaclib);
    jac.set(Mg24, Mg24, scratch);

    scratch = screened_rates(k_Al27_p_to_He4_Mg24_reaclib)*Y(H1)*state.rho;
    jac.set(Mg24, Al27, scratch);

    scratch = screened_rates(k_Si27_n_to_He4_Mg24_reaclib)*Y(N)*state.rho;
    jac.set(Mg24, Si27, scratch);

    scratch = screened_rates(k_Si28_to_He4_Mg24_reaclib);
    jac.set(Mg24, Si28, scratch);

    jac.set(Mg24, Si29, 0.0);

    jac.set(Mg24, Si30, 0.0);

    jac.set(Mg24, Si31, 0.0);

    jac.set(Mg24, P29, 0.0);

    jac.set(Mg24, P31, 0.0);

    jac.set(Mg24, P33, 0.0);

    jac.set(Mg24, S32, 0.0);

    jac.set(Mg24, S33, 0.0);

    jac.set(Mg24, S34, 0.0);

    jac.set(Mg24, S35, 0.0);

    jac.set(Mg24, Cl33, 0.0);

    jac.set(Mg24, Cl35, 0.0);

    jac.set(Mg24, Cl37, 0.0);

    jac.set(Mg24, Ar36, 0.0);

    jac.set(Mg24, Ar37, 0.0);

    jac.set(Mg24, Ar38, 0.0);

    jac.set(Mg24, K37, 0.0);

    jac.set(Mg24, K39, 0.0);

    jac.set(Mg24, Ca40, 0.0);

    jac.set(Mg24, Ca42, 0.0);

    jac.set(Mg24, Sc41, 0.0);

    jac.set(Mg24, Sc43, 0.0);

    jac.set(Mg24, Sc45, 0.0);

    jac.set(Mg24, Ti44, 0.0);

    jac.set(Mg24, Ti46, 0.0);

    jac.set(Mg24, Ti48, 0.0);

    jac.set(Mg24, V47, 0.0);

    jac.set(Mg24, V49, 0.0);

    jac.set(Mg24, V51, 0.0);

    jac.set(Mg24, Cr48, 0.0);

    jac.set(Mg24, Cr50, 0.0);

    jac.set(Mg24, Cr52, 0.0);

    jac.set(Mg24, Cr54, 0.0);

    jac.set(Mg24, Mn51, 0.0);

    jac.set(Mg24, Mn53, 0.0);

    jac.set(Mg24, Mn55, 0.0);

    jac.set(Mg24, Fe52, 0.0);

    jac.set(Mg24, Fe54, 0.0);

    jac.set(Mg24, Fe56, 0.0);

    jac.set(Mg24, Fe58, 0.0);

    jac.set(Mg24, Co55, 0.0);

    jac.set(Mg24, Co56, 0.0);

    jac.set(Mg24, Co57, 0.0);

    jac.set(Mg24, Co59, 0.0);

    jac.set(Mg24, Ni56, 0.0);

    jac.set(Mg24, Ni58, 0.0);

    jac.set(Mg24, Ni60, 0.0);

    jac.set(Mg24, Cu59, 0.0);

    jac.set(Mg24, Zn60, 0.0);


    // Al27 row

    scratch = screened_rates(k_Si27_n_to_p_Al27_reaclib)*Y(Si27)*state.rho;
    jac.set(Al27, N, scratch);

    scratch = -screened_rates(k_Al27_p_to_C12_O16_reaclib)*Y(Al27)*state.rho
              - screened_rates(k_Al27_p_to_He4_Mg24_reaclib)*Y(Al27)*state.rho
              - screened_rates(k_Al27_p_to_Si28_reaclib)*Y(Al27)*state.rho
              - screened_rates(k_Al27_p_to_n_Si27_reaclib)*Y(Al27)*state.rho
              + screened_rates(k_Si30_p_to_He4_Al27_reaclib)*Y(Si30)*state.rho;
    jac.set(Al27, H1, scratch);

    scratch = -screened_rates(k_Al27_He4_to_P31_reaclib)*Y(Al27)*state.rho
              - screened_rates(k_Al27_He4_to_p_Si30_reaclib)*Y(Al27)*state.rho
              + screened_rates(k_Mg24_He4_to_p_Al27_reaclib)*Y(Mg24)*state.rho;
    jac.set(Al27, He4, scratch);

    scratch = screened_rates(k_O16_C12_to_p_Al27_reaclib)*Y(O16)*state.rho;
    jac.set(Al27, C12, scratch);

    scratch = screened_rates(k_O16_C12_to_p_Al27_reaclib)*Y(C12)*state.rho;
    jac.set(Al27, O16, scratch);

    jac.set(Al27, Ne20, 0.0);

    jac.set(Al27, Mg23, 0.0);

    scratch = screened_rates(k_Mg24_He4_to_p_Al27_reaclib)*Y(He4)*state.rho;
    jac.set(Al27, Mg24, scratch);

    scratch = -screened_rates(k_Al27_He4_to_P31_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Al27_He4_to_p_Si30_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Al27_p_to_C12_O16_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Al27_p_to_He4_Mg24_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Al27_p_to_Si28_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Al27_p_to_n_Si27_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Al27_to_Si27_beta_neg_weaktab);
    jac.set(Al27, Al27, scratch);

    scratch = screened_rates(k_Si27_n_to_p_Al27_reaclib)*Y(N)*state.rho
              + screened_rates(k_Si27_to_Al27_electron_capture_weaktab);
    jac.set(Al27, Si27, scratch);

    scratch = screened_rates(k_Si28_to_p_Al27_reaclib);
    jac.set(Al27, Si28, scratch);

    jac.set(Al27, Si29, 0.0);

    scratch = screened_rates(k_Si30_p_to_He4_Al27_reaclib)*Y(H1)*state.rho;
    jac.set(Al27, Si30, scratch);

    jac.set(Al27, Si31, 0.0);

    jac.set(Al27, P29, 0.0);

    scratch = screened_rates(k_P31_to_He4_Al27_reaclib);
    jac.set(Al27, P31, scratch);

    jac.set(Al27, P33, 0.0);

    jac.set(Al27, S32, 0.0);

    jac.set(Al27, S33, 0.0);

    jac.set(Al27, S34, 0.0);

    jac.set(Al27, S35, 0.0);

    jac.set(Al27, Cl33, 0.0);

    jac.set(Al27, Cl35, 0.0);

    jac.set(Al27, Cl37, 0.0);

    jac.set(Al27, Ar36, 0.0);

    jac.set(Al27, Ar37, 0.0);

    jac.set(Al27, Ar38, 0.0);

    jac.set(Al27, K37, 0.0);

    jac.set(Al27, K39, 0.0);

    jac.set(Al27, Ca40, 0.0);

    jac.set(Al27, Ca42, 0.0);

    jac.set(Al27, Sc41, 0.0);

    jac.set(Al27, Sc43, 0.0);

    jac.set(Al27, Sc45, 0.0);

    jac.set(Al27, Ti44, 0.0);

    jac.set(Al27, Ti46, 0.0);

    jac.set(Al27, Ti48, 0.0);

    jac.set(Al27, V47, 0.0);

    jac.set(Al27, V49, 0.0);

    jac.set(Al27, V51, 0.0);

    jac.set(Al27, Cr48, 0.0);

    jac.set(Al27, Cr50, 0.0);

    jac.set(Al27, Cr52, 0.0);

    jac.set(Al27, Cr54, 0.0);

    jac.set(Al27, Mn51, 0.0);

    jac.set(Al27, Mn53, 0.0);

    jac.set(Al27, Mn55, 0.0);

    jac.set(Al27, Fe52, 0.0);

    jac.set(Al27, Fe54, 0.0);

    jac.set(Al27, Fe56, 0.0);

    jac.set(Al27, Fe58, 0.0);

    jac.set(Al27, Co55, 0.0);

    jac.set(Al27, Co56, 0.0);

    jac.set(Al27, Co57, 0.0);

    jac.set(Al27, Co59, 0.0);

    jac.set(Al27, Ni56, 0.0);

    jac.set(Al27, Ni58, 0.0);

    jac.set(Al27, Ni60, 0.0);

    jac.set(Al27, Cu59, 0.0);

    jac.set(Al27, Zn60, 0.0);


    // Si27 row

    scratch = -screened_rates(k_Si27_n_to_C12_O16_reaclib)*Y(Si27)*state.rho
              - screened_rates(k_Si27_n_to_He4_Mg24_reaclib)*Y(Si27)*state.rho
              - screened_rates(k_Si27_n_to_Si28_reaclib)*Y(Si27)*state.rho
              - screened_rates(k_Si27_n_to_p_Al27_reaclib)*Y(Si27)*state.rho;
    jac.set(Si27, N, scratch);

    scratch = screened_rates(k_Al27_p_to_n_Si27_reaclib)*Y(Al27)*state.rho;
    jac.set(Si27, H1, scratch);

    scratch = screened_rates(k_Mg23_He4_to_Si27_reaclib)*Y(Mg23)*state.rho
              + screened_rates(k_Mg24_He4_to_n_Si27_reaclib)*Y(Mg24)*state.rho;
    jac.set(Si27, He4, scratch);

    scratch = screened_rates(k_O16_C12_to_n_Si27_reaclib)*Y(O16)*state.rho;
    jac.set(Si27, C12, scratch);

    scratch = screened_rates(k_O16_C12_to_n_Si27_reaclib)*Y(C12)*state.rho;
    jac.set(Si27, O16, scratch);

    jac.set(Si27, Ne20, 0.0);

    scratch = screened_rates(k_Mg23_He4_to_Si27_reaclib)*Y(He4)*state.rho;
    jac.set(Si27, Mg23, scratch);

    scratch = screened_rates(k_Mg24_He4_to_n_Si27_reaclib)*Y(He4)*state.rho;
    jac.set(Si27, Mg24, scratch);

    scratch = screened_rates(k_Al27_p_to_n_Si27_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Al27_to_Si27_beta_neg_weaktab);
    jac.set(Si27, Al27, scratch);

    scratch = -screened_rates(k_Si27_n_to_C12_O16_reaclib)*Y(N)*state.rho
              - screened_rates(k_Si27_n_to_He4_Mg24_reaclib)*Y(N)*state.rho
              - screened_rates(k_Si27_n_to_Si28_reaclib)*Y(N)*state.rho
              - screened_rates(k_Si27_n_to_p_Al27_reaclib)*Y(N)*state.rho
              - screened_rates(k_Si27_to_Al27_electron_capture_weaktab)
              - screened_rates(k_Si27_to_He4_Mg23_reaclib);
    jac.set(Si27, Si27, scratch);

    scratch = screened_rates(k_Si28_to_n_Si27_reaclib);
    jac.set(Si27, Si28, scratch);

    jac.set(Si27, Si29, 0.0);

    jac.set(Si27, Si30, 0.0);

    jac.set(Si27, Si31, 0.0);

    jac.set(Si27, P29, 0.0);

    jac.set(Si27, P31, 0.0);

    jac.set(Si27, P33, 0.0);

    jac.set(Si27, S32, 0.0);

    jac.set(Si27, S33, 0.0);

    jac.set(Si27, S34, 0.0);

    jac.set(Si27, S35, 0.0);

    jac.set(Si27, Cl33, 0.0);

    jac.set(Si27, Cl35, 0.0);

    jac.set(Si27, Cl37, 0.0);

    jac.set(Si27, Ar36, 0.0);

    jac.set(Si27, Ar37, 0.0);

    jac.set(Si27, Ar38, 0.0);

    jac.set(Si27, K37, 0.0);

    jac.set(Si27, K39, 0.0);

    jac.set(Si27, Ca40, 0.0);

    jac.set(Si27, Ca42, 0.0);

    jac.set(Si27, Sc41, 0.0);

    jac.set(Si27, Sc43, 0.0);

    jac.set(Si27, Sc45, 0.0);

    jac.set(Si27, Ti44, 0.0);

    jac.set(Si27, Ti46, 0.0);

    jac.set(Si27, Ti48, 0.0);

    jac.set(Si27, V47, 0.0);

    jac.set(Si27, V49, 0.0);

    jac.set(Si27, V51, 0.0);

    jac.set(Si27, Cr48, 0.0);

    jac.set(Si27, Cr50, 0.0);

    jac.set(Si27, Cr52, 0.0);

    jac.set(Si27, Cr54, 0.0);

    jac.set(Si27, Mn51, 0.0);

    jac.set(Si27, Mn53, 0.0);

    jac.set(Si27, Mn55, 0.0);

    jac.set(Si27, Fe52, 0.0);

    jac.set(Si27, Fe54, 0.0);

    jac.set(Si27, Fe56, 0.0);

    jac.set(Si27, Fe58, 0.0);

    jac.set(Si27, Co55, 0.0);

    jac.set(Si27, Co56, 0.0);

    jac.set(Si27, Co57, 0.0);

    jac.set(Si27, Co59, 0.0);

    jac.set(Si27, Ni56, 0.0);

    jac.set(Si27, Ni58, 0.0);

    jac.set(Si27, Ni60, 0.0);

    jac.set(Si27, Cu59, 0.0);

    jac.set(Si27, Zn60, 0.0);


    // Si28 row

    scratch = screened_rates(k_Si27_n_to_Si28_reaclib)*Y(Si27)*state.rho
              - screened_rates(k_Si28_n_to_Si29_reaclib)*Y(Si28)*state.rho;
    jac.set(Si28, N, scratch);

    scratch = screened_rates(k_Al27_p_to_Si28_reaclib)*Y(Al27)*state.rho
              + screened_rates(k_P31_p_to_He4_Si28_reaclib)*Y(P31)*state.rho
              - screened_rates(k_Si28_p_to_P29_reaclib)*Y(Si28)*state.rho;
    jac.set(Si28, H1, scratch);

    scratch = screened_rates(k_Mg24_He4_to_Si28_reaclib)*Y(Mg24)*state.rho
              - screened_rates(k_Si28_He4_to_C12_Ne20_reaclib)*Y(Si28)*state.rho
              - screened_rates(k_Si28_He4_to_O16_O16_reaclib)*Y(Si28)*state.rho
              - screened_rates(k_Si28_He4_to_S32_reaclib)*Y(Si28)*state.rho
              - screened_rates(k_Si28_He4_to_p_P31_reaclib)*Y(Si28)*state.rho;
    jac.set(Si28, He4, scratch);

    scratch = screened_rates(k_Ne20_C12_to_He4_Si28_reaclib)*Y(Ne20)*state.rho;
    jac.set(Si28, C12, scratch);

    scratch = 1.0*screened_rates(k_O16_O16_to_He4_Si28_reaclib)*Y(O16)*state.rho;
    jac.set(Si28, O16, scratch);

    scratch = screened_rates(k_Ne20_C12_to_He4_Si28_reaclib)*Y(C12)*state.rho;
    jac.set(Si28, Ne20, scratch);

    jac.set(Si28, Mg23, 0.0);

    scratch = screened_rates(k_Mg24_He4_to_Si28_reaclib)*Y(He4)*state.rho;
    jac.set(Si28, Mg24, scratch);

    scratch = screened_rates(k_Al27_p_to_Si28_reaclib)*Y(H1)*state.rho;
    jac.set(Si28, Al27, scratch);

    scratch = screened_rates(k_Si27_n_to_Si28_reaclib)*Y(N)*state.rho;
    jac.set(Si28, Si27, scratch);

    scratch = -screened_rates(k_Si28_He4_to_C12_Ne20_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Si28_He4_to_O16_O16_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Si28_He4_to_S32_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Si28_He4_to_p_P31_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Si28_n_to_Si29_reaclib)*Y(N)*state.rho
              - screened_rates(k_Si28_p_to_P29_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Si28_to_He4_Mg24_reaclib)
              - screened_rates(k_Si28_to_n_Si27_reaclib)
              - screened_rates(k_Si28_to_p_Al27_reaclib);
    jac.set(Si28, Si28, scratch);

    scratch = screened_rates(k_Si29_to_n_Si28_reaclib);
    jac.set(Si28, Si29, scratch);

    jac.set(Si28, Si30, 0.0);

    jac.set(Si28, Si31, 0.0);

    scratch = screened_rates(k_P29_to_p_Si28_reaclib);
    jac.set(Si28, P29, scratch);

    scratch = screened_rates(k_P31_p_to_He4_Si28_reaclib)*Y(H1)*state.rho;
    jac.set(Si28, P31, scratch);

    jac.set(Si28, P33, 0.0);

    scratch = screened_rates(k_S32_to_He4_Si28_reaclib);
    jac.set(Si28, S32, scratch);

    jac.set(Si28, S33, 0.0);

    jac.set(Si28, S34, 0.0);

    jac.set(Si28, S35, 0.0);

    jac.set(Si28, Cl33, 0.0);

    jac.set(Si28, Cl35, 0.0);

    jac.set(Si28, Cl37, 0.0);

    jac.set(Si28, Ar36, 0.0);

    jac.set(Si28, Ar37, 0.0);

    jac.set(Si28, Ar38, 0.0);

    jac.set(Si28, K37, 0.0);

    jac.set(Si28, K39, 0.0);

    jac.set(Si28, Ca40, 0.0);

    jac.set(Si28, Ca42, 0.0);

    jac.set(Si28, Sc41, 0.0);

    jac.set(Si28, Sc43, 0.0);

    jac.set(Si28, Sc45, 0.0);

    jac.set(Si28, Ti44, 0.0);

    jac.set(Si28, Ti46, 0.0);

    jac.set(Si28, Ti48, 0.0);

    jac.set(Si28, V47, 0.0);

    jac.set(Si28, V49, 0.0);

    jac.set(Si28, V51, 0.0);

    jac.set(Si28, Cr48, 0.0);

    jac.set(Si28, Cr50, 0.0);

    jac.set(Si28, Cr52, 0.0);

    jac.set(Si28, Cr54, 0.0);

    jac.set(Si28, Mn51, 0.0);

    jac.set(Si28, Mn53, 0.0);

    jac.set(Si28, Mn55, 0.0);

    jac.set(Si28, Fe52, 0.0);

    jac.set(Si28, Fe54, 0.0);

    jac.set(Si28, Fe56, 0.0);

    jac.set(Si28, Fe58, 0.0);

    jac.set(Si28, Co55, 0.0);

    jac.set(Si28, Co56, 0.0);

    jac.set(Si28, Co57, 0.0);

    jac.set(Si28, Co59, 0.0);

    jac.set(Si28, Ni56, 0.0);

    jac.set(Si28, Ni58, 0.0);

    jac.set(Si28, Ni60, 0.0);

    jac.set(Si28, Cu59, 0.0);

    jac.set(Si28, Zn60, 0.0);


    // Si29 row

    scratch = screened_rates(k_P29_n_to_p_Si29_reaclib)*Y(P29)*state.rho
              + screened_rates(k_S32_n_to_He4_Si29_reaclib)*Y(S32)*state.rho
              + screened_rates(k_Si28_n_to_Si29_reaclib)*Y(Si28)*state.rho
              - screened_rates(k_Si29_n_to_Si30_reaclib)*Y(Si29)*state.rho;
    jac.set(Si29, N, scratch);

    scratch = -screened_rates(k_Si29_p_to_n_P29_reaclib)*Y(Si29)*state.rho;
    jac.set(Si29, H1, scratch);

    scratch = -screened_rates(k_Si29_He4_to_S33_reaclib)*Y(Si29)*state.rho
              - screened_rates(k_Si29_He4_to_n_S32_reaclib)*Y(Si29)*state.rho;
    jac.set(Si29, He4, scratch);

    jac.set(Si29, C12, 0.0);

    jac.set(Si29, O16, 0.0);

    jac.set(Si29, Ne20, 0.0);

    jac.set(Si29, Mg23, 0.0);

    jac.set(Si29, Mg24, 0.0);

    jac.set(Si29, Al27, 0.0);

    jac.set(Si29, Si27, 0.0);

    scratch = screened_rates(k_Si28_n_to_Si29_reaclib)*Y(N)*state.rho;
    jac.set(Si29, Si28, scratch);

    scratch = -screened_rates(k_Si29_He4_to_S33_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Si29_He4_to_n_S32_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Si29_n_to_Si30_reaclib)*Y(N)*state.rho
              - screened_rates(k_Si29_p_to_n_P29_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Si29_to_P29_beta_neg_weaktab)
              - screened_rates(k_Si29_to_n_Si28_reaclib);
    jac.set(Si29, Si29, scratch);

    scratch = screened_rates(k_Si30_to_n_Si29_reaclib);
    jac.set(Si29, Si30, scratch);

    jac.set(Si29, Si31, 0.0);

    scratch = screened_rates(k_P29_n_to_p_Si29_reaclib)*Y(N)*state.rho
              + screened_rates(k_P29_to_Si29_electron_capture_weaktab);
    jac.set(Si29, P29, scratch);

    jac.set(Si29, P31, 0.0);

    jac.set(Si29, P33, 0.0);

    scratch = screened_rates(k_S32_n_to_He4_Si29_reaclib)*Y(N)*state.rho;
    jac.set(Si29, S32, scratch);

    scratch = screened_rates(k_S33_to_He4_Si29_reaclib);
    jac.set(Si29, S33, scratch);

    jac.set(Si29, S34, 0.0);

    jac.set(Si29, S35, 0.0);

    jac.set(Si29, Cl33, 0.0);

    jac.set(Si29, Cl35, 0.0);

    jac.set(Si29, Cl37, 0.0);

    jac.set(Si29, Ar36, 0.0);

    jac.set(Si29, Ar37, 0.0);

    jac.set(Si29, Ar38, 0.0);

    jac.set(Si29, K37, 0.0);

    jac.set(Si29, K39, 0.0);

    jac.set(Si29, Ca40, 0.0);

    jac.set(Si29, Ca42, 0.0);

    jac.set(Si29, Sc41, 0.0);

    jac.set(Si29, Sc43, 0.0);

    jac.set(Si29, Sc45, 0.0);

    jac.set(Si29, Ti44, 0.0);

    jac.set(Si29, Ti46, 0.0);

    jac.set(Si29, Ti48, 0.0);

    jac.set(Si29, V47, 0.0);

    jac.set(Si29, V49, 0.0);

    jac.set(Si29, V51, 0.0);

    jac.set(Si29, Cr48, 0.0);

    jac.set(Si29, Cr50, 0.0);

    jac.set(Si29, Cr52, 0.0);

    jac.set(Si29, Cr54, 0.0);

    jac.set(Si29, Mn51, 0.0);

    jac.set(Si29, Mn53, 0.0);

    jac.set(Si29, Mn55, 0.0);

    jac.set(Si29, Fe52, 0.0);

    jac.set(Si29, Fe54, 0.0);

    jac.set(Si29, Fe56, 0.0);

    jac.set(Si29, Fe58, 0.0);

    jac.set(Si29, Co55, 0.0);

    jac.set(Si29, Co56, 0.0);

    jac.set(Si29, Co57, 0.0);

    jac.set(Si29, Co59, 0.0);

    jac.set(Si29, Ni56, 0.0);

    jac.set(Si29, Ni58, 0.0);

    jac.set(Si29, Ni60, 0.0);

    jac.set(Si29, Cu59, 0.0);

    jac.set(Si29, Zn60, 0.0);


    // Si30 row

    scratch = screened_rates(k_S33_n_to_He4_Si30_reaclib)*Y(S33)*state.rho
              + screened_rates(k_Si29_n_to_Si30_reaclib)*Y(Si29)*state.rho
              - screened_rates(k_Si30_n_to_Si31_reaclib)*Y(Si30)*state.rho;
    jac.set(Si30, N, scratch);

    scratch = screened_rates(k_P33_p_to_He4_Si30_reaclib)*Y(P33)*state.rho
              - screened_rates(k_Si30_p_to_He4_Al27_reaclib)*Y(Si30)*state.rho
              - screened_rates(k_Si30_p_to_P31_reaclib)*Y(Si30)*state.rho;
    jac.set(Si30, H1, scratch);

    scratch = screened_rates(k_Al27_He4_to_p_Si30_reaclib)*Y(Al27)*state.rho
              - screened_rates(k_Si30_He4_to_S34_reaclib)*Y(Si30)*state.rho
              - screened_rates(k_Si30_He4_to_n_S33_reaclib)*Y(Si30)*state.rho
              - screened_rates(k_Si30_He4_to_p_P33_reaclib)*Y(Si30)*state.rho;
    jac.set(Si30, He4, scratch);

    jac.set(Si30, C12, 0.0);

    jac.set(Si30, O16, 0.0);

    jac.set(Si30, Ne20, 0.0);

    jac.set(Si30, Mg23, 0.0);

    jac.set(Si30, Mg24, 0.0);

    scratch = screened_rates(k_Al27_He4_to_p_Si30_reaclib)*Y(He4)*state.rho;
    jac.set(Si30, Al27, scratch);

    jac.set(Si30, Si27, 0.0);

    jac.set(Si30, Si28, 0.0);

    scratch = screened_rates(k_Si29_n_to_Si30_reaclib)*Y(N)*state.rho;
    jac.set(Si30, Si29, scratch);

    scratch = -screened_rates(k_Si30_He4_to_S34_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Si30_He4_to_n_S33_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Si30_He4_to_p_P33_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Si30_n_to_Si31_reaclib)*Y(N)*state.rho
              - screened_rates(k_Si30_p_to_He4_Al27_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Si30_p_to_P31_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Si30_to_n_Si29_reaclib);
    jac.set(Si30, Si30, scratch);

    scratch = screened_rates(k_Si31_to_n_Si30_reaclib);
    jac.set(Si30, Si31, scratch);

    jac.set(Si30, P29, 0.0);

    scratch = screened_rates(k_P31_to_p_Si30_reaclib);
    jac.set(Si30, P31, scratch);

    scratch = screened_rates(k_P33_p_to_He4_Si30_reaclib)*Y(H1)*state.rho;
    jac.set(Si30, P33, scratch);

    jac.set(Si30, S32, 0.0);

    scratch = screened_rates(k_S33_n_to_He4_Si30_reaclib)*Y(N)*state.rho;
    jac.set(Si30, S33, scratch);

    scratch = screened_rates(k_S34_to_He4_Si30_reaclib);
    jac.set(Si30, S34, scratch);

    jac.set(Si30, S35, 0.0);

    jac.set(Si30, Cl33, 0.0);

    jac.set(Si30, Cl35, 0.0);

    jac.set(Si30, Cl37, 0.0);

    jac.set(Si30, Ar36, 0.0);

    jac.set(Si30, Ar37, 0.0);

    jac.set(Si30, Ar38, 0.0);

    jac.set(Si30, K37, 0.0);

    jac.set(Si30, K39, 0.0);

    jac.set(Si30, Ca40, 0.0);

    jac.set(Si30, Ca42, 0.0);

    jac.set(Si30, Sc41, 0.0);

    jac.set(Si30, Sc43, 0.0);

    jac.set(Si30, Sc45, 0.0);

    jac.set(Si30, Ti44, 0.0);

    jac.set(Si30, Ti46, 0.0);

    jac.set(Si30, Ti48, 0.0);

    jac.set(Si30, V47, 0.0);

    jac.set(Si30, V49, 0.0);

    jac.set(Si30, V51, 0.0);

    jac.set(Si30, Cr48, 0.0);

    jac.set(Si30, Cr50, 0.0);

    jac.set(Si30, Cr52, 0.0);

    jac.set(Si30, Cr54, 0.0);

    jac.set(Si30, Mn51, 0.0);

    jac.set(Si30, Mn53, 0.0);

    jac.set(Si30, Mn55, 0.0);

    jac.set(Si30, Fe52, 0.0);

    jac.set(Si30, Fe54, 0.0);

    jac.set(Si30, Fe56, 0.0);

    jac.set(Si30, Fe58, 0.0);

    jac.set(Si30, Co55, 0.0);

    jac.set(Si30, Co56, 0.0);

    jac.set(Si30, Co57, 0.0);

    jac.set(Si30, Co59, 0.0);

    jac.set(Si30, Ni56, 0.0);

    jac.set(Si30, Ni58, 0.0);

    jac.set(Si30, Ni60, 0.0);

    jac.set(Si30, Cu59, 0.0);

    jac.set(Si30, Zn60, 0.0);


    // Si31 row

    scratch = screened_rates(k_P31_n_to_p_Si31_reaclib)*Y(P31)*state.rho
              + screened_rates(k_S34_n_to_He4_Si31_reaclib)*Y(S34)*state.rho
              + screened_rates(k_Si30_n_to_Si31_reaclib)*Y(Si30)*state.rho;
    jac.set(Si31, N, scratch);

    scratch = -screened_rates(k_Si31_p_to_n_P31_reaclib)*Y(Si31)*state.rho;
    jac.set(Si31, H1, scratch);

    scratch = -screened_rates(k_Si31_He4_to_S35_reaclib)*Y(Si31)*state.rho
              - screened_rates(k_Si31_He4_to_n_S34_reaclib)*Y(Si31)*state.rho;
    jac.set(Si31, He4, scratch);

    jac.set(Si31, C12, 0.0);

    jac.set(Si31, O16, 0.0);

    jac.set(Si31, Ne20, 0.0);

    jac.set(Si31, Mg23, 0.0);

    jac.set(Si31, Mg24, 0.0);

    jac.set(Si31, Al27, 0.0);

    jac.set(Si31, Si27, 0.0);

    jac.set(Si31, Si28, 0.0);

    jac.set(Si31, Si29, 0.0);

    scratch = screened_rates(k_Si30_n_to_Si31_reaclib)*Y(N)*state.rho;
    jac.set(Si31, Si30, scratch);

    scratch = -screened_rates(k_Si31_He4_to_S35_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Si31_He4_to_n_S34_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Si31_p_to_n_P31_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Si31_to_P31_beta_neg_weaktab)
              - screened_rates(k_Si31_to_n_Si30_reaclib);
    jac.set(Si31, Si31, scratch);

    jac.set(Si31, P29, 0.0);

    scratch = screened_rates(k_P31_n_to_p_Si31_reaclib)*Y(N)*state.rho
              + screened_rates(k_P31_to_Si31_electron_capture_weaktab);
    jac.set(Si31, P31, scratch);

    jac.set(Si31, P33, 0.0);

    jac.set(Si31, S32, 0.0);

    jac.set(Si31, S33, 0.0);

    scratch = screened_rates(k_S34_n_to_He4_Si31_reaclib)*Y(N)*state.rho;
    jac.set(Si31, S34, scratch);

    scratch = screened_rates(k_S35_to_He4_Si31_reaclib);
    jac.set(Si31, S35, scratch);

    jac.set(Si31, Cl33, 0.0);

    jac.set(Si31, Cl35, 0.0);

    jac.set(Si31, Cl37, 0.0);

    jac.set(Si31, Ar36, 0.0);

    jac.set(Si31, Ar37, 0.0);

    jac.set(Si31, Ar38, 0.0);

    jac.set(Si31, K37, 0.0);

    jac.set(Si31, K39, 0.0);

    jac.set(Si31, Ca40, 0.0);

    jac.set(Si31, Ca42, 0.0);

    jac.set(Si31, Sc41, 0.0);

    jac.set(Si31, Sc43, 0.0);

    jac.set(Si31, Sc45, 0.0);

    jac.set(Si31, Ti44, 0.0);

    jac.set(Si31, Ti46, 0.0);

    jac.set(Si31, Ti48, 0.0);

    jac.set(Si31, V47, 0.0);

    jac.set(Si31, V49, 0.0);

    jac.set(Si31, V51, 0.0);

    jac.set(Si31, Cr48, 0.0);

    jac.set(Si31, Cr50, 0.0);

    jac.set(Si31, Cr52, 0.0);

    jac.set(Si31, Cr54, 0.0);

    jac.set(Si31, Mn51, 0.0);

    jac.set(Si31, Mn53, 0.0);

    jac.set(Si31, Mn55, 0.0);

    jac.set(Si31, Fe52, 0.0);

    jac.set(Si31, Fe54, 0.0);

    jac.set(Si31, Fe56, 0.0);

    jac.set(Si31, Fe58, 0.0);

    jac.set(Si31, Co55, 0.0);

    jac.set(Si31, Co56, 0.0);

    jac.set(Si31, Co57, 0.0);

    jac.set(Si31, Co59, 0.0);

    jac.set(Si31, Ni56, 0.0);

    jac.set(Si31, Ni58, 0.0);

    jac.set(Si31, Ni60, 0.0);

    jac.set(Si31, Cu59, 0.0);

    jac.set(Si31, Zn60, 0.0);


    // P29 row

    scratch = -screened_rates(k_P29_n_to_p_Si29_reaclib)*Y(P29)*state.rho;
    jac.set(P29, N, scratch);

    scratch = screened_rates(k_S32_p_to_He4_P29_reaclib)*Y(S32)*state.rho
              + screened_rates(k_Si28_p_to_P29_reaclib)*Y(Si28)*state.rho
              + screened_rates(k_Si29_p_to_n_P29_reaclib)*Y(Si29)*state.rho;
    jac.set(P29, H1, scratch);

    scratch = -screened_rates(k_P29_He4_to_Cl33_reaclib)*Y(P29)*state.rho
              - screened_rates(k_P29_He4_to_p_S32_reaclib)*Y(P29)*state.rho;
    jac.set(P29, He4, scratch);

    jac.set(P29, C12, 0.0);

    jac.set(P29, O16, 0.0);

    jac.set(P29, Ne20, 0.0);

    jac.set(P29, Mg23, 0.0);

    jac.set(P29, Mg24, 0.0);

    jac.set(P29, Al27, 0.0);

    jac.set(P29, Si27, 0.0);

    scratch = screened_rates(k_Si28_p_to_P29_reaclib)*Y(H1)*state.rho;
    jac.set(P29, Si28, scratch);

    scratch = screened_rates(k_Si29_p_to_n_P29_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Si29_to_P29_beta_neg_weaktab);
    jac.set(P29, Si29, scratch);

    jac.set(P29, Si30, 0.0);

    jac.set(P29, Si31, 0.0);

    scratch = -screened_rates(k_P29_He4_to_Cl33_reaclib)*Y(He4)*state.rho
              - screened_rates(k_P29_He4_to_p_S32_reaclib)*Y(He4)*state.rho
              - screened_rates(k_P29_n_to_p_Si29_reaclib)*Y(N)*state.rho
              - screened_rates(k_P29_to_Si29_electron_capture_weaktab)
              - screened_rates(k_P29_to_p_Si28_reaclib);
    jac.set(P29, P29, scratch);

    jac.set(P29, P31, 0.0);

    jac.set(P29, P33, 0.0);

    scratch = screened_rates(k_S32_p_to_He4_P29_reaclib)*Y(H1)*state.rho;
    jac.set(P29, S32, scratch);

    jac.set(P29, S33, 0.0);

    jac.set(P29, S34, 0.0);

    jac.set(P29, S35, 0.0);

    scratch = screened_rates(k_Cl33_to_He4_P29_reaclib);
    jac.set(P29, Cl33, scratch);

    jac.set(P29, Cl35, 0.0);

    jac.set(P29, Cl37, 0.0);

    jac.set(P29, Ar36, 0.0);

    jac.set(P29, Ar37, 0.0);

    jac.set(P29, Ar38, 0.0);

    jac.set(P29, K37, 0.0);

    jac.set(P29, K39, 0.0);

    jac.set(P29, Ca40, 0.0);

    jac.set(P29, Ca42, 0.0);

    jac.set(P29, Sc41, 0.0);

    jac.set(P29, Sc43, 0.0);

    jac.set(P29, Sc45, 0.0);

    jac.set(P29, Ti44, 0.0);

    jac.set(P29, Ti46, 0.0);

    jac.set(P29, Ti48, 0.0);

    jac.set(P29, V47, 0.0);

    jac.set(P29, V49, 0.0);

    jac.set(P29, V51, 0.0);

    jac.set(P29, Cr48, 0.0);

    jac.set(P29, Cr50, 0.0);

    jac.set(P29, Cr52, 0.0);

    jac.set(P29, Cr54, 0.0);

    jac.set(P29, Mn51, 0.0);

    jac.set(P29, Mn53, 0.0);

    jac.set(P29, Mn55, 0.0);

    jac.set(P29, Fe52, 0.0);

    jac.set(P29, Fe54, 0.0);

    jac.set(P29, Fe56, 0.0);

    jac.set(P29, Fe58, 0.0);

    jac.set(P29, Co55, 0.0);

    jac.set(P29, Co56, 0.0);

    jac.set(P29, Co57, 0.0);

    jac.set(P29, Co59, 0.0);

    jac.set(P29, Ni56, 0.0);

    jac.set(P29, Ni58, 0.0);

    jac.set(P29, Ni60, 0.0);

    jac.set(P29, Cu59, 0.0);

    jac.set(P29, Zn60, 0.0);


    // P31 row

    scratch = -screened_rates(k_P31_n_to_p_Si31_reaclib)*Y(P31)*state.rho;
    jac.set(P31, N, scratch);

    scratch = -screened_rates(k_P31_p_to_C12_Ne20_reaclib)*Y(P31)*state.rho
              - screened_rates(k_P31_p_to_He4_Si28_reaclib)*Y(P31)*state.rho
              - screened_rates(k_P31_p_to_O16_O16_reaclib)*Y(P31)*state.rho
              - screened_rates(k_P31_p_to_S32_reaclib)*Y(P31)*state.rho
              + screened_rates(k_S34_p_to_He4_P31_reaclib)*Y(S34)*state.rho
              + screened_rates(k_Si30_p_to_P31_reaclib)*Y(Si30)*state.rho
              + screened_rates(k_Si31_p_to_n_P31_reaclib)*Y(Si31)*state.rho;
    jac.set(P31, H1, scratch);

    scratch = screened_rates(k_Al27_He4_to_P31_reaclib)*Y(Al27)*state.rho
              - screened_rates(k_P31_He4_to_Cl35_reaclib)*Y(P31)*state.rho
              - screened_rates(k_P31_He4_to_p_S34_reaclib)*Y(P31)*state.rho
              + screened_rates(k_Si28_He4_to_p_P31_reaclib)*Y(Si28)*state.rho;
    jac.set(P31, He4, scratch);

    scratch = screened_rates(k_Ne20_C12_to_p_P31_reaclib)*Y(Ne20)*state.rho;
    jac.set(P31, C12, scratch);

    scratch = 1.0*screened_rates(k_O16_O16_to_p_P31_reaclib)*Y(O16)*state.rho;
    jac.set(P31, O16, scratch);

    scratch = screened_rates(k_Ne20_C12_to_p_P31_reaclib)*Y(C12)*state.rho;
    jac.set(P31, Ne20, scratch);

    jac.set(P31, Mg23, 0.0);

    jac.set(P31, Mg24, 0.0);

    scratch = screened_rates(k_Al27_He4_to_P31_reaclib)*Y(He4)*state.rho;
    jac.set(P31, Al27, scratch);

    jac.set(P31, Si27, 0.0);

    scratch = screened_rates(k_Si28_He4_to_p_P31_reaclib)*Y(He4)*state.rho;
    jac.set(P31, Si28, scratch);

    jac.set(P31, Si29, 0.0);

    scratch = screened_rates(k_Si30_p_to_P31_reaclib)*Y(H1)*state.rho;
    jac.set(P31, Si30, scratch);

    scratch = screened_rates(k_Si31_p_to_n_P31_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Si31_to_P31_beta_neg_weaktab);
    jac.set(P31, Si31, scratch);

    jac.set(P31, P29, 0.0);

    scratch = -screened_rates(k_P31_He4_to_Cl35_reaclib)*Y(He4)*state.rho
              - screened_rates(k_P31_He4_to_p_S34_reaclib)*Y(He4)*state.rho
              - screened_rates(k_P31_n_to_p_Si31_reaclib)*Y(N)*state.rho
              - screened_rates(k_P31_p_to_C12_Ne20_reaclib)*Y(H1)*state.rho
              - screened_rates(k_P31_p_to_He4_Si28_reaclib)*Y(H1)*state.rho
              - screened_rates(k_P31_p_to_O16_O16_reaclib)*Y(H1)*state.rho
              - screened_rates(k_P31_p_to_S32_reaclib)*Y(H1)*state.rho
              - screened_rates(k_P31_to_He4_Al27_reaclib)
              - screened_rates(k_P31_to_Si31_electron_capture_weaktab)
              - screened_rates(k_P31_to_p_Si30_reaclib);
    jac.set(P31, P31, scratch);

    jac.set(P31, P33, 0.0);

    scratch = screened_rates(k_S32_to_p_P31_reaclib);
    jac.set(P31, S32, scratch);

    jac.set(P31, S33, 0.0);

    scratch = screened_rates(k_S34_p_to_He4_P31_reaclib)*Y(H1)*state.rho;
    jac.set(P31, S34, scratch);

    jac.set(P31, S35, 0.0);

    jac.set(P31, Cl33, 0.0);

    scratch = screened_rates(k_Cl35_to_He4_P31_reaclib);
    jac.set(P31, Cl35, scratch);

    jac.set(P31, Cl37, 0.0);

    jac.set(P31, Ar36, 0.0);

    jac.set(P31, Ar37, 0.0);

    jac.set(P31, Ar38, 0.0);

    jac.set(P31, K37, 0.0);

    jac.set(P31, K39, 0.0);

    jac.set(P31, Ca40, 0.0);

    jac.set(P31, Ca42, 0.0);

    jac.set(P31, Sc41, 0.0);

    jac.set(P31, Sc43, 0.0);

    jac.set(P31, Sc45, 0.0);

    jac.set(P31, Ti44, 0.0);

    jac.set(P31, Ti46, 0.0);

    jac.set(P31, Ti48, 0.0);

    jac.set(P31, V47, 0.0);

    jac.set(P31, V49, 0.0);

    jac.set(P31, V51, 0.0);

    jac.set(P31, Cr48, 0.0);

    jac.set(P31, Cr50, 0.0);

    jac.set(P31, Cr52, 0.0);

    jac.set(P31, Cr54, 0.0);

    jac.set(P31, Mn51, 0.0);

    jac.set(P31, Mn53, 0.0);

    jac.set(P31, Mn55, 0.0);

    jac.set(P31, Fe52, 0.0);

    jac.set(P31, Fe54, 0.0);

    jac.set(P31, Fe56, 0.0);

    jac.set(P31, Fe58, 0.0);

    jac.set(P31, Co55, 0.0);

    jac.set(P31, Co56, 0.0);

    jac.set(P31, Co57, 0.0);

    jac.set(P31, Co59, 0.0);

    jac.set(P31, Ni56, 0.0);

    jac.set(P31, Ni58, 0.0);

    jac.set(P31, Ni60, 0.0);

    jac.set(P31, Cu59, 0.0);

    jac.set(P31, Zn60, 0.0);


    // P33 row

    scratch = screened_rates(k_S33_n_to_p_P33_reaclib)*Y(S33)*state.rho;
    jac.set(P33, N, scratch);

    scratch = -screened_rates(k_P33_p_to_He4_Si30_reaclib)*Y(P33)*state.rho
              - screened_rates(k_P33_p_to_S34_reaclib)*Y(P33)*state.rho
              - screened_rates(k_P33_p_to_n_S33_reaclib)*Y(P33)*state.rho;
    jac.set(P33, H1, scratch);

    scratch = -screened_rates(k_P33_He4_to_Cl37_reaclib)*Y(P33)*state.rho
              + screened_rates(k_Si30_He4_to_p_P33_reaclib)*Y(Si30)*state.rho;
    jac.set(P33, He4, scratch);

    jac.set(P33, C12, 0.0);

    jac.set(P33, O16, 0.0);

    jac.set(P33, Ne20, 0.0);

    jac.set(P33, Mg23, 0.0);

    jac.set(P33, Mg24, 0.0);

    jac.set(P33, Al27, 0.0);

    jac.set(P33, Si27, 0.0);

    jac.set(P33, Si28, 0.0);

    jac.set(P33, Si29, 0.0);

    scratch = screened_rates(k_Si30_He4_to_p_P33_reaclib)*Y(He4)*state.rho;
    jac.set(P33, Si30, scratch);

    jac.set(P33, Si31, 0.0);

    jac.set(P33, P29, 0.0);

    jac.set(P33, P31, 0.0);

    scratch = -screened_rates(k_P33_He4_to_Cl37_reaclib)*Y(He4)*state.rho
              - screened_rates(k_P33_p_to_He4_Si30_reaclib)*Y(H1)*state.rho
              - screened_rates(k_P33_p_to_S34_reaclib)*Y(H1)*state.rho
              - screened_rates(k_P33_p_to_n_S33_reaclib)*Y(H1)*state.rho
              - screened_rates(k_P33_to_S33_beta_neg_weaktab);
    jac.set(P33, P33, scratch);

    jac.set(P33, S32, 0.0);

    scratch = screened_rates(k_S33_n_to_p_P33_reaclib)*Y(N)*state.rho
              + screened_rates(k_S33_to_P33_electron_capture_weaktab);
    jac.set(P33, S33, scratch);

    scratch = screened_rates(k_S34_to_p_P33_reaclib);
    jac.set(P33, S34, scratch);

    jac.set(P33, S35, 0.0);

    jac.set(P33, Cl33, 0.0);

    jac.set(P33, Cl35, 0.0);

    scratch = screened_rates(k_Cl37_to_He4_P33_reaclib);
    jac.set(P33, Cl37, scratch);

    jac.set(P33, Ar36, 0.0);

    jac.set(P33, Ar37, 0.0);

    jac.set(P33, Ar38, 0.0);

    jac.set(P33, K37, 0.0);

    jac.set(P33, K39, 0.0);

    jac.set(P33, Ca40, 0.0);

    jac.set(P33, Ca42, 0.0);

    jac.set(P33, Sc41, 0.0);

    jac.set(P33, Sc43, 0.0);

    jac.set(P33, Sc45, 0.0);

    jac.set(P33, Ti44, 0.0);

    jac.set(P33, Ti46, 0.0);

    jac.set(P33, Ti48, 0.0);

    jac.set(P33, V47, 0.0);

    jac.set(P33, V49, 0.0);

    jac.set(P33, V51, 0.0);

    jac.set(P33, Cr48, 0.0);

    jac.set(P33, Cr50, 0.0);

    jac.set(P33, Cr52, 0.0);

    jac.set(P33, Cr54, 0.0);

    jac.set(P33, Mn51, 0.0);

    jac.set(P33, Mn53, 0.0);

    jac.set(P33, Mn55, 0.0);

    jac.set(P33, Fe52, 0.0);

    jac.set(P33, Fe54, 0.0);

    jac.set(P33, Fe56, 0.0);

    jac.set(P33, Fe58, 0.0);

    jac.set(P33, Co55, 0.0);

    jac.set(P33, Co56, 0.0);

    jac.set(P33, Co57, 0.0);

    jac.set(P33, Co59, 0.0);

    jac.set(P33, Ni56, 0.0);

    jac.set(P33, Ni58, 0.0);

    jac.set(P33, Ni60, 0.0);

    jac.set(P33, Cu59, 0.0);

    jac.set(P33, Zn60, 0.0);


    // S32 row

    scratch = -screened_rates(k_S32_n_to_He4_Si29_reaclib)*Y(S32)*state.rho
              - screened_rates(k_S32_n_to_S33_reaclib)*Y(S32)*state.rho;
    jac.set(S32, N, scratch);

    scratch = screened_rates(k_Cl35_p_to_He4_S32_reaclib)*Y(Cl35)*state.rho
              + screened_rates(k_P31_p_to_S32_reaclib)*Y(P31)*state.rho
              - screened_rates(k_S32_p_to_Cl33_reaclib)*Y(S32)*state.rho
              - screened_rates(k_S32_p_to_He4_P29_reaclib)*Y(S32)*state.rho;
    jac.set(S32, H1, scratch);

    scratch = screened_rates(k_P29_He4_to_p_S32_reaclib)*Y(P29)*state.rho
              - screened_rates(k_S32_He4_to_Ar36_reaclib)*Y(S32)*state.rho
              - screened_rates(k_S32_He4_to_p_Cl35_reaclib)*Y(S32)*state.rho
              + screened_rates(k_Si28_He4_to_S32_reaclib)*Y(Si28)*state.rho
              + screened_rates(k_Si29_He4_to_n_S32_reaclib)*Y(Si29)*state.rho;
    jac.set(S32, He4, scratch);

    jac.set(S32, C12, 0.0);

    jac.set(S32, O16, 0.0);

    jac.set(S32, Ne20, 0.0);

    jac.set(S32, Mg23, 0.0);

    jac.set(S32, Mg24, 0.0);

    jac.set(S32, Al27, 0.0);

    jac.set(S32, Si27, 0.0);

    scratch = screened_rates(k_Si28_He4_to_S32_reaclib)*Y(He4)*state.rho;
    jac.set(S32, Si28, scratch);

    scratch = screened_rates(k_Si29_He4_to_n_S32_reaclib)*Y(He4)*state.rho;
    jac.set(S32, Si29, scratch);

    jac.set(S32, Si30, 0.0);

    jac.set(S32, Si31, 0.0);

    scratch = screened_rates(k_P29_He4_to_p_S32_reaclib)*Y(He4)*state.rho;
    jac.set(S32, P29, scratch);

    scratch = screened_rates(k_P31_p_to_S32_reaclib)*Y(H1)*state.rho;
    jac.set(S32, P31, scratch);

    jac.set(S32, P33, 0.0);

    scratch = -screened_rates(k_S32_He4_to_Ar36_reaclib)*Y(He4)*state.rho
              - screened_rates(k_S32_He4_to_p_Cl35_reaclib)*Y(He4)*state.rho
              - screened_rates(k_S32_n_to_He4_Si29_reaclib)*Y(N)*state.rho
              - screened_rates(k_S32_n_to_S33_reaclib)*Y(N)*state.rho
              - screened_rates(k_S32_p_to_Cl33_reaclib)*Y(H1)*state.rho
              - screened_rates(k_S32_p_to_He4_P29_reaclib)*Y(H1)*state.rho
              - screened_rates(k_S32_to_He4_Si28_reaclib)
              - screened_rates(k_S32_to_p_P31_reaclib);
    jac.set(S32, S32, scratch);

    scratch = screened_rates(k_S33_to_n_S32_reaclib);
    jac.set(S32, S33, scratch);

    jac.set(S32, S34, 0.0);

    jac.set(S32, S35, 0.0);

    scratch = screened_rates(k_Cl33_to_p_S32_reaclib);
    jac.set(S32, Cl33, scratch);

    scratch = screened_rates(k_Cl35_p_to_He4_S32_reaclib)*Y(H1)*state.rho;
    jac.set(S32, Cl35, scratch);

    jac.set(S32, Cl37, 0.0);

    scratch = screened_rates(k_Ar36_to_He4_S32_reaclib);
    jac.set(S32, Ar36, scratch);

    jac.set(S32, Ar37, 0.0);

    jac.set(S32, Ar38, 0.0);

    jac.set(S32, K37, 0.0);

    jac.set(S32, K39, 0.0);

    jac.set(S32, Ca40, 0.0);

    jac.set(S32, Ca42, 0.0);

    jac.set(S32, Sc41, 0.0);

    jac.set(S32, Sc43, 0.0);

    jac.set(S32, Sc45, 0.0);

    jac.set(S32, Ti44, 0.0);

    jac.set(S32, Ti46, 0.0);

    jac.set(S32, Ti48, 0.0);

    jac.set(S32, V47, 0.0);

    jac.set(S32, V49, 0.0);

    jac.set(S32, V51, 0.0);

    jac.set(S32, Cr48, 0.0);

    jac.set(S32, Cr50, 0.0);

    jac.set(S32, Cr52, 0.0);

    jac.set(S32, Cr54, 0.0);

    jac.set(S32, Mn51, 0.0);

    jac.set(S32, Mn53, 0.0);

    jac.set(S32, Mn55, 0.0);

    jac.set(S32, Fe52, 0.0);

    jac.set(S32, Fe54, 0.0);

    jac.set(S32, Fe56, 0.0);

    jac.set(S32, Fe58, 0.0);

    jac.set(S32, Co55, 0.0);

    jac.set(S32, Co56, 0.0);

    jac.set(S32, Co57, 0.0);

    jac.set(S32, Co59, 0.0);

    jac.set(S32, Ni56, 0.0);

    jac.set(S32, Ni58, 0.0);

    jac.set(S32, Ni60, 0.0);

    jac.set(S32, Cu59, 0.0);

    jac.set(S32, Zn60, 0.0);


    // S33 row

    scratch = screened_rates(k_Ar36_n_to_He4_S33_reaclib)*Y(Ar36)*state.rho
              + screened_rates(k_Cl33_n_to_p_S33_reaclib)*Y(Cl33)*state.rho
              + screened_rates(k_S32_n_to_S33_reaclib)*Y(S32)*state.rho
              - screened_rates(k_S33_n_to_He4_Si30_reaclib)*Y(S33)*state.rho
              - screened_rates(k_S33_n_to_S34_reaclib)*Y(S33)*state.rho
              - screened_rates(k_S33_n_to_p_P33_reaclib)*Y(S33)*state.rho;
    jac.set(S33, N, scratch);

    scratch = screened_rates(k_P33_p_to_n_S33_reaclib)*Y(P33)*state.rho
              - screened_rates(k_S33_p_to_n_Cl33_reaclib)*Y(S33)*state.rho;
    jac.set(S33, H1, scratch);

    scratch = -screened_rates(k_S33_He4_to_Ar37_reaclib)*Y(S33)*state.rho
              - screened_rates(k_S33_He4_to_n_Ar36_reaclib)*Y(S33)*state.rho
              + screened_rates(k_Si29_He4_to_S33_reaclib)*Y(Si29)*state.rho
              + screened_rates(k_Si30_He4_to_n_S33_reaclib)*Y(Si30)*state.rho;
    jac.set(S33, He4, scratch);

    jac.set(S33, C12, 0.0);

    jac.set(S33, O16, 0.0);

    jac.set(S33, Ne20, 0.0);

    jac.set(S33, Mg23, 0.0);

    jac.set(S33, Mg24, 0.0);

    jac.set(S33, Al27, 0.0);

    jac.set(S33, Si27, 0.0);

    jac.set(S33, Si28, 0.0);

    scratch = screened_rates(k_Si29_He4_to_S33_reaclib)*Y(He4)*state.rho;
    jac.set(S33, Si29, scratch);

    scratch = screened_rates(k_Si30_He4_to_n_S33_reaclib)*Y(He4)*state.rho;
    jac.set(S33, Si30, scratch);

    jac.set(S33, Si31, 0.0);

    jac.set(S33, P29, 0.0);

    jac.set(S33, P31, 0.0);

    scratch = screened_rates(k_P33_p_to_n_S33_reaclib)*Y(H1)*state.rho
              + screened_rates(k_P33_to_S33_beta_neg_weaktab);
    jac.set(S33, P33, scratch);

    scratch = screened_rates(k_S32_n_to_S33_reaclib)*Y(N)*state.rho;
    jac.set(S33, S32, scratch);

    scratch = -screened_rates(k_S33_He4_to_Ar37_reaclib)*Y(He4)*state.rho
              - screened_rates(k_S33_He4_to_n_Ar36_reaclib)*Y(He4)*state.rho
              - screened_rates(k_S33_n_to_He4_Si30_reaclib)*Y(N)*state.rho
              - screened_rates(k_S33_n_to_S34_reaclib)*Y(N)*state.rho
              - screened_rates(k_S33_n_to_p_P33_reaclib)*Y(N)*state.rho
              - screened_rates(k_S33_p_to_n_Cl33_reaclib)*Y(H1)*state.rho
              - screened_rates(k_S33_to_Cl33_beta_neg_weaktab)
              - screened_rates(k_S33_to_He4_Si29_reaclib)
              - screened_rates(k_S33_to_P33_electron_capture_weaktab)
              - screened_rates(k_S33_to_n_S32_reaclib);
    jac.set(S33, S33, scratch);

    scratch = screened_rates(k_S34_to_n_S33_reaclib);
    jac.set(S33, S34, scratch);

    jac.set(S33, S35, 0.0);

    scratch = screened_rates(k_Cl33_n_to_p_S33_reaclib)*Y(N)*state.rho
              + screened_rates(k_Cl33_to_S33_electron_capture_weaktab);
    jac.set(S33, Cl33, scratch);

    jac.set(S33, Cl35, 0.0);

    jac.set(S33, Cl37, 0.0);

    scratch = screened_rates(k_Ar36_n_to_He4_S33_reaclib)*Y(N)*state.rho;
    jac.set(S33, Ar36, scratch);

    scratch = screened_rates(k_Ar37_to_He4_S33_reaclib);
    jac.set(S33, Ar37, scratch);

    jac.set(S33, Ar38, 0.0);

    jac.set(S33, K37, 0.0);

    jac.set(S33, K39, 0.0);

    jac.set(S33, Ca40, 0.0);

    jac.set(S33, Ca42, 0.0);

    jac.set(S33, Sc41, 0.0);

    jac.set(S33, Sc43, 0.0);

    jac.set(S33, Sc45, 0.0);

    jac.set(S33, Ti44, 0.0);

    jac.set(S33, Ti46, 0.0);

    jac.set(S33, Ti48, 0.0);

    jac.set(S33, V47, 0.0);

    jac.set(S33, V49, 0.0);

    jac.set(S33, V51, 0.0);

    jac.set(S33, Cr48, 0.0);

    jac.set(S33, Cr50, 0.0);

    jac.set(S33, Cr52, 0.0);

    jac.set(S33, Cr54, 0.0);

    jac.set(S33, Mn51, 0.0);

    jac.set(S33, Mn53, 0.0);

    jac.set(S33, Mn55, 0.0);

    jac.set(S33, Fe52, 0.0);

    jac.set(S33, Fe54, 0.0);

    jac.set(S33, Fe56, 0.0);

    jac.set(S33, Fe58, 0.0);

    jac.set(S33, Co55, 0.0);

    jac.set(S33, Co56, 0.0);

    jac.set(S33, Co57, 0.0);

    jac.set(S33, Co59, 0.0);

    jac.set(S33, Ni56, 0.0);

    jac.set(S33, Ni58, 0.0);

    jac.set(S33, Ni60, 0.0);

    jac.set(S33, Cu59, 0.0);

    jac.set(S33, Zn60, 0.0);


    // S34 row

    scratch = screened_rates(k_Ar37_n_to_He4_S34_reaclib)*Y(Ar37)*state.rho
              + screened_rates(k_S33_n_to_S34_reaclib)*Y(S33)*state.rho
              - screened_rates(k_S34_n_to_He4_Si31_reaclib)*Y(S34)*state.rho
              - screened_rates(k_S34_n_to_S35_reaclib)*Y(S34)*state.rho;
    jac.set(S34, N, scratch);

    scratch = screened_rates(k_Cl37_p_to_He4_S34_reaclib)*Y(Cl37)*state.rho
              + screened_rates(k_P33_p_to_S34_reaclib)*Y(P33)*state.rho
              - screened_rates(k_S34_p_to_Cl35_reaclib)*Y(S34)*state.rho
              - screened_rates(k_S34_p_to_He4_P31_reaclib)*Y(S34)*state.rho;
    jac.set(S34, H1, scratch);

    scratch = screened_rates(k_P31_He4_to_p_S34_reaclib)*Y(P31)*state.rho
              - screened_rates(k_S34_He4_to_Ar38_reaclib)*Y(S34)*state.rho
              - screened_rates(k_S34_He4_to_n_Ar37_reaclib)*Y(S34)*state.rho
              - screened_rates(k_S34_He4_to_p_Cl37_reaclib)*Y(S34)*state.rho
              + screened_rates(k_Si30_He4_to_S34_reaclib)*Y(Si30)*state.rho
              + screened_rates(k_Si31_He4_to_n_S34_reaclib)*Y(Si31)*state.rho;
    jac.set(S34, He4, scratch);

    jac.set(S34, C12, 0.0);

    jac.set(S34, O16, 0.0);

    jac.set(S34, Ne20, 0.0);

    jac.set(S34, Mg23, 0.0);

    jac.set(S34, Mg24, 0.0);

    jac.set(S34, Al27, 0.0);

    jac.set(S34, Si27, 0.0);

    jac.set(S34, Si28, 0.0);

    jac.set(S34, Si29, 0.0);

    scratch = screened_rates(k_Si30_He4_to_S34_reaclib)*Y(He4)*state.rho;
    jac.set(S34, Si30, scratch);

    scratch = screened_rates(k_Si31_He4_to_n_S34_reaclib)*Y(He4)*state.rho;
    jac.set(S34, Si31, scratch);

    jac.set(S34, P29, 0.0);

    scratch = screened_rates(k_P31_He4_to_p_S34_reaclib)*Y(He4)*state.rho;
    jac.set(S34, P31, scratch);

    scratch = screened_rates(k_P33_p_to_S34_reaclib)*Y(H1)*state.rho;
    jac.set(S34, P33, scratch);

    jac.set(S34, S32, 0.0);

    scratch = screened_rates(k_S33_n_to_S34_reaclib)*Y(N)*state.rho;
    jac.set(S34, S33, scratch);

    scratch = -screened_rates(k_S34_He4_to_Ar38_reaclib)*Y(He4)*state.rho
              - screened_rates(k_S34_He4_to_n_Ar37_reaclib)*Y(He4)*state.rho
              - screened_rates(k_S34_He4_to_p_Cl37_reaclib)*Y(He4)*state.rho
              - screened_rates(k_S34_n_to_He4_Si31_reaclib)*Y(N)*state.rho
              - screened_rates(k_S34_n_to_S35_reaclib)*Y(N)*state.rho
              - screened_rates(k_S34_p_to_Cl35_reaclib)*Y(H1)*state.rho
              - screened_rates(k_S34_p_to_He4_P31_reaclib)*Y(H1)*state.rho
              - screened_rates(k_S34_to_He4_Si30_reaclib)
              - screened_rates(k_S34_to_n_S33_reaclib)
              - screened_rates(k_S34_to_p_P33_reaclib);
    jac.set(S34, S34, scratch);

    scratch = screened_rates(k_S35_to_n_S34_reaclib);
    jac.set(S34, S35, scratch);

    jac.set(S34, Cl33, 0.0);

    scratch = screened_rates(k_Cl35_to_p_S34_reaclib);
    jac.set(S34, Cl35, scratch);

    scratch = screened_rates(k_Cl37_p_to_He4_S34_reaclib)*Y(H1)*state.rho;
    jac.set(S34, Cl37, scratch);

    jac.set(S34, Ar36, 0.0);

    scratch = screened_rates(k_Ar37_n_to_He4_S34_reaclib)*Y(N)*state.rho;
    jac.set(S34, Ar37, scratch);

    scratch = screened_rates(k_Ar38_to_He4_S34_reaclib);
    jac.set(S34, Ar38, scratch);

    jac.set(S34, K37, 0.0);

    jac.set(S34, K39, 0.0);

    jac.set(S34, Ca40, 0.0);

    jac.set(S34, Ca42, 0.0);

    jac.set(S34, Sc41, 0.0);

    jac.set(S34, Sc43, 0.0);

    jac.set(S34, Sc45, 0.0);

    jac.set(S34, Ti44, 0.0);

    jac.set(S34, Ti46, 0.0);

    jac.set(S34, Ti48, 0.0);

    jac.set(S34, V47, 0.0);

    jac.set(S34, V49, 0.0);

    jac.set(S34, V51, 0.0);

    jac.set(S34, Cr48, 0.0);

    jac.set(S34, Cr50, 0.0);

    jac.set(S34, Cr52, 0.0);

    jac.set(S34, Cr54, 0.0);

    jac.set(S34, Mn51, 0.0);

    jac.set(S34, Mn53, 0.0);

    jac.set(S34, Mn55, 0.0);

    jac.set(S34, Fe52, 0.0);

    jac.set(S34, Fe54, 0.0);

    jac.set(S34, Fe56, 0.0);

    jac.set(S34, Fe58, 0.0);

    jac.set(S34, Co55, 0.0);

    jac.set(S34, Co56, 0.0);

    jac.set(S34, Co57, 0.0);

    jac.set(S34, Co59, 0.0);

    jac.set(S34, Ni56, 0.0);

    jac.set(S34, Ni58, 0.0);

    jac.set(S34, Ni60, 0.0);

    jac.set(S34, Cu59, 0.0);

    jac.set(S34, Zn60, 0.0);


    // S35 row

    scratch = screened_rates(k_Ar38_n_to_He4_S35_reaclib)*Y(Ar38)*state.rho
              + screened_rates(k_Cl35_n_to_p_S35_reaclib)*Y(Cl35)*state.rho
              + screened_rates(k_S34_n_to_S35_reaclib)*Y(S34)*state.rho;
    jac.set(S35, N, scratch);

    scratch = -screened_rates(k_S35_p_to_n_Cl35_reaclib)*Y(S35)*state.rho;
    jac.set(S35, H1, scratch);

    scratch = -screened_rates(k_S35_He4_to_n_Ar38_reaclib)*Y(S35)*state.rho
              + screened_rates(k_Si31_He4_to_S35_reaclib)*Y(Si31)*state.rho;
    jac.set(S35, He4, scratch);

    jac.set(S35, C12, 0.0);

    jac.set(S35, O16, 0.0);

    jac.set(S35, Ne20, 0.0);

    jac.set(S35, Mg23, 0.0);

    jac.set(S35, Mg24, 0.0);

    jac.set(S35, Al27, 0.0);

    jac.set(S35, Si27, 0.0);

    jac.set(S35, Si28, 0.0);

    jac.set(S35, Si29, 0.0);

    jac.set(S35, Si30, 0.0);

    scratch = screened_rates(k_Si31_He4_to_S35_reaclib)*Y(He4)*state.rho;
    jac.set(S35, Si31, scratch);

    jac.set(S35, P29, 0.0);

    jac.set(S35, P31, 0.0);

    jac.set(S35, P33, 0.0);

    jac.set(S35, S32, 0.0);

    jac.set(S35, S33, 0.0);

    scratch = screened_rates(k_S34_n_to_S35_reaclib)*Y(N)*state.rho;
    jac.set(S35, S34, scratch);

    scratch = -screened_rates(k_S35_He4_to_n_Ar38_reaclib)*Y(He4)*state.rho
              - screened_rates(k_S35_p_to_n_Cl35_reaclib)*Y(H1)*state.rho
              - screened_rates(k_S35_to_Cl35_beta_neg_weaktab)
              - screened_rates(k_S35_to_He4_Si31_reaclib)
              - screened_rates(k_S35_to_n_S34_reaclib);
    jac.set(S35, S35, scratch);

    jac.set(S35, Cl33, 0.0);

    scratch = screened_rates(k_Cl35_n_to_p_S35_reaclib)*Y(N)*state.rho
              + screened_rates(k_Cl35_to_S35_electron_capture_weaktab);
    jac.set(S35, Cl35, scratch);

    jac.set(S35, Cl37, 0.0);

    jac.set(S35, Ar36, 0.0);

    jac.set(S35, Ar37, 0.0);

    scratch = screened_rates(k_Ar38_n_to_He4_S35_reaclib)*Y(N)*state.rho;
    jac.set(S35, Ar38, scratch);

    jac.set(S35, K37, 0.0);

    jac.set(S35, K39, 0.0);

    jac.set(S35, Ca40, 0.0);

    jac.set(S35, Ca42, 0.0);

    jac.set(S35, Sc41, 0.0);

    jac.set(S35, Sc43, 0.0);

    jac.set(S35, Sc45, 0.0);

    jac.set(S35, Ti44, 0.0);

    jac.set(S35, Ti46, 0.0);

    jac.set(S35, Ti48, 0.0);

    jac.set(S35, V47, 0.0);

    jac.set(S35, V49, 0.0);

    jac.set(S35, V51, 0.0);

    jac.set(S35, Cr48, 0.0);

    jac.set(S35, Cr50, 0.0);

    jac.set(S35, Cr52, 0.0);

    jac.set(S35, Cr54, 0.0);

    jac.set(S35, Mn51, 0.0);

    jac.set(S35, Mn53, 0.0);

    jac.set(S35, Mn55, 0.0);

    jac.set(S35, Fe52, 0.0);

    jac.set(S35, Fe54, 0.0);

    jac.set(S35, Fe56, 0.0);

    jac.set(S35, Fe58, 0.0);

    jac.set(S35, Co55, 0.0);

    jac.set(S35, Co56, 0.0);

    jac.set(S35, Co57, 0.0);

    jac.set(S35, Co59, 0.0);

    jac.set(S35, Ni56, 0.0);

    jac.set(S35, Ni58, 0.0);

    jac.set(S35, Ni60, 0.0);

    jac.set(S35, Cu59, 0.0);

    jac.set(S35, Zn60, 0.0);


    // Cl33 row

    scratch = -screened_rates(k_Cl33_n_to_p_S33_reaclib)*Y(Cl33)*state.rho;
    jac.set(Cl33, N, scratch);

    scratch = screened_rates(k_Ar36_p_to_He4_Cl33_reaclib)*Y(Ar36)*state.rho
              + screened_rates(k_S32_p_to_Cl33_reaclib)*Y(S32)*state.rho
              + screened_rates(k_S33_p_to_n_Cl33_reaclib)*Y(S33)*state.rho;
    jac.set(Cl33, H1, scratch);

    scratch = -screened_rates(k_Cl33_He4_to_K37_reaclib)*Y(Cl33)*state.rho
              - screened_rates(k_Cl33_He4_to_p_Ar36_reaclib)*Y(Cl33)*state.rho
              + screened_rates(k_P29_He4_to_Cl33_reaclib)*Y(P29)*state.rho;
    jac.set(Cl33, He4, scratch);

    jac.set(Cl33, C12, 0.0);

    jac.set(Cl33, O16, 0.0);

    jac.set(Cl33, Ne20, 0.0);

    jac.set(Cl33, Mg23, 0.0);

    jac.set(Cl33, Mg24, 0.0);

    jac.set(Cl33, Al27, 0.0);

    jac.set(Cl33, Si27, 0.0);

    jac.set(Cl33, Si28, 0.0);

    jac.set(Cl33, Si29, 0.0);

    jac.set(Cl33, Si30, 0.0);

    jac.set(Cl33, Si31, 0.0);

    scratch = screened_rates(k_P29_He4_to_Cl33_reaclib)*Y(He4)*state.rho;
    jac.set(Cl33, P29, scratch);

    jac.set(Cl33, P31, 0.0);

    jac.set(Cl33, P33, 0.0);

    scratch = screened_rates(k_S32_p_to_Cl33_reaclib)*Y(H1)*state.rho;
    jac.set(Cl33, S32, scratch);

    scratch = screened_rates(k_S33_p_to_n_Cl33_reaclib)*Y(H1)*state.rho
              + screened_rates(k_S33_to_Cl33_beta_neg_weaktab);
    jac.set(Cl33, S33, scratch);

    jac.set(Cl33, S34, 0.0);

    jac.set(Cl33, S35, 0.0);

    scratch = -screened_rates(k_Cl33_He4_to_K37_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Cl33_He4_to_p_Ar36_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Cl33_n_to_p_S33_reaclib)*Y(N)*state.rho
              - screened_rates(k_Cl33_to_He4_P29_reaclib)
              - screened_rates(k_Cl33_to_S33_electron_capture_weaktab)
              - screened_rates(k_Cl33_to_p_S32_reaclib);
    jac.set(Cl33, Cl33, scratch);

    jac.set(Cl33, Cl35, 0.0);

    jac.set(Cl33, Cl37, 0.0);

    scratch = screened_rates(k_Ar36_p_to_He4_Cl33_reaclib)*Y(H1)*state.rho;
    jac.set(Cl33, Ar36, scratch);

    jac.set(Cl33, Ar37, 0.0);

    jac.set(Cl33, Ar38, 0.0);

    scratch = screened_rates(k_K37_to_He4_Cl33_reaclib);
    jac.set(Cl33, K37, scratch);

    jac.set(Cl33, K39, 0.0);

    jac.set(Cl33, Ca40, 0.0);

    jac.set(Cl33, Ca42, 0.0);

    jac.set(Cl33, Sc41, 0.0);

    jac.set(Cl33, Sc43, 0.0);

    jac.set(Cl33, Sc45, 0.0);

    jac.set(Cl33, Ti44, 0.0);

    jac.set(Cl33, Ti46, 0.0);

    jac.set(Cl33, Ti48, 0.0);

    jac.set(Cl33, V47, 0.0);

    jac.set(Cl33, V49, 0.0);

    jac.set(Cl33, V51, 0.0);

    jac.set(Cl33, Cr48, 0.0);

    jac.set(Cl33, Cr50, 0.0);

    jac.set(Cl33, Cr52, 0.0);

    jac.set(Cl33, Cr54, 0.0);

    jac.set(Cl33, Mn51, 0.0);

    jac.set(Cl33, Mn53, 0.0);

    jac.set(Cl33, Mn55, 0.0);

    jac.set(Cl33, Fe52, 0.0);

    jac.set(Cl33, Fe54, 0.0);

    jac.set(Cl33, Fe56, 0.0);

    jac.set(Cl33, Fe58, 0.0);

    jac.set(Cl33, Co55, 0.0);

    jac.set(Cl33, Co56, 0.0);

    jac.set(Cl33, Co57, 0.0);

    jac.set(Cl33, Co59, 0.0);

    jac.set(Cl33, Ni56, 0.0);

    jac.set(Cl33, Ni58, 0.0);

    jac.set(Cl33, Ni60, 0.0);

    jac.set(Cl33, Cu59, 0.0);

    jac.set(Cl33, Zn60, 0.0);


    // Cl35 row

    scratch = -screened_rates(k_Cl35_n_to_p_S35_reaclib)*Y(Cl35)*state.rho;
    jac.set(Cl35, N, scratch);

    scratch = screened_rates(k_Ar38_p_to_He4_Cl35_reaclib)*Y(Ar38)*state.rho
              - screened_rates(k_Cl35_p_to_Ar36_reaclib)*Y(Cl35)*state.rho
              - screened_rates(k_Cl35_p_to_He4_S32_reaclib)*Y(Cl35)*state.rho
              + screened_rates(k_S34_p_to_Cl35_reaclib)*Y(S34)*state.rho
              + screened_rates(k_S35_p_to_n_Cl35_reaclib)*Y(S35)*state.rho;
    jac.set(Cl35, H1, scratch);

    scratch = -screened_rates(k_Cl35_He4_to_K39_reaclib)*Y(Cl35)*state.rho
              - screened_rates(k_Cl35_He4_to_p_Ar38_reaclib)*Y(Cl35)*state.rho
              + screened_rates(k_P31_He4_to_Cl35_reaclib)*Y(P31)*state.rho
              + screened_rates(k_S32_He4_to_p_Cl35_reaclib)*Y(S32)*state.rho;
    jac.set(Cl35, He4, scratch);

    jac.set(Cl35, C12, 0.0);

    jac.set(Cl35, O16, 0.0);

    jac.set(Cl35, Ne20, 0.0);

    jac.set(Cl35, Mg23, 0.0);

    jac.set(Cl35, Mg24, 0.0);

    jac.set(Cl35, Al27, 0.0);

    jac.set(Cl35, Si27, 0.0);

    jac.set(Cl35, Si28, 0.0);

    jac.set(Cl35, Si29, 0.0);

    jac.set(Cl35, Si30, 0.0);

    jac.set(Cl35, Si31, 0.0);

    jac.set(Cl35, P29, 0.0);

    scratch = screened_rates(k_P31_He4_to_Cl35_reaclib)*Y(He4)*state.rho;
    jac.set(Cl35, P31, scratch);

    jac.set(Cl35, P33, 0.0);

    scratch = screened_rates(k_S32_He4_to_p_Cl35_reaclib)*Y(He4)*state.rho;
    jac.set(Cl35, S32, scratch);

    jac.set(Cl35, S33, 0.0);

    scratch = screened_rates(k_S34_p_to_Cl35_reaclib)*Y(H1)*state.rho;
    jac.set(Cl35, S34, scratch);

    scratch = screened_rates(k_S35_p_to_n_Cl35_reaclib)*Y(H1)*state.rho
              + screened_rates(k_S35_to_Cl35_beta_neg_weaktab);
    jac.set(Cl35, S35, scratch);

    jac.set(Cl35, Cl33, 0.0);

    scratch = -screened_rates(k_Cl35_He4_to_K39_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Cl35_He4_to_p_Ar38_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Cl35_n_to_p_S35_reaclib)*Y(N)*state.rho
              - screened_rates(k_Cl35_p_to_Ar36_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Cl35_p_to_He4_S32_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Cl35_to_He4_P31_reaclib)
              - screened_rates(k_Cl35_to_S35_electron_capture_weaktab)
              - screened_rates(k_Cl35_to_p_S34_reaclib);
    jac.set(Cl35, Cl35, scratch);

    jac.set(Cl35, Cl37, 0.0);

    scratch = screened_rates(k_Ar36_to_p_Cl35_reaclib);
    jac.set(Cl35, Ar36, scratch);

    jac.set(Cl35, Ar37, 0.0);

    scratch = screened_rates(k_Ar38_p_to_He4_Cl35_reaclib)*Y(H1)*state.rho;
    jac.set(Cl35, Ar38, scratch);

    jac.set(Cl35, K37, 0.0);

    scratch = screened_rates(k_K39_to_He4_Cl35_reaclib);
    jac.set(Cl35, K39, scratch);

    jac.set(Cl35, Ca40, 0.0);

    jac.set(Cl35, Ca42, 0.0);

    jac.set(Cl35, Sc41, 0.0);

    jac.set(Cl35, Sc43, 0.0);

    jac.set(Cl35, Sc45, 0.0);

    jac.set(Cl35, Ti44, 0.0);

    jac.set(Cl35, Ti46, 0.0);

    jac.set(Cl35, Ti48, 0.0);

    jac.set(Cl35, V47, 0.0);

    jac.set(Cl35, V49, 0.0);

    jac.set(Cl35, V51, 0.0);

    jac.set(Cl35, Cr48, 0.0);

    jac.set(Cl35, Cr50, 0.0);

    jac.set(Cl35, Cr52, 0.0);

    jac.set(Cl35, Cr54, 0.0);

    jac.set(Cl35, Mn51, 0.0);

    jac.set(Cl35, Mn53, 0.0);

    jac.set(Cl35, Mn55, 0.0);

    jac.set(Cl35, Fe52, 0.0);

    jac.set(Cl35, Fe54, 0.0);

    jac.set(Cl35, Fe56, 0.0);

    jac.set(Cl35, Fe58, 0.0);

    jac.set(Cl35, Co55, 0.0);

    jac.set(Cl35, Co56, 0.0);

    jac.set(Cl35, Co57, 0.0);

    jac.set(Cl35, Co59, 0.0);

    jac.set(Cl35, Ni56, 0.0);

    jac.set(Cl35, Ni58, 0.0);

    jac.set(Cl35, Ni60, 0.0);

    jac.set(Cl35, Cu59, 0.0);

    jac.set(Cl35, Zn60, 0.0);


    // Cl37 row

    scratch = screened_rates(k_Ar37_n_to_p_Cl37_reaclib)*Y(Ar37)*state.rho;
    jac.set(Cl37, N, scratch);

    scratch = -screened_rates(k_Cl37_p_to_Ar38_reaclib)*Y(Cl37)*state.rho
              - screened_rates(k_Cl37_p_to_He4_S34_reaclib)*Y(Cl37)*state.rho
              - screened_rates(k_Cl37_p_to_n_Ar37_reaclib)*Y(Cl37)*state.rho;
    jac.set(Cl37, H1, scratch);

    scratch = screened_rates(k_P33_He4_to_Cl37_reaclib)*Y(P33)*state.rho
              + screened_rates(k_S34_He4_to_p_Cl37_reaclib)*Y(S34)*state.rho;
    jac.set(Cl37, He4, scratch);

    jac.set(Cl37, C12, 0.0);

    jac.set(Cl37, O16, 0.0);

    jac.set(Cl37, Ne20, 0.0);

    jac.set(Cl37, Mg23, 0.0);

    jac.set(Cl37, Mg24, 0.0);

    jac.set(Cl37, Al27, 0.0);

    jac.set(Cl37, Si27, 0.0);

    jac.set(Cl37, Si28, 0.0);

    jac.set(Cl37, Si29, 0.0);

    jac.set(Cl37, Si30, 0.0);

    jac.set(Cl37, Si31, 0.0);

    jac.set(Cl37, P29, 0.0);

    jac.set(Cl37, P31, 0.0);

    scratch = screened_rates(k_P33_He4_to_Cl37_reaclib)*Y(He4)*state.rho;
    jac.set(Cl37, P33, scratch);

    jac.set(Cl37, S32, 0.0);

    jac.set(Cl37, S33, 0.0);

    scratch = screened_rates(k_S34_He4_to_p_Cl37_reaclib)*Y(He4)*state.rho;
    jac.set(Cl37, S34, scratch);

    jac.set(Cl37, S35, 0.0);

    jac.set(Cl37, Cl33, 0.0);

    jac.set(Cl37, Cl35, 0.0);

    scratch = -screened_rates(k_Cl37_p_to_Ar38_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Cl37_p_to_He4_S34_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Cl37_p_to_n_Ar37_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Cl37_to_Ar37_beta_neg_weaktab)
              - screened_rates(k_Cl37_to_He4_P33_reaclib);
    jac.set(Cl37, Cl37, scratch);

    jac.set(Cl37, Ar36, 0.0);

    scratch = screened_rates(k_Ar37_n_to_p_Cl37_reaclib)*Y(N)*state.rho
              + screened_rates(k_Ar37_to_Cl37_electron_capture_weaktab);
    jac.set(Cl37, Ar37, scratch);

    scratch = screened_rates(k_Ar38_to_p_Cl37_reaclib);
    jac.set(Cl37, Ar38, scratch);

    jac.set(Cl37, K37, 0.0);

    jac.set(Cl37, K39, 0.0);

    jac.set(Cl37, Ca40, 0.0);

    jac.set(Cl37, Ca42, 0.0);

    jac.set(Cl37, Sc41, 0.0);

    jac.set(Cl37, Sc43, 0.0);

    jac.set(Cl37, Sc45, 0.0);

    jac.set(Cl37, Ti44, 0.0);

    jac.set(Cl37, Ti46, 0.0);

    jac.set(Cl37, Ti48, 0.0);

    jac.set(Cl37, V47, 0.0);

    jac.set(Cl37, V49, 0.0);

    jac.set(Cl37, V51, 0.0);

    jac.set(Cl37, Cr48, 0.0);

    jac.set(Cl37, Cr50, 0.0);

    jac.set(Cl37, Cr52, 0.0);

    jac.set(Cl37, Cr54, 0.0);

    jac.set(Cl37, Mn51, 0.0);

    jac.set(Cl37, Mn53, 0.0);

    jac.set(Cl37, Mn55, 0.0);

    jac.set(Cl37, Fe52, 0.0);

    jac.set(Cl37, Fe54, 0.0);

    jac.set(Cl37, Fe56, 0.0);

    jac.set(Cl37, Fe58, 0.0);

    jac.set(Cl37, Co55, 0.0);

    jac.set(Cl37, Co56, 0.0);

    jac.set(Cl37, Co57, 0.0);

    jac.set(Cl37, Co59, 0.0);

    jac.set(Cl37, Ni56, 0.0);

    jac.set(Cl37, Ni58, 0.0);

    jac.set(Cl37, Ni60, 0.0);

    jac.set(Cl37, Cu59, 0.0);

    jac.set(Cl37, Zn60, 0.0);


    // Ar36 row

    scratch = -screened_rates(k_Ar36_n_to_Ar37_reaclib)*Y(Ar36)*state.rho
              - screened_rates(k_Ar36_n_to_He4_S33_reaclib)*Y(Ar36)*state.rho;
    jac.set(Ar36, N, scratch);

    scratch = -screened_rates(k_Ar36_p_to_He4_Cl33_reaclib)*Y(Ar36)*state.rho
              - screened_rates(k_Ar36_p_to_K37_reaclib)*Y(Ar36)*state.rho
              + screened_rates(k_Cl35_p_to_Ar36_reaclib)*Y(Cl35)*state.rho
              + screened_rates(k_K39_p_to_He4_Ar36_reaclib)*Y(K39)*state.rho;
    jac.set(Ar36, H1, scratch);

    scratch = -screened_rates(k_Ar36_He4_to_Ca40_reaclib)*Y(Ar36)*state.rho
              - screened_rates(k_Ar36_He4_to_p_K39_reaclib)*Y(Ar36)*state.rho
              + screened_rates(k_Cl33_He4_to_p_Ar36_reaclib)*Y(Cl33)*state.rho
              + screened_rates(k_S32_He4_to_Ar36_reaclib)*Y(S32)*state.rho
              + screened_rates(k_S33_He4_to_n_Ar36_reaclib)*Y(S33)*state.rho;
    jac.set(Ar36, He4, scratch);

    jac.set(Ar36, C12, 0.0);

    jac.set(Ar36, O16, 0.0);

    jac.set(Ar36, Ne20, 0.0);

    jac.set(Ar36, Mg23, 0.0);

    jac.set(Ar36, Mg24, 0.0);

    jac.set(Ar36, Al27, 0.0);

    jac.set(Ar36, Si27, 0.0);

    jac.set(Ar36, Si28, 0.0);

    jac.set(Ar36, Si29, 0.0);

    jac.set(Ar36, Si30, 0.0);

    jac.set(Ar36, Si31, 0.0);

    jac.set(Ar36, P29, 0.0);

    jac.set(Ar36, P31, 0.0);

    jac.set(Ar36, P33, 0.0);

    scratch = screened_rates(k_S32_He4_to_Ar36_reaclib)*Y(He4)*state.rho;
    jac.set(Ar36, S32, scratch);

    scratch = screened_rates(k_S33_He4_to_n_Ar36_reaclib)*Y(He4)*state.rho;
    jac.set(Ar36, S33, scratch);

    jac.set(Ar36, S34, 0.0);

    jac.set(Ar36, S35, 0.0);

    scratch = screened_rates(k_Cl33_He4_to_p_Ar36_reaclib)*Y(He4)*state.rho;
    jac.set(Ar36, Cl33, scratch);

    scratch = screened_rates(k_Cl35_p_to_Ar36_reaclib)*Y(H1)*state.rho;
    jac.set(Ar36, Cl35, scratch);

    jac.set(Ar36, Cl37, 0.0);

    scratch = -screened_rates(k_Ar36_He4_to_Ca40_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ar36_He4_to_p_K39_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ar36_n_to_Ar37_reaclib)*Y(N)*state.rho
              - screened_rates(k_Ar36_n_to_He4_S33_reaclib)*Y(N)*state.rho
              - screened_rates(k_Ar36_p_to_He4_Cl33_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Ar36_p_to_K37_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Ar36_to_He4_S32_reaclib)
              - screened_rates(k_Ar36_to_p_Cl35_reaclib);
    jac.set(Ar36, Ar36, scratch);

    scratch = screened_rates(k_Ar37_to_n_Ar36_reaclib);
    jac.set(Ar36, Ar37, scratch);

    jac.set(Ar36, Ar38, 0.0);

    scratch = screened_rates(k_K37_to_p_Ar36_reaclib);
    jac.set(Ar36, K37, scratch);

    scratch = screened_rates(k_K39_p_to_He4_Ar36_reaclib)*Y(H1)*state.rho;
    jac.set(Ar36, K39, scratch);

    scratch = screened_rates(k_Ca40_to_He4_Ar36_reaclib);
    jac.set(Ar36, Ca40, scratch);

    jac.set(Ar36, Ca42, 0.0);

    jac.set(Ar36, Sc41, 0.0);

    jac.set(Ar36, Sc43, 0.0);

    jac.set(Ar36, Sc45, 0.0);

    jac.set(Ar36, Ti44, 0.0);

    jac.set(Ar36, Ti46, 0.0);

    jac.set(Ar36, Ti48, 0.0);

    jac.set(Ar36, V47, 0.0);

    jac.set(Ar36, V49, 0.0);

    jac.set(Ar36, V51, 0.0);

    jac.set(Ar36, Cr48, 0.0);

    jac.set(Ar36, Cr50, 0.0);

    jac.set(Ar36, Cr52, 0.0);

    jac.set(Ar36, Cr54, 0.0);

    jac.set(Ar36, Mn51, 0.0);

    jac.set(Ar36, Mn53, 0.0);

    jac.set(Ar36, Mn55, 0.0);

    jac.set(Ar36, Fe52, 0.0);

    jac.set(Ar36, Fe54, 0.0);

    jac.set(Ar36, Fe56, 0.0);

    jac.set(Ar36, Fe58, 0.0);

    jac.set(Ar36, Co55, 0.0);

    jac.set(Ar36, Co56, 0.0);

    jac.set(Ar36, Co57, 0.0);

    jac.set(Ar36, Co59, 0.0);

    jac.set(Ar36, Ni56, 0.0);

    jac.set(Ar36, Ni58, 0.0);

    jac.set(Ar36, Ni60, 0.0);

    jac.set(Ar36, Cu59, 0.0);

    jac.set(Ar36, Zn60, 0.0);


    // Ar37 row

    scratch = screened_rates(k_Ar36_n_to_Ar37_reaclib)*Y(Ar36)*state.rho
              - screened_rates(k_Ar37_n_to_Ar38_reaclib)*Y(Ar37)*state.rho
              - screened_rates(k_Ar37_n_to_He4_S34_reaclib)*Y(Ar37)*state.rho
              - screened_rates(k_Ar37_n_to_p_Cl37_reaclib)*Y(Ar37)*state.rho
              + screened_rates(k_Ca40_n_to_He4_Ar37_reaclib)*Y(Ca40)*state.rho
              + screened_rates(k_K37_n_to_p_Ar37_reaclib)*Y(K37)*state.rho;
    jac.set(Ar37, N, scratch);

    scratch = -screened_rates(k_Ar37_p_to_n_K37_reaclib)*Y(Ar37)*state.rho
              + screened_rates(k_Cl37_p_to_n_Ar37_reaclib)*Y(Cl37)*state.rho;
    jac.set(Ar37, H1, scratch);

    scratch = -screened_rates(k_Ar37_He4_to_n_Ca40_reaclib)*Y(Ar37)*state.rho
              + screened_rates(k_S33_He4_to_Ar37_reaclib)*Y(S33)*state.rho
              + screened_rates(k_S34_He4_to_n_Ar37_reaclib)*Y(S34)*state.rho;
    jac.set(Ar37, He4, scratch);

    jac.set(Ar37, C12, 0.0);

    jac.set(Ar37, O16, 0.0);

    jac.set(Ar37, Ne20, 0.0);

    jac.set(Ar37, Mg23, 0.0);

    jac.set(Ar37, Mg24, 0.0);

    jac.set(Ar37, Al27, 0.0);

    jac.set(Ar37, Si27, 0.0);

    jac.set(Ar37, Si28, 0.0);

    jac.set(Ar37, Si29, 0.0);

    jac.set(Ar37, Si30, 0.0);

    jac.set(Ar37, Si31, 0.0);

    jac.set(Ar37, P29, 0.0);

    jac.set(Ar37, P31, 0.0);

    jac.set(Ar37, P33, 0.0);

    jac.set(Ar37, S32, 0.0);

    scratch = screened_rates(k_S33_He4_to_Ar37_reaclib)*Y(He4)*state.rho;
    jac.set(Ar37, S33, scratch);

    scratch = screened_rates(k_S34_He4_to_n_Ar37_reaclib)*Y(He4)*state.rho;
    jac.set(Ar37, S34, scratch);

    jac.set(Ar37, S35, 0.0);

    jac.set(Ar37, Cl33, 0.0);

    jac.set(Ar37, Cl35, 0.0);

    scratch = screened_rates(k_Cl37_p_to_n_Ar37_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Cl37_to_Ar37_beta_neg_weaktab);
    jac.set(Ar37, Cl37, scratch);

    scratch = screened_rates(k_Ar36_n_to_Ar37_reaclib)*Y(N)*state.rho;
    jac.set(Ar37, Ar36, scratch);

    scratch = -screened_rates(k_Ar37_He4_to_n_Ca40_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ar37_n_to_Ar38_reaclib)*Y(N)*state.rho
              - screened_rates(k_Ar37_n_to_He4_S34_reaclib)*Y(N)*state.rho
              - screened_rates(k_Ar37_n_to_p_Cl37_reaclib)*Y(N)*state.rho
              - screened_rates(k_Ar37_p_to_n_K37_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Ar37_to_Cl37_electron_capture_weaktab)
              - screened_rates(k_Ar37_to_He4_S33_reaclib)
              - screened_rates(k_Ar37_to_K37_beta_neg_weaktab)
              - screened_rates(k_Ar37_to_n_Ar36_reaclib);
    jac.set(Ar37, Ar37, scratch);

    scratch = screened_rates(k_Ar38_to_n_Ar37_reaclib);
    jac.set(Ar37, Ar38, scratch);

    scratch = screened_rates(k_K37_n_to_p_Ar37_reaclib)*Y(N)*state.rho
              + screened_rates(k_K37_to_Ar37_electron_capture_weaktab);
    jac.set(Ar37, K37, scratch);

    jac.set(Ar37, K39, 0.0);

    scratch = screened_rates(k_Ca40_n_to_He4_Ar37_reaclib)*Y(N)*state.rho;
    jac.set(Ar37, Ca40, scratch);

    jac.set(Ar37, Ca42, 0.0);

    jac.set(Ar37, Sc41, 0.0);

    jac.set(Ar37, Sc43, 0.0);

    jac.set(Ar37, Sc45, 0.0);

    jac.set(Ar37, Ti44, 0.0);

    jac.set(Ar37, Ti46, 0.0);

    jac.set(Ar37, Ti48, 0.0);

    jac.set(Ar37, V47, 0.0);

    jac.set(Ar37, V49, 0.0);

    jac.set(Ar37, V51, 0.0);

    jac.set(Ar37, Cr48, 0.0);

    jac.set(Ar37, Cr50, 0.0);

    jac.set(Ar37, Cr52, 0.0);

    jac.set(Ar37, Cr54, 0.0);

    jac.set(Ar37, Mn51, 0.0);

    jac.set(Ar37, Mn53, 0.0);

    jac.set(Ar37, Mn55, 0.0);

    jac.set(Ar37, Fe52, 0.0);

    jac.set(Ar37, Fe54, 0.0);

    jac.set(Ar37, Fe56, 0.0);

    jac.set(Ar37, Fe58, 0.0);

    jac.set(Ar37, Co55, 0.0);

    jac.set(Ar37, Co56, 0.0);

    jac.set(Ar37, Co57, 0.0);

    jac.set(Ar37, Co59, 0.0);

    jac.set(Ar37, Ni56, 0.0);

    jac.set(Ar37, Ni58, 0.0);

    jac.set(Ar37, Ni60, 0.0);

    jac.set(Ar37, Cu59, 0.0);

    jac.set(Ar37, Zn60, 0.0);


    // Ar38 row

    scratch = screened_rates(k_Ar37_n_to_Ar38_reaclib)*Y(Ar37)*state.rho
              - screened_rates(k_Ar38_n_to_He4_S35_reaclib)*Y(Ar38)*state.rho;
    jac.set(Ar38, N, scratch);

    scratch = -screened_rates(k_Ar38_p_to_He4_Cl35_reaclib)*Y(Ar38)*state.rho
              - screened_rates(k_Ar38_p_to_K39_reaclib)*Y(Ar38)*state.rho
              + screened_rates(k_Cl37_p_to_Ar38_reaclib)*Y(Cl37)*state.rho;
    jac.set(Ar38, H1, scratch);

    scratch = -screened_rates(k_Ar38_He4_to_Ca42_reaclib)*Y(Ar38)*state.rho
              + screened_rates(k_Cl35_He4_to_p_Ar38_reaclib)*Y(Cl35)*state.rho
              + screened_rates(k_S34_He4_to_Ar38_reaclib)*Y(S34)*state.rho
              + screened_rates(k_S35_He4_to_n_Ar38_reaclib)*Y(S35)*state.rho;
    jac.set(Ar38, He4, scratch);

    jac.set(Ar38, C12, 0.0);

    jac.set(Ar38, O16, 0.0);

    jac.set(Ar38, Ne20, 0.0);

    jac.set(Ar38, Mg23, 0.0);

    jac.set(Ar38, Mg24, 0.0);

    jac.set(Ar38, Al27, 0.0);

    jac.set(Ar38, Si27, 0.0);

    jac.set(Ar38, Si28, 0.0);

    jac.set(Ar38, Si29, 0.0);

    jac.set(Ar38, Si30, 0.0);

    jac.set(Ar38, Si31, 0.0);

    jac.set(Ar38, P29, 0.0);

    jac.set(Ar38, P31, 0.0);

    jac.set(Ar38, P33, 0.0);

    jac.set(Ar38, S32, 0.0);

    jac.set(Ar38, S33, 0.0);

    scratch = screened_rates(k_S34_He4_to_Ar38_reaclib)*Y(He4)*state.rho;
    jac.set(Ar38, S34, scratch);

    scratch = screened_rates(k_S35_He4_to_n_Ar38_reaclib)*Y(He4)*state.rho;
    jac.set(Ar38, S35, scratch);

    jac.set(Ar38, Cl33, 0.0);

    scratch = screened_rates(k_Cl35_He4_to_p_Ar38_reaclib)*Y(He4)*state.rho;
    jac.set(Ar38, Cl35, scratch);

    scratch = screened_rates(k_Cl37_p_to_Ar38_reaclib)*Y(H1)*state.rho;
    jac.set(Ar38, Cl37, scratch);

    jac.set(Ar38, Ar36, 0.0);

    scratch = screened_rates(k_Ar37_n_to_Ar38_reaclib)*Y(N)*state.rho;
    jac.set(Ar38, Ar37, scratch);

    scratch = -screened_rates(k_Ar38_He4_to_Ca42_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ar38_n_to_He4_S35_reaclib)*Y(N)*state.rho
              - screened_rates(k_Ar38_p_to_He4_Cl35_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Ar38_p_to_K39_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Ar38_to_He4_S34_reaclib)
              - screened_rates(k_Ar38_to_n_Ar37_reaclib)
              - screened_rates(k_Ar38_to_p_Cl37_reaclib);
    jac.set(Ar38, Ar38, scratch);

    jac.set(Ar38, K37, 0.0);

    scratch = screened_rates(k_K39_to_p_Ar38_reaclib);
    jac.set(Ar38, K39, scratch);

    jac.set(Ar38, Ca40, 0.0);

    scratch = screened_rates(k_Ca42_to_He4_Ar38_reaclib);
    jac.set(Ar38, Ca42, scratch);

    jac.set(Ar38, Sc41, 0.0);

    jac.set(Ar38, Sc43, 0.0);

    jac.set(Ar38, Sc45, 0.0);

    jac.set(Ar38, Ti44, 0.0);

    jac.set(Ar38, Ti46, 0.0);

    jac.set(Ar38, Ti48, 0.0);

    jac.set(Ar38, V47, 0.0);

    jac.set(Ar38, V49, 0.0);

    jac.set(Ar38, V51, 0.0);

    jac.set(Ar38, Cr48, 0.0);

    jac.set(Ar38, Cr50, 0.0);

    jac.set(Ar38, Cr52, 0.0);

    jac.set(Ar38, Cr54, 0.0);

    jac.set(Ar38, Mn51, 0.0);

    jac.set(Ar38, Mn53, 0.0);

    jac.set(Ar38, Mn55, 0.0);

    jac.set(Ar38, Fe52, 0.0);

    jac.set(Ar38, Fe54, 0.0);

    jac.set(Ar38, Fe56, 0.0);

    jac.set(Ar38, Fe58, 0.0);

    jac.set(Ar38, Co55, 0.0);

    jac.set(Ar38, Co56, 0.0);

    jac.set(Ar38, Co57, 0.0);

    jac.set(Ar38, Co59, 0.0);

    jac.set(Ar38, Ni56, 0.0);

    jac.set(Ar38, Ni58, 0.0);

    jac.set(Ar38, Ni60, 0.0);

    jac.set(Ar38, Cu59, 0.0);

    jac.set(Ar38, Zn60, 0.0);


    // K37 row

    scratch = -screened_rates(k_K37_n_to_p_Ar37_reaclib)*Y(K37)*state.rho;
    jac.set(K37, N, scratch);

    scratch = screened_rates(k_Ar36_p_to_K37_reaclib)*Y(Ar36)*state.rho
              + screened_rates(k_Ar37_p_to_n_K37_reaclib)*Y(Ar37)*state.rho
              + screened_rates(k_Ca40_p_to_He4_K37_reaclib)*Y(Ca40)*state.rho;
    jac.set(K37, H1, scratch);

    scratch = screened_rates(k_Cl33_He4_to_K37_reaclib)*Y(Cl33)*state.rho
              - screened_rates(k_K37_He4_to_Sc41_reaclib)*Y(K37)*state.rho
              - screened_rates(k_K37_He4_to_p_Ca40_reaclib)*Y(K37)*state.rho;
    jac.set(K37, He4, scratch);

    jac.set(K37, C12, 0.0);

    jac.set(K37, O16, 0.0);

    jac.set(K37, Ne20, 0.0);

    jac.set(K37, Mg23, 0.0);

    jac.set(K37, Mg24, 0.0);

    jac.set(K37, Al27, 0.0);

    jac.set(K37, Si27, 0.0);

    jac.set(K37, Si28, 0.0);

    jac.set(K37, Si29, 0.0);

    jac.set(K37, Si30, 0.0);

    jac.set(K37, Si31, 0.0);

    jac.set(K37, P29, 0.0);

    jac.set(K37, P31, 0.0);

    jac.set(K37, P33, 0.0);

    jac.set(K37, S32, 0.0);

    jac.set(K37, S33, 0.0);

    jac.set(K37, S34, 0.0);

    jac.set(K37, S35, 0.0);

    scratch = screened_rates(k_Cl33_He4_to_K37_reaclib)*Y(He4)*state.rho;
    jac.set(K37, Cl33, scratch);

    jac.set(K37, Cl35, 0.0);

    jac.set(K37, Cl37, 0.0);

    scratch = screened_rates(k_Ar36_p_to_K37_reaclib)*Y(H1)*state.rho;
    jac.set(K37, Ar36, scratch);

    scratch = screened_rates(k_Ar37_p_to_n_K37_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Ar37_to_K37_beta_neg_weaktab);
    jac.set(K37, Ar37, scratch);

    jac.set(K37, Ar38, 0.0);

    scratch = -screened_rates(k_K37_He4_to_Sc41_reaclib)*Y(He4)*state.rho
              - screened_rates(k_K37_He4_to_p_Ca40_reaclib)*Y(He4)*state.rho
              - screened_rates(k_K37_n_to_p_Ar37_reaclib)*Y(N)*state.rho
              - screened_rates(k_K37_to_Ar37_electron_capture_weaktab)
              - screened_rates(k_K37_to_He4_Cl33_reaclib)
              - screened_rates(k_K37_to_p_Ar36_reaclib);
    jac.set(K37, K37, scratch);

    jac.set(K37, K39, 0.0);

    scratch = screened_rates(k_Ca40_p_to_He4_K37_reaclib)*Y(H1)*state.rho;
    jac.set(K37, Ca40, scratch);

    jac.set(K37, Ca42, 0.0);

    scratch = screened_rates(k_Sc41_to_He4_K37_reaclib);
    jac.set(K37, Sc41, scratch);

    jac.set(K37, Sc43, 0.0);

    jac.set(K37, Sc45, 0.0);

    jac.set(K37, Ti44, 0.0);

    jac.set(K37, Ti46, 0.0);

    jac.set(K37, Ti48, 0.0);

    jac.set(K37, V47, 0.0);

    jac.set(K37, V49, 0.0);

    jac.set(K37, V51, 0.0);

    jac.set(K37, Cr48, 0.0);

    jac.set(K37, Cr50, 0.0);

    jac.set(K37, Cr52, 0.0);

    jac.set(K37, Cr54, 0.0);

    jac.set(K37, Mn51, 0.0);

    jac.set(K37, Mn53, 0.0);

    jac.set(K37, Mn55, 0.0);

    jac.set(K37, Fe52, 0.0);

    jac.set(K37, Fe54, 0.0);

    jac.set(K37, Fe56, 0.0);

    jac.set(K37, Fe58, 0.0);

    jac.set(K37, Co55, 0.0);

    jac.set(K37, Co56, 0.0);

    jac.set(K37, Co57, 0.0);

    jac.set(K37, Co59, 0.0);

    jac.set(K37, Ni56, 0.0);

    jac.set(K37, Ni58, 0.0);

    jac.set(K37, Ni60, 0.0);

    jac.set(K37, Cu59, 0.0);

    jac.set(K37, Zn60, 0.0);


    // K39 row

    jac.set(K39, N, 0.0);

    scratch = screened_rates(k_Ar38_p_to_K39_reaclib)*Y(Ar38)*state.rho
              + screened_rates(k_Ca42_p_to_He4_K39_reaclib)*Y(Ca42)*state.rho
              - screened_rates(k_K39_p_to_Ca40_reaclib)*Y(K39)*state.rho
              - screened_rates(k_K39_p_to_He4_Ar36_reaclib)*Y(K39)*state.rho;
    jac.set(K39, H1, scratch);

    scratch = screened_rates(k_Ar36_He4_to_p_K39_reaclib)*Y(Ar36)*state.rho
              + screened_rates(k_Cl35_He4_to_K39_reaclib)*Y(Cl35)*state.rho
              - screened_rates(k_K39_He4_to_Sc43_reaclib)*Y(K39)*state.rho
              - screened_rates(k_K39_He4_to_p_Ca42_reaclib)*Y(K39)*state.rho;
    jac.set(K39, He4, scratch);

    jac.set(K39, C12, 0.0);

    jac.set(K39, O16, 0.0);

    jac.set(K39, Ne20, 0.0);

    jac.set(K39, Mg23, 0.0);

    jac.set(K39, Mg24, 0.0);

    jac.set(K39, Al27, 0.0);

    jac.set(K39, Si27, 0.0);

    jac.set(K39, Si28, 0.0);

    jac.set(K39, Si29, 0.0);

    jac.set(K39, Si30, 0.0);

    jac.set(K39, Si31, 0.0);

    jac.set(K39, P29, 0.0);

    jac.set(K39, P31, 0.0);

    jac.set(K39, P33, 0.0);

    jac.set(K39, S32, 0.0);

    jac.set(K39, S33, 0.0);

    jac.set(K39, S34, 0.0);

    jac.set(K39, S35, 0.0);

    jac.set(K39, Cl33, 0.0);

    scratch = screened_rates(k_Cl35_He4_to_K39_reaclib)*Y(He4)*state.rho;
    jac.set(K39, Cl35, scratch);

    jac.set(K39, Cl37, 0.0);

    scratch = screened_rates(k_Ar36_He4_to_p_K39_reaclib)*Y(He4)*state.rho;
    jac.set(K39, Ar36, scratch);

    jac.set(K39, Ar37, 0.0);

    scratch = screened_rates(k_Ar38_p_to_K39_reaclib)*Y(H1)*state.rho;
    jac.set(K39, Ar38, scratch);

    jac.set(K39, K37, 0.0);

    scratch = -screened_rates(k_K39_He4_to_Sc43_reaclib)*Y(He4)*state.rho
              - screened_rates(k_K39_He4_to_p_Ca42_reaclib)*Y(He4)*state.rho
              - screened_rates(k_K39_p_to_Ca40_reaclib)*Y(H1)*state.rho
              - screened_rates(k_K39_p_to_He4_Ar36_reaclib)*Y(H1)*state.rho
              - screened_rates(k_K39_to_He4_Cl35_reaclib)
              - screened_rates(k_K39_to_p_Ar38_reaclib);
    jac.set(K39, K39, scratch);

    scratch = screened_rates(k_Ca40_to_p_K39_reaclib);
    jac.set(K39, Ca40, scratch);

    scratch = screened_rates(k_Ca42_p_to_He4_K39_reaclib)*Y(H1)*state.rho;
    jac.set(K39, Ca42, scratch);

    jac.set(K39, Sc41, 0.0);

    scratch = screened_rates(k_Sc43_to_He4_K39_reaclib);
    jac.set(K39, Sc43, scratch);

    jac.set(K39, Sc45, 0.0);

    jac.set(K39, Ti44, 0.0);

    jac.set(K39, Ti46, 0.0);

    jac.set(K39, Ti48, 0.0);

    jac.set(K39, V47, 0.0);

    jac.set(K39, V49, 0.0);

    jac.set(K39, V51, 0.0);

    jac.set(K39, Cr48, 0.0);

    jac.set(K39, Cr50, 0.0);

    jac.set(K39, Cr52, 0.0);

    jac.set(K39, Cr54, 0.0);

    jac.set(K39, Mn51, 0.0);

    jac.set(K39, Mn53, 0.0);

    jac.set(K39, Mn55, 0.0);

    jac.set(K39, Fe52, 0.0);

    jac.set(K39, Fe54, 0.0);

    jac.set(K39, Fe56, 0.0);

    jac.set(K39, Fe58, 0.0);

    jac.set(K39, Co55, 0.0);

    jac.set(K39, Co56, 0.0);

    jac.set(K39, Co57, 0.0);

    jac.set(K39, Co59, 0.0);

    jac.set(K39, Ni56, 0.0);

    jac.set(K39, Ni58, 0.0);

    jac.set(K39, Ni60, 0.0);

    jac.set(K39, Cu59, 0.0);

    jac.set(K39, Zn60, 0.0);


    // Ca40 row

    scratch = -1.0*screened_rates(k_Ca40_n_n_to_Ca42_approx)*Y(Ca40)*Y(N)*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_Ca40_n_to_He4_Ar37_reaclib)*Y(Ca40)*state.rho
              - 0.5*Y(Ca40)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Ca40_n_n_to_Ca42_approx_dYN
              + Y(Ca42)*rate_eval.drate_Ca42_to_Ca40_n_n_approx_dYN;
    jac.set(Ca40, N, scratch);

    scratch = -screened_rates(k_Ca40_p_to_He4_K37_reaclib)*Y(Ca40)*state.rho
              - screened_rates(k_Ca40_p_to_Sc41_reaclib)*Y(Ca40)*state.rho
              + screened_rates(k_K39_p_to_Ca40_reaclib)*Y(K39)*state.rho
              + screened_rates(k_Sc43_p_to_He4_Ca40_reaclib)*Y(Sc43)*state.rho;
    jac.set(Ca40, H1, scratch);

    scratch = screened_rates(k_Ar36_He4_to_Ca40_reaclib)*Y(Ar36)*state.rho
              + screened_rates(k_Ar37_He4_to_n_Ca40_reaclib)*Y(Ar37)*state.rho
              - screened_rates(k_Ca40_He4_to_Ti44_reaclib)*Y(Ca40)*state.rho
              - screened_rates(k_Ca40_He4_to_p_Sc43_reaclib)*Y(Ca40)*state.rho
              + screened_rates(k_K37_He4_to_p_Ca40_reaclib)*Y(K37)*state.rho;
    jac.set(Ca40, He4, scratch);

    jac.set(Ca40, C12, 0.0);

    jac.set(Ca40, O16, 0.0);

    jac.set(Ca40, Ne20, 0.0);

    jac.set(Ca40, Mg23, 0.0);

    jac.set(Ca40, Mg24, 0.0);

    jac.set(Ca40, Al27, 0.0);

    jac.set(Ca40, Si27, 0.0);

    jac.set(Ca40, Si28, 0.0);

    jac.set(Ca40, Si29, 0.0);

    jac.set(Ca40, Si30, 0.0);

    jac.set(Ca40, Si31, 0.0);

    jac.set(Ca40, P29, 0.0);

    jac.set(Ca40, P31, 0.0);

    jac.set(Ca40, P33, 0.0);

    jac.set(Ca40, S32, 0.0);

    jac.set(Ca40, S33, 0.0);

    jac.set(Ca40, S34, 0.0);

    jac.set(Ca40, S35, 0.0);

    jac.set(Ca40, Cl33, 0.0);

    jac.set(Ca40, Cl35, 0.0);

    jac.set(Ca40, Cl37, 0.0);

    scratch = screened_rates(k_Ar36_He4_to_Ca40_reaclib)*Y(He4)*state.rho;
    jac.set(Ca40, Ar36, scratch);

    scratch = screened_rates(k_Ar37_He4_to_n_Ca40_reaclib)*Y(He4)*state.rho;
    jac.set(Ca40, Ar37, scratch);

    jac.set(Ca40, Ar38, 0.0);

    scratch = screened_rates(k_K37_He4_to_p_Ca40_reaclib)*Y(He4)*state.rho;
    jac.set(Ca40, K37, scratch);

    scratch = screened_rates(k_K39_p_to_Ca40_reaclib)*Y(H1)*state.rho;
    jac.set(Ca40, K39, scratch);

    scratch = -screened_rates(k_Ca40_He4_to_Ti44_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ca40_He4_to_p_Sc43_reaclib)*Y(He4)*state.rho
              - 0.5*screened_rates(k_Ca40_n_n_to_Ca42_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_Ca40_n_to_He4_Ar37_reaclib)*Y(N)*state.rho
              - screened_rates(k_Ca40_p_to_He4_K37_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Ca40_p_to_Sc41_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Ca40_to_He4_Ar36_reaclib)
              - screened_rates(k_Ca40_to_p_K39_reaclib);
    jac.set(Ca40, Ca40, scratch);

    scratch = screened_rates(k_Ca42_to_Ca40_n_n_approx);
    jac.set(Ca40, Ca42, scratch);

    scratch = screened_rates(k_Sc41_to_p_Ca40_reaclib);
    jac.set(Ca40, Sc41, scratch);

    scratch = screened_rates(k_Sc43_p_to_He4_Ca40_reaclib)*Y(H1)*state.rho;
    jac.set(Ca40, Sc43, scratch);

    jac.set(Ca40, Sc45, 0.0);

    scratch = screened_rates(k_Ti44_to_He4_Ca40_reaclib);
    jac.set(Ca40, Ti44, scratch);

    jac.set(Ca40, Ti46, 0.0);

    jac.set(Ca40, Ti48, 0.0);

    jac.set(Ca40, V47, 0.0);

    jac.set(Ca40, V49, 0.0);

    jac.set(Ca40, V51, 0.0);

    jac.set(Ca40, Cr48, 0.0);

    jac.set(Ca40, Cr50, 0.0);

    jac.set(Ca40, Cr52, 0.0);

    jac.set(Ca40, Cr54, 0.0);

    jac.set(Ca40, Mn51, 0.0);

    jac.set(Ca40, Mn53, 0.0);

    jac.set(Ca40, Mn55, 0.0);

    jac.set(Ca40, Fe52, 0.0);

    jac.set(Ca40, Fe54, 0.0);

    jac.set(Ca40, Fe56, 0.0);

    jac.set(Ca40, Fe58, 0.0);

    jac.set(Ca40, Co55, 0.0);

    jac.set(Ca40, Co56, 0.0);

    jac.set(Ca40, Co57, 0.0);

    jac.set(Ca40, Co59, 0.0);

    jac.set(Ca40, Ni56, 0.0);

    jac.set(Ca40, Ni58, 0.0);

    jac.set(Ca40, Ni60, 0.0);

    jac.set(Ca40, Cu59, 0.0);

    jac.set(Ca40, Zn60, 0.0);


    // Ca42 row

    scratch = 1.0*screened_rates(k_Ca40_n_n_to_Ca42_approx)*Y(Ca40)*Y(N)*amrex::Math::powi<2>(state.rho)
              + 0.5*Y(Ca40)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Ca40_n_n_to_Ca42_approx_dYN
              - Y(Ca42)*rate_eval.drate_Ca42_to_Ca40_n_n_approx_dYN;
    jac.set(Ca42, N, scratch);

    scratch = -screened_rates(k_Ca42_p_to_He4_K39_reaclib)*Y(Ca42)*state.rho
              - screened_rates(k_Ca42_p_to_Sc43_reaclib)*Y(Ca42)*state.rho
              + screened_rates(k_Sc45_p_to_He4_Ca42_reaclib)*Y(Sc45)*state.rho;
    jac.set(Ca42, H1, scratch);

    scratch = screened_rates(k_Ar38_He4_to_Ca42_reaclib)*Y(Ar38)*state.rho
              - screened_rates(k_Ca42_He4_to_Ti46_reaclib)*Y(Ca42)*state.rho
              - screened_rates(k_Ca42_He4_to_p_Sc45_reaclib)*Y(Ca42)*state.rho
              + screened_rates(k_K39_He4_to_p_Ca42_reaclib)*Y(K39)*state.rho;
    jac.set(Ca42, He4, scratch);

    jac.set(Ca42, C12, 0.0);

    jac.set(Ca42, O16, 0.0);

    jac.set(Ca42, Ne20, 0.0);

    jac.set(Ca42, Mg23, 0.0);

    jac.set(Ca42, Mg24, 0.0);

    jac.set(Ca42, Al27, 0.0);

    jac.set(Ca42, Si27, 0.0);

    jac.set(Ca42, Si28, 0.0);

    jac.set(Ca42, Si29, 0.0);

    jac.set(Ca42, Si30, 0.0);

    jac.set(Ca42, Si31, 0.0);

    jac.set(Ca42, P29, 0.0);

    jac.set(Ca42, P31, 0.0);

    jac.set(Ca42, P33, 0.0);

    jac.set(Ca42, S32, 0.0);

    jac.set(Ca42, S33, 0.0);

    jac.set(Ca42, S34, 0.0);

    jac.set(Ca42, S35, 0.0);

    jac.set(Ca42, Cl33, 0.0);

    jac.set(Ca42, Cl35, 0.0);

    jac.set(Ca42, Cl37, 0.0);

    jac.set(Ca42, Ar36, 0.0);

    jac.set(Ca42, Ar37, 0.0);

    scratch = screened_rates(k_Ar38_He4_to_Ca42_reaclib)*Y(He4)*state.rho;
    jac.set(Ca42, Ar38, scratch);

    jac.set(Ca42, K37, 0.0);

    scratch = screened_rates(k_K39_He4_to_p_Ca42_reaclib)*Y(He4)*state.rho;
    jac.set(Ca42, K39, scratch);

    scratch = 0.5*screened_rates(k_Ca40_n_n_to_Ca42_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(Ca42, Ca40, scratch);

    scratch = -screened_rates(k_Ca42_He4_to_Ti46_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ca42_He4_to_p_Sc45_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ca42_p_to_He4_K39_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Ca42_p_to_Sc43_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Ca42_to_Ca40_n_n_approx)
              - screened_rates(k_Ca42_to_He4_Ar38_reaclib);
    jac.set(Ca42, Ca42, scratch);

    jac.set(Ca42, Sc41, 0.0);

    scratch = screened_rates(k_Sc43_to_p_Ca42_reaclib);
    jac.set(Ca42, Sc43, scratch);

    scratch = screened_rates(k_Sc45_p_to_He4_Ca42_reaclib)*Y(H1)*state.rho;
    jac.set(Ca42, Sc45, scratch);

    jac.set(Ca42, Ti44, 0.0);

    scratch = screened_rates(k_Ti46_to_He4_Ca42_reaclib);
    jac.set(Ca42, Ti46, scratch);

    jac.set(Ca42, Ti48, 0.0);

    jac.set(Ca42, V47, 0.0);

    jac.set(Ca42, V49, 0.0);

    jac.set(Ca42, V51, 0.0);

    jac.set(Ca42, Cr48, 0.0);

    jac.set(Ca42, Cr50, 0.0);

    jac.set(Ca42, Cr52, 0.0);

    jac.set(Ca42, Cr54, 0.0);

    jac.set(Ca42, Mn51, 0.0);

    jac.set(Ca42, Mn53, 0.0);

    jac.set(Ca42, Mn55, 0.0);

    jac.set(Ca42, Fe52, 0.0);

    jac.set(Ca42, Fe54, 0.0);

    jac.set(Ca42, Fe56, 0.0);

    jac.set(Ca42, Fe58, 0.0);

    jac.set(Ca42, Co55, 0.0);

    jac.set(Ca42, Co56, 0.0);

    jac.set(Ca42, Co57, 0.0);

    jac.set(Ca42, Co59, 0.0);

    jac.set(Ca42, Ni56, 0.0);

    jac.set(Ca42, Ni58, 0.0);

    jac.set(Ca42, Ni60, 0.0);

    jac.set(Ca42, Cu59, 0.0);

    jac.set(Ca42, Zn60, 0.0);


    // Sc41 row

    scratch = -1.0*screened_rates(k_Sc41_n_n_to_Sc43_approx)*Y(Sc41)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 0.5*Y(Sc41)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Sc41_n_n_to_Sc43_approx_dYN
              + Y(Sc43)*rate_eval.drate_Sc43_to_Sc41_n_n_approx_dYN;
    jac.set(Sc41, N, scratch);

    scratch = screened_rates(k_Ca40_p_to_Sc41_reaclib)*Y(Ca40)*state.rho
              + screened_rates(k_Ti44_p_to_He4_Sc41_reaclib)*Y(Ti44)*state.rho;
    jac.set(Sc41, H1, scratch);

    scratch = screened_rates(k_K37_He4_to_Sc41_reaclib)*Y(K37)*state.rho
              - screened_rates(k_Sc41_He4_to_p_Ti44_reaclib)*Y(Sc41)*state.rho;
    jac.set(Sc41, He4, scratch);

    jac.set(Sc41, C12, 0.0);

    jac.set(Sc41, O16, 0.0);

    jac.set(Sc41, Ne20, 0.0);

    jac.set(Sc41, Mg23, 0.0);

    jac.set(Sc41, Mg24, 0.0);

    jac.set(Sc41, Al27, 0.0);

    jac.set(Sc41, Si27, 0.0);

    jac.set(Sc41, Si28, 0.0);

    jac.set(Sc41, Si29, 0.0);

    jac.set(Sc41, Si30, 0.0);

    jac.set(Sc41, Si31, 0.0);

    jac.set(Sc41, P29, 0.0);

    jac.set(Sc41, P31, 0.0);

    jac.set(Sc41, P33, 0.0);

    jac.set(Sc41, S32, 0.0);

    jac.set(Sc41, S33, 0.0);

    jac.set(Sc41, S34, 0.0);

    jac.set(Sc41, S35, 0.0);

    jac.set(Sc41, Cl33, 0.0);

    jac.set(Sc41, Cl35, 0.0);

    jac.set(Sc41, Cl37, 0.0);

    jac.set(Sc41, Ar36, 0.0);

    jac.set(Sc41, Ar37, 0.0);

    jac.set(Sc41, Ar38, 0.0);

    scratch = screened_rates(k_K37_He4_to_Sc41_reaclib)*Y(He4)*state.rho;
    jac.set(Sc41, K37, scratch);

    jac.set(Sc41, K39, 0.0);

    scratch = screened_rates(k_Ca40_p_to_Sc41_reaclib)*Y(H1)*state.rho;
    jac.set(Sc41, Ca40, scratch);

    jac.set(Sc41, Ca42, 0.0);

    scratch = -screened_rates(k_Sc41_He4_to_p_Ti44_reaclib)*Y(He4)*state.rho
              - 0.5*screened_rates(k_Sc41_n_n_to_Sc43_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_Sc41_to_He4_K37_reaclib)
              - screened_rates(k_Sc41_to_p_Ca40_reaclib);
    jac.set(Sc41, Sc41, scratch);

    scratch = screened_rates(k_Sc43_to_Sc41_n_n_approx);
    jac.set(Sc41, Sc43, scratch);

    jac.set(Sc41, Sc45, 0.0);

    scratch = screened_rates(k_Ti44_p_to_He4_Sc41_reaclib)*Y(H1)*state.rho;
    jac.set(Sc41, Ti44, scratch);

    jac.set(Sc41, Ti46, 0.0);

    jac.set(Sc41, Ti48, 0.0);

    jac.set(Sc41, V47, 0.0);

    jac.set(Sc41, V49, 0.0);

    jac.set(Sc41, V51, 0.0);

    jac.set(Sc41, Cr48, 0.0);

    jac.set(Sc41, Cr50, 0.0);

    jac.set(Sc41, Cr52, 0.0);

    jac.set(Sc41, Cr54, 0.0);

    jac.set(Sc41, Mn51, 0.0);

    jac.set(Sc41, Mn53, 0.0);

    jac.set(Sc41, Mn55, 0.0);

    jac.set(Sc41, Fe52, 0.0);

    jac.set(Sc41, Fe54, 0.0);

    jac.set(Sc41, Fe56, 0.0);

    jac.set(Sc41, Fe58, 0.0);

    jac.set(Sc41, Co55, 0.0);

    jac.set(Sc41, Co56, 0.0);

    jac.set(Sc41, Co57, 0.0);

    jac.set(Sc41, Co59, 0.0);

    jac.set(Sc41, Ni56, 0.0);

    jac.set(Sc41, Ni58, 0.0);

    jac.set(Sc41, Ni60, 0.0);

    jac.set(Sc41, Cu59, 0.0);

    jac.set(Sc41, Zn60, 0.0);


    // Sc43 row

    scratch = 1.0*screened_rates(k_Sc41_n_n_to_Sc43_approx)*Y(Sc41)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 1.0*screened_rates(k_Sc43_n_n_to_Sc45_approx)*Y(Sc43)*Y(N)*amrex::Math::powi<2>(state.rho)
              + 0.5*Y(Sc41)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Sc41_n_n_to_Sc43_approx_dYN
              - 0.5*Y(Sc43)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Sc43_n_n_to_Sc45_approx_dYN
              - Y(Sc43)*rate_eval.drate_Sc43_to_Sc41_n_n_approx_dYN
              + Y(Sc45)*rate_eval.drate_Sc45_to_Sc43_n_n_approx_dYN;
    jac.set(Sc43, N, scratch);

    scratch = screened_rates(k_Ca42_p_to_Sc43_reaclib)*Y(Ca42)*state.rho
              - screened_rates(k_Sc43_p_to_He4_Ca40_reaclib)*Y(Sc43)*state.rho
              - screened_rates(k_Sc43_p_to_Ti44_reaclib)*Y(Sc43)*state.rho
              + screened_rates(k_Ti46_p_to_He4_Sc43_reaclib)*Y(Ti46)*state.rho;
    jac.set(Sc43, H1, scratch);

    scratch = screened_rates(k_Ca40_He4_to_p_Sc43_reaclib)*Y(Ca40)*state.rho
              + screened_rates(k_K39_He4_to_Sc43_reaclib)*Y(K39)*state.rho
              - screened_rates(k_Sc43_He4_to_V47_reaclib)*Y(Sc43)*state.rho
              - screened_rates(k_Sc43_He4_to_p_Ti46_reaclib)*Y(Sc43)*state.rho;
    jac.set(Sc43, He4, scratch);

    jac.set(Sc43, C12, 0.0);

    jac.set(Sc43, O16, 0.0);

    jac.set(Sc43, Ne20, 0.0);

    jac.set(Sc43, Mg23, 0.0);

    jac.set(Sc43, Mg24, 0.0);

    jac.set(Sc43, Al27, 0.0);

    jac.set(Sc43, Si27, 0.0);

    jac.set(Sc43, Si28, 0.0);

    jac.set(Sc43, Si29, 0.0);

    jac.set(Sc43, Si30, 0.0);

    jac.set(Sc43, Si31, 0.0);

    jac.set(Sc43, P29, 0.0);

    jac.set(Sc43, P31, 0.0);

    jac.set(Sc43, P33, 0.0);

    jac.set(Sc43, S32, 0.0);

    jac.set(Sc43, S33, 0.0);

    jac.set(Sc43, S34, 0.0);

    jac.set(Sc43, S35, 0.0);

    jac.set(Sc43, Cl33, 0.0);

    jac.set(Sc43, Cl35, 0.0);

    jac.set(Sc43, Cl37, 0.0);

    jac.set(Sc43, Ar36, 0.0);

    jac.set(Sc43, Ar37, 0.0);

    jac.set(Sc43, Ar38, 0.0);

    jac.set(Sc43, K37, 0.0);

    scratch = screened_rates(k_K39_He4_to_Sc43_reaclib)*Y(He4)*state.rho;
    jac.set(Sc43, K39, scratch);

    scratch = screened_rates(k_Ca40_He4_to_p_Sc43_reaclib)*Y(He4)*state.rho;
    jac.set(Sc43, Ca40, scratch);

    scratch = screened_rates(k_Ca42_p_to_Sc43_reaclib)*Y(H1)*state.rho;
    jac.set(Sc43, Ca42, scratch);

    scratch = 0.5*screened_rates(k_Sc41_n_n_to_Sc43_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(Sc43, Sc41, scratch);

    scratch = -screened_rates(k_Sc43_He4_to_V47_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Sc43_He4_to_p_Ti46_reaclib)*Y(He4)*state.rho
              - 0.5*screened_rates(k_Sc43_n_n_to_Sc45_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_Sc43_p_to_He4_Ca40_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Sc43_p_to_Ti44_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Sc43_to_He4_K39_reaclib)
              - screened_rates(k_Sc43_to_Sc41_n_n_approx)
              - screened_rates(k_Sc43_to_p_Ca42_reaclib);
    jac.set(Sc43, Sc43, scratch);

    scratch = screened_rates(k_Sc45_to_Sc43_n_n_approx);
    jac.set(Sc43, Sc45, scratch);

    scratch = screened_rates(k_Ti44_to_p_Sc43_reaclib);
    jac.set(Sc43, Ti44, scratch);

    scratch = screened_rates(k_Ti46_p_to_He4_Sc43_reaclib)*Y(H1)*state.rho;
    jac.set(Sc43, Ti46, scratch);

    jac.set(Sc43, Ti48, 0.0);

    scratch = screened_rates(k_V47_to_He4_Sc43_reaclib);
    jac.set(Sc43, V47, scratch);

    jac.set(Sc43, V49, 0.0);

    jac.set(Sc43, V51, 0.0);

    jac.set(Sc43, Cr48, 0.0);

    jac.set(Sc43, Cr50, 0.0);

    jac.set(Sc43, Cr52, 0.0);

    jac.set(Sc43, Cr54, 0.0);

    jac.set(Sc43, Mn51, 0.0);

    jac.set(Sc43, Mn53, 0.0);

    jac.set(Sc43, Mn55, 0.0);

    jac.set(Sc43, Fe52, 0.0);

    jac.set(Sc43, Fe54, 0.0);

    jac.set(Sc43, Fe56, 0.0);

    jac.set(Sc43, Fe58, 0.0);

    jac.set(Sc43, Co55, 0.0);

    jac.set(Sc43, Co56, 0.0);

    jac.set(Sc43, Co57, 0.0);

    jac.set(Sc43, Co59, 0.0);

    jac.set(Sc43, Ni56, 0.0);

    jac.set(Sc43, Ni58, 0.0);

    jac.set(Sc43, Ni60, 0.0);

    jac.set(Sc43, Cu59, 0.0);

    jac.set(Sc43, Zn60, 0.0);


    // Sc45 row

    scratch = 1.0*screened_rates(k_Sc43_n_n_to_Sc45_approx)*Y(Sc43)*Y(N)*amrex::Math::powi<2>(state.rho)
              + 0.5*Y(Sc43)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Sc43_n_n_to_Sc45_approx_dYN
              - Y(Sc45)*rate_eval.drate_Sc45_to_Sc43_n_n_approx_dYN;
    jac.set(Sc45, N, scratch);

    scratch = -screened_rates(k_Sc45_p_to_He4_Ca42_reaclib)*Y(Sc45)*state.rho
              - screened_rates(k_Sc45_p_to_Ti46_reaclib)*Y(Sc45)*state.rho
              + screened_rates(k_Ti48_p_to_He4_Sc45_reaclib)*Y(Ti48)*state.rho;
    jac.set(Sc45, H1, scratch);

    scratch = screened_rates(k_Ca42_He4_to_p_Sc45_reaclib)*Y(Ca42)*state.rho
              - screened_rates(k_Sc45_He4_to_V49_reaclib)*Y(Sc45)*state.rho
              - screened_rates(k_Sc45_He4_to_p_Ti48_reaclib)*Y(Sc45)*state.rho;
    jac.set(Sc45, He4, scratch);

    jac.set(Sc45, C12, 0.0);

    jac.set(Sc45, O16, 0.0);

    jac.set(Sc45, Ne20, 0.0);

    jac.set(Sc45, Mg23, 0.0);

    jac.set(Sc45, Mg24, 0.0);

    jac.set(Sc45, Al27, 0.0);

    jac.set(Sc45, Si27, 0.0);

    jac.set(Sc45, Si28, 0.0);

    jac.set(Sc45, Si29, 0.0);

    jac.set(Sc45, Si30, 0.0);

    jac.set(Sc45, Si31, 0.0);

    jac.set(Sc45, P29, 0.0);

    jac.set(Sc45, P31, 0.0);

    jac.set(Sc45, P33, 0.0);

    jac.set(Sc45, S32, 0.0);

    jac.set(Sc45, S33, 0.0);

    jac.set(Sc45, S34, 0.0);

    jac.set(Sc45, S35, 0.0);

    jac.set(Sc45, Cl33, 0.0);

    jac.set(Sc45, Cl35, 0.0);

    jac.set(Sc45, Cl37, 0.0);

    jac.set(Sc45, Ar36, 0.0);

    jac.set(Sc45, Ar37, 0.0);

    jac.set(Sc45, Ar38, 0.0);

    jac.set(Sc45, K37, 0.0);

    jac.set(Sc45, K39, 0.0);

    jac.set(Sc45, Ca40, 0.0);

    scratch = screened_rates(k_Ca42_He4_to_p_Sc45_reaclib)*Y(He4)*state.rho;
    jac.set(Sc45, Ca42, scratch);

    jac.set(Sc45, Sc41, 0.0);

    scratch = 0.5*screened_rates(k_Sc43_n_n_to_Sc45_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(Sc45, Sc43, scratch);

    scratch = -screened_rates(k_Sc45_He4_to_V49_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Sc45_He4_to_p_Ti48_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Sc45_p_to_He4_Ca42_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Sc45_p_to_Ti46_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Sc45_to_Sc43_n_n_approx);
    jac.set(Sc45, Sc45, scratch);

    jac.set(Sc45, Ti44, 0.0);

    scratch = screened_rates(k_Ti46_to_p_Sc45_reaclib);
    jac.set(Sc45, Ti46, scratch);

    scratch = screened_rates(k_Ti48_p_to_He4_Sc45_reaclib)*Y(H1)*state.rho;
    jac.set(Sc45, Ti48, scratch);

    jac.set(Sc45, V47, 0.0);

    scratch = screened_rates(k_V49_to_He4_Sc45_reaclib);
    jac.set(Sc45, V49, scratch);

    jac.set(Sc45, V51, 0.0);

    jac.set(Sc45, Cr48, 0.0);

    jac.set(Sc45, Cr50, 0.0);

    jac.set(Sc45, Cr52, 0.0);

    jac.set(Sc45, Cr54, 0.0);

    jac.set(Sc45, Mn51, 0.0);

    jac.set(Sc45, Mn53, 0.0);

    jac.set(Sc45, Mn55, 0.0);

    jac.set(Sc45, Fe52, 0.0);

    jac.set(Sc45, Fe54, 0.0);

    jac.set(Sc45, Fe56, 0.0);

    jac.set(Sc45, Fe58, 0.0);

    jac.set(Sc45, Co55, 0.0);

    jac.set(Sc45, Co56, 0.0);

    jac.set(Sc45, Co57, 0.0);

    jac.set(Sc45, Co59, 0.0);

    jac.set(Sc45, Ni56, 0.0);

    jac.set(Sc45, Ni58, 0.0);

    jac.set(Sc45, Ni60, 0.0);

    jac.set(Sc45, Cu59, 0.0);

    jac.set(Sc45, Zn60, 0.0);


    // Ti44 row

    scratch = -1.0*screened_rates(k_Ti44_n_n_to_Ti46_approx)*Y(Ti44)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 0.5*Y(Ti44)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Ti44_n_n_to_Ti46_approx_dYN
              + Y(Ti46)*rate_eval.drate_Ti46_to_Ti44_n_n_approx_dYN;
    jac.set(Ti44, N, scratch);

    scratch = screened_rates(k_Sc43_p_to_Ti44_reaclib)*Y(Sc43)*state.rho
              - screened_rates(k_Ti44_p_to_He4_Sc41_reaclib)*Y(Ti44)*state.rho
              + screened_rates(k_V47_p_to_He4_Ti44_reaclib)*Y(V47)*state.rho;
    jac.set(Ti44, H1, scratch);

    scratch = screened_rates(k_Ca40_He4_to_Ti44_reaclib)*Y(Ca40)*state.rho
              + screened_rates(k_Sc41_He4_to_p_Ti44_reaclib)*Y(Sc41)*state.rho
              - screened_rates(k_Ti44_He4_to_Cr48_reaclib)*Y(Ti44)*state.rho
              - screened_rates(k_Ti44_He4_to_p_V47_reaclib)*Y(Ti44)*state.rho;
    jac.set(Ti44, He4, scratch);

    jac.set(Ti44, C12, 0.0);

    jac.set(Ti44, O16, 0.0);

    jac.set(Ti44, Ne20, 0.0);

    jac.set(Ti44, Mg23, 0.0);

    jac.set(Ti44, Mg24, 0.0);

    jac.set(Ti44, Al27, 0.0);

    jac.set(Ti44, Si27, 0.0);

    jac.set(Ti44, Si28, 0.0);

    jac.set(Ti44, Si29, 0.0);

    jac.set(Ti44, Si30, 0.0);

    jac.set(Ti44, Si31, 0.0);

    jac.set(Ti44, P29, 0.0);

    jac.set(Ti44, P31, 0.0);

    jac.set(Ti44, P33, 0.0);

    jac.set(Ti44, S32, 0.0);

    jac.set(Ti44, S33, 0.0);

    jac.set(Ti44, S34, 0.0);

    jac.set(Ti44, S35, 0.0);

    jac.set(Ti44, Cl33, 0.0);

    jac.set(Ti44, Cl35, 0.0);

    jac.set(Ti44, Cl37, 0.0);

    jac.set(Ti44, Ar36, 0.0);

    jac.set(Ti44, Ar37, 0.0);

    jac.set(Ti44, Ar38, 0.0);

    jac.set(Ti44, K37, 0.0);

    jac.set(Ti44, K39, 0.0);

    scratch = screened_rates(k_Ca40_He4_to_Ti44_reaclib)*Y(He4)*state.rho;
    jac.set(Ti44, Ca40, scratch);

    jac.set(Ti44, Ca42, 0.0);

    scratch = screened_rates(k_Sc41_He4_to_p_Ti44_reaclib)*Y(He4)*state.rho;
    jac.set(Ti44, Sc41, scratch);

    scratch = screened_rates(k_Sc43_p_to_Ti44_reaclib)*Y(H1)*state.rho;
    jac.set(Ti44, Sc43, scratch);

    jac.set(Ti44, Sc45, 0.0);

    scratch = -screened_rates(k_Ti44_He4_to_Cr48_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ti44_He4_to_p_V47_reaclib)*Y(He4)*state.rho
              - 0.5*screened_rates(k_Ti44_n_n_to_Ti46_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_Ti44_p_to_He4_Sc41_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Ti44_to_He4_Ca40_reaclib)
              - screened_rates(k_Ti44_to_p_Sc43_reaclib);
    jac.set(Ti44, Ti44, scratch);

    scratch = screened_rates(k_Ti46_to_Ti44_n_n_approx);
    jac.set(Ti44, Ti46, scratch);

    jac.set(Ti44, Ti48, 0.0);

    scratch = screened_rates(k_V47_p_to_He4_Ti44_reaclib)*Y(H1)*state.rho;
    jac.set(Ti44, V47, scratch);

    jac.set(Ti44, V49, 0.0);

    jac.set(Ti44, V51, 0.0);

    scratch = screened_rates(k_Cr48_to_He4_Ti44_reaclib);
    jac.set(Ti44, Cr48, scratch);

    jac.set(Ti44, Cr50, 0.0);

    jac.set(Ti44, Cr52, 0.0);

    jac.set(Ti44, Cr54, 0.0);

    jac.set(Ti44, Mn51, 0.0);

    jac.set(Ti44, Mn53, 0.0);

    jac.set(Ti44, Mn55, 0.0);

    jac.set(Ti44, Fe52, 0.0);

    jac.set(Ti44, Fe54, 0.0);

    jac.set(Ti44, Fe56, 0.0);

    jac.set(Ti44, Fe58, 0.0);

    jac.set(Ti44, Co55, 0.0);

    jac.set(Ti44, Co56, 0.0);

    jac.set(Ti44, Co57, 0.0);

    jac.set(Ti44, Co59, 0.0);

    jac.set(Ti44, Ni56, 0.0);

    jac.set(Ti44, Ni58, 0.0);

    jac.set(Ti44, Ni60, 0.0);

    jac.set(Ti44, Cu59, 0.0);

    jac.set(Ti44, Zn60, 0.0);


    // Ti46 row

    scratch = 1.0*screened_rates(k_Ti44_n_n_to_Ti46_approx)*Y(Ti44)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 1.0*screened_rates(k_Ti46_n_n_to_Ti48_approx)*Y(Ti46)*Y(N)*amrex::Math::powi<2>(state.rho)
              + 0.5*Y(Ti44)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Ti44_n_n_to_Ti46_approx_dYN
              - 0.5*Y(Ti46)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Ti46_n_n_to_Ti48_approx_dYN
              - Y(Ti46)*rate_eval.drate_Ti46_to_Ti44_n_n_approx_dYN
              + Y(Ti48)*rate_eval.drate_Ti48_to_Ti46_n_n_approx_dYN;
    jac.set(Ti46, N, scratch);

    scratch = screened_rates(k_Sc45_p_to_Ti46_reaclib)*Y(Sc45)*state.rho
              - screened_rates(k_Ti46_p_to_He4_Sc43_reaclib)*Y(Ti46)*state.rho
              - screened_rates(k_Ti46_p_to_V47_reaclib)*Y(Ti46)*state.rho
              + screened_rates(k_V49_p_to_He4_Ti46_reaclib)*Y(V49)*state.rho;
    jac.set(Ti46, H1, scratch);

    scratch = screened_rates(k_Ca42_He4_to_Ti46_reaclib)*Y(Ca42)*state.rho
              + screened_rates(k_Sc43_He4_to_p_Ti46_reaclib)*Y(Sc43)*state.rho
              - screened_rates(k_Ti46_He4_to_Cr50_reaclib)*Y(Ti46)*state.rho
              - screened_rates(k_Ti46_He4_to_p_V49_reaclib)*Y(Ti46)*state.rho;
    jac.set(Ti46, He4, scratch);

    jac.set(Ti46, C12, 0.0);

    jac.set(Ti46, O16, 0.0);

    jac.set(Ti46, Ne20, 0.0);

    jac.set(Ti46, Mg23, 0.0);

    jac.set(Ti46, Mg24, 0.0);

    jac.set(Ti46, Al27, 0.0);

    jac.set(Ti46, Si27, 0.0);

    jac.set(Ti46, Si28, 0.0);

    jac.set(Ti46, Si29, 0.0);

    jac.set(Ti46, Si30, 0.0);

    jac.set(Ti46, Si31, 0.0);

    jac.set(Ti46, P29, 0.0);

    jac.set(Ti46, P31, 0.0);

    jac.set(Ti46, P33, 0.0);

    jac.set(Ti46, S32, 0.0);

    jac.set(Ti46, S33, 0.0);

    jac.set(Ti46, S34, 0.0);

    jac.set(Ti46, S35, 0.0);

    jac.set(Ti46, Cl33, 0.0);

    jac.set(Ti46, Cl35, 0.0);

    jac.set(Ti46, Cl37, 0.0);

    jac.set(Ti46, Ar36, 0.0);

    jac.set(Ti46, Ar37, 0.0);

    jac.set(Ti46, Ar38, 0.0);

    jac.set(Ti46, K37, 0.0);

    jac.set(Ti46, K39, 0.0);

    jac.set(Ti46, Ca40, 0.0);

    scratch = screened_rates(k_Ca42_He4_to_Ti46_reaclib)*Y(He4)*state.rho;
    jac.set(Ti46, Ca42, scratch);

    jac.set(Ti46, Sc41, 0.0);

    scratch = screened_rates(k_Sc43_He4_to_p_Ti46_reaclib)*Y(He4)*state.rho;
    jac.set(Ti46, Sc43, scratch);

    scratch = screened_rates(k_Sc45_p_to_Ti46_reaclib)*Y(H1)*state.rho;
    jac.set(Ti46, Sc45, scratch);

    scratch = 0.5*screened_rates(k_Ti44_n_n_to_Ti46_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(Ti46, Ti44, scratch);

    scratch = -screened_rates(k_Ti46_He4_to_Cr50_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ti46_He4_to_p_V49_reaclib)*Y(He4)*state.rho
              - 0.5*screened_rates(k_Ti46_n_n_to_Ti48_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_Ti46_p_to_He4_Sc43_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Ti46_p_to_V47_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Ti46_to_He4_Ca42_reaclib)
              - screened_rates(k_Ti46_to_Ti44_n_n_approx)
              - screened_rates(k_Ti46_to_p_Sc45_reaclib);
    jac.set(Ti46, Ti46, scratch);

    scratch = screened_rates(k_Ti48_to_Ti46_n_n_approx);
    jac.set(Ti46, Ti48, scratch);

    scratch = screened_rates(k_V47_to_p_Ti46_reaclib);
    jac.set(Ti46, V47, scratch);

    scratch = screened_rates(k_V49_p_to_He4_Ti46_reaclib)*Y(H1)*state.rho;
    jac.set(Ti46, V49, scratch);

    jac.set(Ti46, V51, 0.0);

    jac.set(Ti46, Cr48, 0.0);

    scratch = screened_rates(k_Cr50_to_He4_Ti46_reaclib);
    jac.set(Ti46, Cr50, scratch);

    jac.set(Ti46, Cr52, 0.0);

    jac.set(Ti46, Cr54, 0.0);

    jac.set(Ti46, Mn51, 0.0);

    jac.set(Ti46, Mn53, 0.0);

    jac.set(Ti46, Mn55, 0.0);

    jac.set(Ti46, Fe52, 0.0);

    jac.set(Ti46, Fe54, 0.0);

    jac.set(Ti46, Fe56, 0.0);

    jac.set(Ti46, Fe58, 0.0);

    jac.set(Ti46, Co55, 0.0);

    jac.set(Ti46, Co56, 0.0);

    jac.set(Ti46, Co57, 0.0);

    jac.set(Ti46, Co59, 0.0);

    jac.set(Ti46, Ni56, 0.0);

    jac.set(Ti46, Ni58, 0.0);

    jac.set(Ti46, Ni60, 0.0);

    jac.set(Ti46, Cu59, 0.0);

    jac.set(Ti46, Zn60, 0.0);


    // Ti48 row

    scratch = 1.0*screened_rates(k_Ti46_n_n_to_Ti48_approx)*Y(Ti46)*Y(N)*amrex::Math::powi<2>(state.rho)
              + 0.5*Y(Ti46)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Ti46_n_n_to_Ti48_approx_dYN
              - Y(Ti48)*rate_eval.drate_Ti48_to_Ti46_n_n_approx_dYN;
    jac.set(Ti48, N, scratch);

    scratch = -screened_rates(k_Ti48_p_to_He4_Sc45_reaclib)*Y(Ti48)*state.rho
              - screened_rates(k_Ti48_p_to_V49_reaclib)*Y(Ti48)*state.rho
              + screened_rates(k_V51_p_to_He4_Ti48_reaclib)*Y(V51)*state.rho;
    jac.set(Ti48, H1, scratch);

    scratch = screened_rates(k_Sc45_He4_to_p_Ti48_reaclib)*Y(Sc45)*state.rho
              - screened_rates(k_Ti48_He4_to_Cr52_reaclib)*Y(Ti48)*state.rho
              - screened_rates(k_Ti48_He4_to_p_V51_reaclib)*Y(Ti48)*state.rho;
    jac.set(Ti48, He4, scratch);

    jac.set(Ti48, C12, 0.0);

    jac.set(Ti48, O16, 0.0);

    jac.set(Ti48, Ne20, 0.0);

    jac.set(Ti48, Mg23, 0.0);

    jac.set(Ti48, Mg24, 0.0);

    jac.set(Ti48, Al27, 0.0);

    jac.set(Ti48, Si27, 0.0);

    jac.set(Ti48, Si28, 0.0);

    jac.set(Ti48, Si29, 0.0);

    jac.set(Ti48, Si30, 0.0);

    jac.set(Ti48, Si31, 0.0);

    jac.set(Ti48, P29, 0.0);

    jac.set(Ti48, P31, 0.0);

    jac.set(Ti48, P33, 0.0);

    jac.set(Ti48, S32, 0.0);

    jac.set(Ti48, S33, 0.0);

    jac.set(Ti48, S34, 0.0);

    jac.set(Ti48, S35, 0.0);

    jac.set(Ti48, Cl33, 0.0);

    jac.set(Ti48, Cl35, 0.0);

    jac.set(Ti48, Cl37, 0.0);

    jac.set(Ti48, Ar36, 0.0);

    jac.set(Ti48, Ar37, 0.0);

    jac.set(Ti48, Ar38, 0.0);

    jac.set(Ti48, K37, 0.0);

    jac.set(Ti48, K39, 0.0);

    jac.set(Ti48, Ca40, 0.0);

    jac.set(Ti48, Ca42, 0.0);

    jac.set(Ti48, Sc41, 0.0);

    jac.set(Ti48, Sc43, 0.0);

    scratch = screened_rates(k_Sc45_He4_to_p_Ti48_reaclib)*Y(He4)*state.rho;
    jac.set(Ti48, Sc45, scratch);

    jac.set(Ti48, Ti44, 0.0);

    scratch = 0.5*screened_rates(k_Ti46_n_n_to_Ti48_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(Ti48, Ti46, scratch);

    scratch = -screened_rates(k_Ti48_He4_to_Cr52_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ti48_He4_to_p_V51_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ti48_p_to_He4_Sc45_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Ti48_p_to_V49_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Ti48_to_Ti46_n_n_approx);
    jac.set(Ti48, Ti48, scratch);

    jac.set(Ti48, V47, 0.0);

    scratch = screened_rates(k_V49_to_p_Ti48_reaclib);
    jac.set(Ti48, V49, scratch);

    scratch = screened_rates(k_V51_p_to_He4_Ti48_reaclib)*Y(H1)*state.rho;
    jac.set(Ti48, V51, scratch);

    jac.set(Ti48, Cr48, 0.0);

    jac.set(Ti48, Cr50, 0.0);

    scratch = screened_rates(k_Cr52_to_He4_Ti48_reaclib);
    jac.set(Ti48, Cr52, scratch);

    jac.set(Ti48, Cr54, 0.0);

    jac.set(Ti48, Mn51, 0.0);

    jac.set(Ti48, Mn53, 0.0);

    jac.set(Ti48, Mn55, 0.0);

    jac.set(Ti48, Fe52, 0.0);

    jac.set(Ti48, Fe54, 0.0);

    jac.set(Ti48, Fe56, 0.0);

    jac.set(Ti48, Fe58, 0.0);

    jac.set(Ti48, Co55, 0.0);

    jac.set(Ti48, Co56, 0.0);

    jac.set(Ti48, Co57, 0.0);

    jac.set(Ti48, Co59, 0.0);

    jac.set(Ti48, Ni56, 0.0);

    jac.set(Ti48, Ni58, 0.0);

    jac.set(Ti48, Ni60, 0.0);

    jac.set(Ti48, Cu59, 0.0);

    jac.set(Ti48, Zn60, 0.0);


    // V47 row

    scratch = -1.0*screened_rates(k_V47_n_n_to_V49_approx)*Y(V47)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 0.5*Y(V47)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_V47_n_n_to_V49_approx_dYN
              + Y(V49)*rate_eval.drate_V49_to_V47_n_n_approx_dYN;
    jac.set(V47, N, scratch);

    scratch = screened_rates(k_Cr50_p_to_He4_V47_reaclib)*Y(Cr50)*state.rho
              + screened_rates(k_Ti46_p_to_V47_reaclib)*Y(Ti46)*state.rho
              - screened_rates(k_V47_p_to_Cr48_reaclib)*Y(V47)*state.rho
              - screened_rates(k_V47_p_to_He4_Ti44_reaclib)*Y(V47)*state.rho;
    jac.set(V47, H1, scratch);

    scratch = screened_rates(k_Sc43_He4_to_V47_reaclib)*Y(Sc43)*state.rho
              + screened_rates(k_Ti44_He4_to_p_V47_reaclib)*Y(Ti44)*state.rho
              - screened_rates(k_V47_He4_to_Mn51_reaclib)*Y(V47)*state.rho
              - screened_rates(k_V47_He4_to_p_Cr50_reaclib)*Y(V47)*state.rho;
    jac.set(V47, He4, scratch);

    jac.set(V47, C12, 0.0);

    jac.set(V47, O16, 0.0);

    jac.set(V47, Ne20, 0.0);

    jac.set(V47, Mg23, 0.0);

    jac.set(V47, Mg24, 0.0);

    jac.set(V47, Al27, 0.0);

    jac.set(V47, Si27, 0.0);

    jac.set(V47, Si28, 0.0);

    jac.set(V47, Si29, 0.0);

    jac.set(V47, Si30, 0.0);

    jac.set(V47, Si31, 0.0);

    jac.set(V47, P29, 0.0);

    jac.set(V47, P31, 0.0);

    jac.set(V47, P33, 0.0);

    jac.set(V47, S32, 0.0);

    jac.set(V47, S33, 0.0);

    jac.set(V47, S34, 0.0);

    jac.set(V47, S35, 0.0);

    jac.set(V47, Cl33, 0.0);

    jac.set(V47, Cl35, 0.0);

    jac.set(V47, Cl37, 0.0);

    jac.set(V47, Ar36, 0.0);

    jac.set(V47, Ar37, 0.0);

    jac.set(V47, Ar38, 0.0);

    jac.set(V47, K37, 0.0);

    jac.set(V47, K39, 0.0);

    jac.set(V47, Ca40, 0.0);

    jac.set(V47, Ca42, 0.0);

    jac.set(V47, Sc41, 0.0);

    scratch = screened_rates(k_Sc43_He4_to_V47_reaclib)*Y(He4)*state.rho;
    jac.set(V47, Sc43, scratch);

    jac.set(V47, Sc45, 0.0);

    scratch = screened_rates(k_Ti44_He4_to_p_V47_reaclib)*Y(He4)*state.rho;
    jac.set(V47, Ti44, scratch);

    scratch = screened_rates(k_Ti46_p_to_V47_reaclib)*Y(H1)*state.rho;
    jac.set(V47, Ti46, scratch);

    jac.set(V47, Ti48, 0.0);

    scratch = -screened_rates(k_V47_He4_to_Mn51_reaclib)*Y(He4)*state.rho
              - screened_rates(k_V47_He4_to_p_Cr50_reaclib)*Y(He4)*state.rho
              - 0.5*screened_rates(k_V47_n_n_to_V49_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_V47_p_to_Cr48_reaclib)*Y(H1)*state.rho
              - screened_rates(k_V47_p_to_He4_Ti44_reaclib)*Y(H1)*state.rho
              - screened_rates(k_V47_to_He4_Sc43_reaclib)
              - screened_rates(k_V47_to_p_Ti46_reaclib);
    jac.set(V47, V47, scratch);

    scratch = screened_rates(k_V49_to_V47_n_n_approx);
    jac.set(V47, V49, scratch);

    jac.set(V47, V51, 0.0);

    scratch = screened_rates(k_Cr48_to_p_V47_reaclib);
    jac.set(V47, Cr48, scratch);

    scratch = screened_rates(k_Cr50_p_to_He4_V47_reaclib)*Y(H1)*state.rho;
    jac.set(V47, Cr50, scratch);

    jac.set(V47, Cr52, 0.0);

    jac.set(V47, Cr54, 0.0);

    scratch = screened_rates(k_Mn51_to_He4_V47_reaclib);
    jac.set(V47, Mn51, scratch);

    jac.set(V47, Mn53, 0.0);

    jac.set(V47, Mn55, 0.0);

    jac.set(V47, Fe52, 0.0);

    jac.set(V47, Fe54, 0.0);

    jac.set(V47, Fe56, 0.0);

    jac.set(V47, Fe58, 0.0);

    jac.set(V47, Co55, 0.0);

    jac.set(V47, Co56, 0.0);

    jac.set(V47, Co57, 0.0);

    jac.set(V47, Co59, 0.0);

    jac.set(V47, Ni56, 0.0);

    jac.set(V47, Ni58, 0.0);

    jac.set(V47, Ni60, 0.0);

    jac.set(V47, Cu59, 0.0);

    jac.set(V47, Zn60, 0.0);


    // V49 row

    scratch = 1.0*screened_rates(k_V47_n_n_to_V49_approx)*Y(V47)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 1.0*screened_rates(k_V49_n_n_to_V51_approx)*Y(V49)*Y(N)*amrex::Math::powi<2>(state.rho)
              + 0.5*Y(V47)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_V47_n_n_to_V49_approx_dYN
              - 0.5*Y(V49)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_V49_n_n_to_V51_approx_dYN
              - Y(V49)*rate_eval.drate_V49_to_V47_n_n_approx_dYN
              + Y(V51)*rate_eval.drate_V51_to_V49_n_n_approx_dYN;
    jac.set(V49, N, scratch);

    scratch = screened_rates(k_Cr52_p_to_He4_V49_reaclib)*Y(Cr52)*state.rho
              + screened_rates(k_Ti48_p_to_V49_reaclib)*Y(Ti48)*state.rho
              - screened_rates(k_V49_p_to_Cr50_reaclib)*Y(V49)*state.rho
              - screened_rates(k_V49_p_to_He4_Ti46_reaclib)*Y(V49)*state.rho;
    jac.set(V49, H1, scratch);

    scratch = screened_rates(k_Sc45_He4_to_V49_reaclib)*Y(Sc45)*state.rho
              + screened_rates(k_Ti46_He4_to_p_V49_reaclib)*Y(Ti46)*state.rho
              - screened_rates(k_V49_He4_to_Mn53_reaclib)*Y(V49)*state.rho
              - screened_rates(k_V49_He4_to_p_Cr52_reaclib)*Y(V49)*state.rho;
    jac.set(V49, He4, scratch);

    jac.set(V49, C12, 0.0);

    jac.set(V49, O16, 0.0);

    jac.set(V49, Ne20, 0.0);

    jac.set(V49, Mg23, 0.0);

    jac.set(V49, Mg24, 0.0);

    jac.set(V49, Al27, 0.0);

    jac.set(V49, Si27, 0.0);

    jac.set(V49, Si28, 0.0);

    jac.set(V49, Si29, 0.0);

    jac.set(V49, Si30, 0.0);

    jac.set(V49, Si31, 0.0);

    jac.set(V49, P29, 0.0);

    jac.set(V49, P31, 0.0);

    jac.set(V49, P33, 0.0);

    jac.set(V49, S32, 0.0);

    jac.set(V49, S33, 0.0);

    jac.set(V49, S34, 0.0);

    jac.set(V49, S35, 0.0);

    jac.set(V49, Cl33, 0.0);

    jac.set(V49, Cl35, 0.0);

    jac.set(V49, Cl37, 0.0);

    jac.set(V49, Ar36, 0.0);

    jac.set(V49, Ar37, 0.0);

    jac.set(V49, Ar38, 0.0);

    jac.set(V49, K37, 0.0);

    jac.set(V49, K39, 0.0);

    jac.set(V49, Ca40, 0.0);

    jac.set(V49, Ca42, 0.0);

    jac.set(V49, Sc41, 0.0);

    jac.set(V49, Sc43, 0.0);

    scratch = screened_rates(k_Sc45_He4_to_V49_reaclib)*Y(He4)*state.rho;
    jac.set(V49, Sc45, scratch);

    jac.set(V49, Ti44, 0.0);

    scratch = screened_rates(k_Ti46_He4_to_p_V49_reaclib)*Y(He4)*state.rho;
    jac.set(V49, Ti46, scratch);

    scratch = screened_rates(k_Ti48_p_to_V49_reaclib)*Y(H1)*state.rho;
    jac.set(V49, Ti48, scratch);

    scratch = 0.5*screened_rates(k_V47_n_n_to_V49_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(V49, V47, scratch);

    scratch = -screened_rates(k_V49_He4_to_Mn53_reaclib)*Y(He4)*state.rho
              - screened_rates(k_V49_He4_to_p_Cr52_reaclib)*Y(He4)*state.rho
              - 0.5*screened_rates(k_V49_n_n_to_V51_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_V49_p_to_Cr50_reaclib)*Y(H1)*state.rho
              - screened_rates(k_V49_p_to_He4_Ti46_reaclib)*Y(H1)*state.rho
              - screened_rates(k_V49_to_He4_Sc45_reaclib)
              - screened_rates(k_V49_to_V47_n_n_approx)
              - screened_rates(k_V49_to_p_Ti48_reaclib);
    jac.set(V49, V49, scratch);

    scratch = screened_rates(k_V51_to_V49_n_n_approx);
    jac.set(V49, V51, scratch);

    jac.set(V49, Cr48, 0.0);

    scratch = screened_rates(k_Cr50_to_p_V49_reaclib);
    jac.set(V49, Cr50, scratch);

    scratch = screened_rates(k_Cr52_p_to_He4_V49_reaclib)*Y(H1)*state.rho;
    jac.set(V49, Cr52, scratch);

    jac.set(V49, Cr54, 0.0);

    jac.set(V49, Mn51, 0.0);

    scratch = screened_rates(k_Mn53_to_He4_V49_reaclib);
    jac.set(V49, Mn53, scratch);

    jac.set(V49, Mn55, 0.0);

    jac.set(V49, Fe52, 0.0);

    jac.set(V49, Fe54, 0.0);

    jac.set(V49, Fe56, 0.0);

    jac.set(V49, Fe58, 0.0);

    jac.set(V49, Co55, 0.0);

    jac.set(V49, Co56, 0.0);

    jac.set(V49, Co57, 0.0);

    jac.set(V49, Co59, 0.0);

    jac.set(V49, Ni56, 0.0);

    jac.set(V49, Ni58, 0.0);

    jac.set(V49, Ni60, 0.0);

    jac.set(V49, Cu59, 0.0);

    jac.set(V49, Zn60, 0.0);


    // V51 row

    scratch = 1.0*screened_rates(k_V49_n_n_to_V51_approx)*Y(V49)*Y(N)*amrex::Math::powi<2>(state.rho)
              + 0.5*Y(V49)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_V49_n_n_to_V51_approx_dYN
              - Y(V51)*rate_eval.drate_V51_to_V49_n_n_approx_dYN;
    jac.set(V51, N, scratch);

    scratch = screened_rates(k_Cr54_p_to_He4_V51_reaclib)*Y(Cr54)*state.rho
              - screened_rates(k_V51_p_to_Cr52_reaclib)*Y(V51)*state.rho
              - screened_rates(k_V51_p_to_He4_Ti48_reaclib)*Y(V51)*state.rho;
    jac.set(V51, H1, scratch);

    scratch = screened_rates(k_Ti48_He4_to_p_V51_reaclib)*Y(Ti48)*state.rho
              - screened_rates(k_V51_He4_to_Mn55_reaclib)*Y(V51)*state.rho
              - screened_rates(k_V51_He4_to_p_Cr54_reaclib)*Y(V51)*state.rho;
    jac.set(V51, He4, scratch);

    jac.set(V51, C12, 0.0);

    jac.set(V51, O16, 0.0);

    jac.set(V51, Ne20, 0.0);

    jac.set(V51, Mg23, 0.0);

    jac.set(V51, Mg24, 0.0);

    jac.set(V51, Al27, 0.0);

    jac.set(V51, Si27, 0.0);

    jac.set(V51, Si28, 0.0);

    jac.set(V51, Si29, 0.0);

    jac.set(V51, Si30, 0.0);

    jac.set(V51, Si31, 0.0);

    jac.set(V51, P29, 0.0);

    jac.set(V51, P31, 0.0);

    jac.set(V51, P33, 0.0);

    jac.set(V51, S32, 0.0);

    jac.set(V51, S33, 0.0);

    jac.set(V51, S34, 0.0);

    jac.set(V51, S35, 0.0);

    jac.set(V51, Cl33, 0.0);

    jac.set(V51, Cl35, 0.0);

    jac.set(V51, Cl37, 0.0);

    jac.set(V51, Ar36, 0.0);

    jac.set(V51, Ar37, 0.0);

    jac.set(V51, Ar38, 0.0);

    jac.set(V51, K37, 0.0);

    jac.set(V51, K39, 0.0);

    jac.set(V51, Ca40, 0.0);

    jac.set(V51, Ca42, 0.0);

    jac.set(V51, Sc41, 0.0);

    jac.set(V51, Sc43, 0.0);

    jac.set(V51, Sc45, 0.0);

    jac.set(V51, Ti44, 0.0);

    jac.set(V51, Ti46, 0.0);

    scratch = screened_rates(k_Ti48_He4_to_p_V51_reaclib)*Y(He4)*state.rho;
    jac.set(V51, Ti48, scratch);

    jac.set(V51, V47, 0.0);

    scratch = 0.5*screened_rates(k_V49_n_n_to_V51_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(V51, V49, scratch);

    scratch = -screened_rates(k_V51_He4_to_Mn55_reaclib)*Y(He4)*state.rho
              - screened_rates(k_V51_He4_to_p_Cr54_reaclib)*Y(He4)*state.rho
              - screened_rates(k_V51_p_to_Cr52_reaclib)*Y(H1)*state.rho
              - screened_rates(k_V51_p_to_He4_Ti48_reaclib)*Y(H1)*state.rho
              - screened_rates(k_V51_to_V49_n_n_approx);
    jac.set(V51, V51, scratch);

    jac.set(V51, Cr48, 0.0);

    jac.set(V51, Cr50, 0.0);

    scratch = screened_rates(k_Cr52_to_p_V51_reaclib);
    jac.set(V51, Cr52, scratch);

    scratch = screened_rates(k_Cr54_p_to_He4_V51_reaclib)*Y(H1)*state.rho;
    jac.set(V51, Cr54, scratch);

    jac.set(V51, Mn51, 0.0);

    jac.set(V51, Mn53, 0.0);

    scratch = screened_rates(k_Mn55_to_He4_V51_reaclib);
    jac.set(V51, Mn55, scratch);

    jac.set(V51, Fe52, 0.0);

    jac.set(V51, Fe54, 0.0);

    jac.set(V51, Fe56, 0.0);

    jac.set(V51, Fe58, 0.0);

    jac.set(V51, Co55, 0.0);

    jac.set(V51, Co56, 0.0);

    jac.set(V51, Co57, 0.0);

    jac.set(V51, Co59, 0.0);

    jac.set(V51, Ni56, 0.0);

    jac.set(V51, Ni58, 0.0);

    jac.set(V51, Ni60, 0.0);

    jac.set(V51, Cu59, 0.0);

    jac.set(V51, Zn60, 0.0);


    // Cr48 row

    scratch = -1.0*screened_rates(k_Cr48_n_n_to_Cr50_approx)*Y(Cr48)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 0.5*Y(Cr48)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Cr48_n_n_to_Cr50_approx_dYN
              + Y(Cr50)*rate_eval.drate_Cr50_to_Cr48_n_n_approx_dYN;
    jac.set(Cr48, N, scratch);

    scratch = screened_rates(k_Mn51_p_to_He4_Cr48_reaclib)*Y(Mn51)*state.rho
              + screened_rates(k_V47_p_to_Cr48_reaclib)*Y(V47)*state.rho;
    jac.set(Cr48, H1, scratch);

    scratch = -screened_rates(k_Cr48_He4_to_Fe52_reaclib)*Y(Cr48)*state.rho
              - screened_rates(k_Cr48_He4_to_p_Mn51_reaclib)*Y(Cr48)*state.rho
              + screened_rates(k_Ti44_He4_to_Cr48_reaclib)*Y(Ti44)*state.rho;
    jac.set(Cr48, He4, scratch);

    jac.set(Cr48, C12, 0.0);

    jac.set(Cr48, O16, 0.0);

    jac.set(Cr48, Ne20, 0.0);

    jac.set(Cr48, Mg23, 0.0);

    jac.set(Cr48, Mg24, 0.0);

    jac.set(Cr48, Al27, 0.0);

    jac.set(Cr48, Si27, 0.0);

    jac.set(Cr48, Si28, 0.0);

    jac.set(Cr48, Si29, 0.0);

    jac.set(Cr48, Si30, 0.0);

    jac.set(Cr48, Si31, 0.0);

    jac.set(Cr48, P29, 0.0);

    jac.set(Cr48, P31, 0.0);

    jac.set(Cr48, P33, 0.0);

    jac.set(Cr48, S32, 0.0);

    jac.set(Cr48, S33, 0.0);

    jac.set(Cr48, S34, 0.0);

    jac.set(Cr48, S35, 0.0);

    jac.set(Cr48, Cl33, 0.0);

    jac.set(Cr48, Cl35, 0.0);

    jac.set(Cr48, Cl37, 0.0);

    jac.set(Cr48, Ar36, 0.0);

    jac.set(Cr48, Ar37, 0.0);

    jac.set(Cr48, Ar38, 0.0);

    jac.set(Cr48, K37, 0.0);

    jac.set(Cr48, K39, 0.0);

    jac.set(Cr48, Ca40, 0.0);

    jac.set(Cr48, Ca42, 0.0);

    jac.set(Cr48, Sc41, 0.0);

    jac.set(Cr48, Sc43, 0.0);

    jac.set(Cr48, Sc45, 0.0);

    scratch = screened_rates(k_Ti44_He4_to_Cr48_reaclib)*Y(He4)*state.rho;
    jac.set(Cr48, Ti44, scratch);

    jac.set(Cr48, Ti46, 0.0);

    jac.set(Cr48, Ti48, 0.0);

    scratch = screened_rates(k_V47_p_to_Cr48_reaclib)*Y(H1)*state.rho;
    jac.set(Cr48, V47, scratch);

    jac.set(Cr48, V49, 0.0);

    jac.set(Cr48, V51, 0.0);

    scratch = -screened_rates(k_Cr48_He4_to_Fe52_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Cr48_He4_to_p_Mn51_reaclib)*Y(He4)*state.rho
              - 0.5*screened_rates(k_Cr48_n_n_to_Cr50_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_Cr48_to_He4_Ti44_reaclib)
              - screened_rates(k_Cr48_to_p_V47_reaclib);
    jac.set(Cr48, Cr48, scratch);

    scratch = screened_rates(k_Cr50_to_Cr48_n_n_approx);
    jac.set(Cr48, Cr50, scratch);

    jac.set(Cr48, Cr52, 0.0);

    jac.set(Cr48, Cr54, 0.0);

    scratch = screened_rates(k_Mn51_p_to_He4_Cr48_reaclib)*Y(H1)*state.rho;
    jac.set(Cr48, Mn51, scratch);

    jac.set(Cr48, Mn53, 0.0);

    jac.set(Cr48, Mn55, 0.0);

    scratch = screened_rates(k_Fe52_to_He4_Cr48_reaclib);
    jac.set(Cr48, Fe52, scratch);

    jac.set(Cr48, Fe54, 0.0);

    jac.set(Cr48, Fe56, 0.0);

    jac.set(Cr48, Fe58, 0.0);

    jac.set(Cr48, Co55, 0.0);

    jac.set(Cr48, Co56, 0.0);

    jac.set(Cr48, Co57, 0.0);

    jac.set(Cr48, Co59, 0.0);

    jac.set(Cr48, Ni56, 0.0);

    jac.set(Cr48, Ni58, 0.0);

    jac.set(Cr48, Ni60, 0.0);

    jac.set(Cr48, Cu59, 0.0);

    jac.set(Cr48, Zn60, 0.0);


    // Cr50 row

    scratch = 1.0*screened_rates(k_Cr48_n_n_to_Cr50_approx)*Y(Cr48)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 1.0*screened_rates(k_Cr50_n_n_to_Cr52_approx)*Y(Cr50)*Y(N)*amrex::Math::powi<2>(state.rho)
              + 0.5*Y(Cr48)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Cr48_n_n_to_Cr50_approx_dYN
              - 0.5*Y(Cr50)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Cr50_n_n_to_Cr52_approx_dYN
              - Y(Cr50)*rate_eval.drate_Cr50_to_Cr48_n_n_approx_dYN
              + Y(Cr52)*rate_eval.drate_Cr52_to_Cr50_n_n_approx_dYN;
    jac.set(Cr50, N, scratch);

    scratch = -screened_rates(k_Cr50_p_to_He4_V47_reaclib)*Y(Cr50)*state.rho
              - screened_rates(k_Cr50_p_to_Mn51_reaclib)*Y(Cr50)*state.rho
              + screened_rates(k_Mn53_p_to_He4_Cr50_reaclib)*Y(Mn53)*state.rho
              + screened_rates(k_V49_p_to_Cr50_reaclib)*Y(V49)*state.rho;
    jac.set(Cr50, H1, scratch);

    scratch = -screened_rates(k_Cr50_He4_to_Fe54_reaclib)*Y(Cr50)*state.rho
              - screened_rates(k_Cr50_He4_to_p_Mn53_reaclib)*Y(Cr50)*state.rho
              + screened_rates(k_Ti46_He4_to_Cr50_reaclib)*Y(Ti46)*state.rho
              + screened_rates(k_V47_He4_to_p_Cr50_reaclib)*Y(V47)*state.rho;
    jac.set(Cr50, He4, scratch);

    jac.set(Cr50, C12, 0.0);

    jac.set(Cr50, O16, 0.0);

    jac.set(Cr50, Ne20, 0.0);

    jac.set(Cr50, Mg23, 0.0);

    jac.set(Cr50, Mg24, 0.0);

    jac.set(Cr50, Al27, 0.0);

    jac.set(Cr50, Si27, 0.0);

    jac.set(Cr50, Si28, 0.0);

    jac.set(Cr50, Si29, 0.0);

    jac.set(Cr50, Si30, 0.0);

    jac.set(Cr50, Si31, 0.0);

    jac.set(Cr50, P29, 0.0);

    jac.set(Cr50, P31, 0.0);

    jac.set(Cr50, P33, 0.0);

    jac.set(Cr50, S32, 0.0);

    jac.set(Cr50, S33, 0.0);

    jac.set(Cr50, S34, 0.0);

    jac.set(Cr50, S35, 0.0);

    jac.set(Cr50, Cl33, 0.0);

    jac.set(Cr50, Cl35, 0.0);

    jac.set(Cr50, Cl37, 0.0);

    jac.set(Cr50, Ar36, 0.0);

    jac.set(Cr50, Ar37, 0.0);

    jac.set(Cr50, Ar38, 0.0);

    jac.set(Cr50, K37, 0.0);

    jac.set(Cr50, K39, 0.0);

    jac.set(Cr50, Ca40, 0.0);

    jac.set(Cr50, Ca42, 0.0);

    jac.set(Cr50, Sc41, 0.0);

    jac.set(Cr50, Sc43, 0.0);

    jac.set(Cr50, Sc45, 0.0);

    jac.set(Cr50, Ti44, 0.0);

    scratch = screened_rates(k_Ti46_He4_to_Cr50_reaclib)*Y(He4)*state.rho;
    jac.set(Cr50, Ti46, scratch);

    jac.set(Cr50, Ti48, 0.0);

    scratch = screened_rates(k_V47_He4_to_p_Cr50_reaclib)*Y(He4)*state.rho;
    jac.set(Cr50, V47, scratch);

    scratch = screened_rates(k_V49_p_to_Cr50_reaclib)*Y(H1)*state.rho;
    jac.set(Cr50, V49, scratch);

    jac.set(Cr50, V51, 0.0);

    scratch = 0.5*screened_rates(k_Cr48_n_n_to_Cr50_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(Cr50, Cr48, scratch);

    scratch = -screened_rates(k_Cr50_He4_to_Fe54_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Cr50_He4_to_p_Mn53_reaclib)*Y(He4)*state.rho
              - 0.5*screened_rates(k_Cr50_n_n_to_Cr52_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_Cr50_p_to_He4_V47_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Cr50_p_to_Mn51_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Cr50_to_Cr48_n_n_approx)
              - screened_rates(k_Cr50_to_He4_Ti46_reaclib)
              - screened_rates(k_Cr50_to_p_V49_reaclib);
    jac.set(Cr50, Cr50, scratch);

    scratch = screened_rates(k_Cr52_to_Cr50_n_n_approx);
    jac.set(Cr50, Cr52, scratch);

    jac.set(Cr50, Cr54, 0.0);

    scratch = screened_rates(k_Mn51_to_p_Cr50_reaclib);
    jac.set(Cr50, Mn51, scratch);

    scratch = screened_rates(k_Mn53_p_to_He4_Cr50_reaclib)*Y(H1)*state.rho;
    jac.set(Cr50, Mn53, scratch);

    jac.set(Cr50, Mn55, 0.0);

    jac.set(Cr50, Fe52, 0.0);

    scratch = screened_rates(k_Fe54_to_He4_Cr50_reaclib);
    jac.set(Cr50, Fe54, scratch);

    jac.set(Cr50, Fe56, 0.0);

    jac.set(Cr50, Fe58, 0.0);

    jac.set(Cr50, Co55, 0.0);

    jac.set(Cr50, Co56, 0.0);

    jac.set(Cr50, Co57, 0.0);

    jac.set(Cr50, Co59, 0.0);

    jac.set(Cr50, Ni56, 0.0);

    jac.set(Cr50, Ni58, 0.0);

    jac.set(Cr50, Ni60, 0.0);

    jac.set(Cr50, Cu59, 0.0);

    jac.set(Cr50, Zn60, 0.0);


    // Cr52 row

    scratch = 1.0*screened_rates(k_Cr50_n_n_to_Cr52_approx)*Y(Cr50)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 1.0*screened_rates(k_Cr52_n_n_to_Cr54_approx)*Y(Cr52)*Y(N)*amrex::Math::powi<2>(state.rho)
              + 0.5*Y(Cr50)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Cr50_n_n_to_Cr52_approx_dYN
              - 0.5*Y(Cr52)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Cr52_n_n_to_Cr54_approx_dYN
              - Y(Cr52)*rate_eval.drate_Cr52_to_Cr50_n_n_approx_dYN
              + Y(Cr54)*rate_eval.drate_Cr54_to_Cr52_n_n_approx_dYN;
    jac.set(Cr52, N, scratch);

    scratch = -screened_rates(k_Cr52_p_to_He4_V49_reaclib)*Y(Cr52)*state.rho
              - screened_rates(k_Cr52_p_to_Mn53_reaclib)*Y(Cr52)*state.rho
              + screened_rates(k_Mn55_p_to_He4_Cr52_reaclib)*Y(Mn55)*state.rho
              + screened_rates(k_V51_p_to_Cr52_reaclib)*Y(V51)*state.rho;
    jac.set(Cr52, H1, scratch);

    scratch = -screened_rates(k_Cr52_He4_to_Fe56_reaclib)*Y(Cr52)*state.rho
              - screened_rates(k_Cr52_He4_to_p_Mn55_reaclib)*Y(Cr52)*state.rho
              + screened_rates(k_Ti48_He4_to_Cr52_reaclib)*Y(Ti48)*state.rho
              + screened_rates(k_V49_He4_to_p_Cr52_reaclib)*Y(V49)*state.rho;
    jac.set(Cr52, He4, scratch);

    jac.set(Cr52, C12, 0.0);

    jac.set(Cr52, O16, 0.0);

    jac.set(Cr52, Ne20, 0.0);

    jac.set(Cr52, Mg23, 0.0);

    jac.set(Cr52, Mg24, 0.0);

    jac.set(Cr52, Al27, 0.0);

    jac.set(Cr52, Si27, 0.0);

    jac.set(Cr52, Si28, 0.0);

    jac.set(Cr52, Si29, 0.0);

    jac.set(Cr52, Si30, 0.0);

    jac.set(Cr52, Si31, 0.0);

    jac.set(Cr52, P29, 0.0);

    jac.set(Cr52, P31, 0.0);

    jac.set(Cr52, P33, 0.0);

    jac.set(Cr52, S32, 0.0);

    jac.set(Cr52, S33, 0.0);

    jac.set(Cr52, S34, 0.0);

    jac.set(Cr52, S35, 0.0);

    jac.set(Cr52, Cl33, 0.0);

    jac.set(Cr52, Cl35, 0.0);

    jac.set(Cr52, Cl37, 0.0);

    jac.set(Cr52, Ar36, 0.0);

    jac.set(Cr52, Ar37, 0.0);

    jac.set(Cr52, Ar38, 0.0);

    jac.set(Cr52, K37, 0.0);

    jac.set(Cr52, K39, 0.0);

    jac.set(Cr52, Ca40, 0.0);

    jac.set(Cr52, Ca42, 0.0);

    jac.set(Cr52, Sc41, 0.0);

    jac.set(Cr52, Sc43, 0.0);

    jac.set(Cr52, Sc45, 0.0);

    jac.set(Cr52, Ti44, 0.0);

    jac.set(Cr52, Ti46, 0.0);

    scratch = screened_rates(k_Ti48_He4_to_Cr52_reaclib)*Y(He4)*state.rho;
    jac.set(Cr52, Ti48, scratch);

    jac.set(Cr52, V47, 0.0);

    scratch = screened_rates(k_V49_He4_to_p_Cr52_reaclib)*Y(He4)*state.rho;
    jac.set(Cr52, V49, scratch);

    scratch = screened_rates(k_V51_p_to_Cr52_reaclib)*Y(H1)*state.rho;
    jac.set(Cr52, V51, scratch);

    jac.set(Cr52, Cr48, 0.0);

    scratch = 0.5*screened_rates(k_Cr50_n_n_to_Cr52_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(Cr52, Cr50, scratch);

    scratch = -screened_rates(k_Cr52_He4_to_Fe56_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Cr52_He4_to_p_Mn55_reaclib)*Y(He4)*state.rho
              - 0.5*screened_rates(k_Cr52_n_n_to_Cr54_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_Cr52_p_to_He4_V49_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Cr52_p_to_Mn53_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Cr52_to_Cr50_n_n_approx)
              - screened_rates(k_Cr52_to_He4_Ti48_reaclib)
              - screened_rates(k_Cr52_to_p_V51_reaclib);
    jac.set(Cr52, Cr52, scratch);

    scratch = screened_rates(k_Cr54_to_Cr52_n_n_approx);
    jac.set(Cr52, Cr54, scratch);

    jac.set(Cr52, Mn51, 0.0);

    scratch = screened_rates(k_Mn53_to_p_Cr52_reaclib);
    jac.set(Cr52, Mn53, scratch);

    scratch = screened_rates(k_Mn55_p_to_He4_Cr52_reaclib)*Y(H1)*state.rho;
    jac.set(Cr52, Mn55, scratch);

    jac.set(Cr52, Fe52, 0.0);

    jac.set(Cr52, Fe54, 0.0);

    scratch = screened_rates(k_Fe56_to_He4_Cr52_reaclib);
    jac.set(Cr52, Fe56, scratch);

    jac.set(Cr52, Fe58, 0.0);

    jac.set(Cr52, Co55, 0.0);

    jac.set(Cr52, Co56, 0.0);

    jac.set(Cr52, Co57, 0.0);

    jac.set(Cr52, Co59, 0.0);

    jac.set(Cr52, Ni56, 0.0);

    jac.set(Cr52, Ni58, 0.0);

    jac.set(Cr52, Ni60, 0.0);

    jac.set(Cr52, Cu59, 0.0);

    jac.set(Cr52, Zn60, 0.0);


    // Cr54 row

    scratch = 1.0*screened_rates(k_Cr52_n_n_to_Cr54_approx)*Y(Cr52)*Y(N)*amrex::Math::powi<2>(state.rho)
              + 0.5*Y(Cr52)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Cr52_n_n_to_Cr54_approx_dYN
              - Y(Cr54)*rate_eval.drate_Cr54_to_Cr52_n_n_approx_dYN;
    jac.set(Cr54, N, scratch);

    scratch = -screened_rates(k_Cr54_p_to_He4_V51_reaclib)*Y(Cr54)*state.rho
              - screened_rates(k_Cr54_p_to_Mn55_reaclib)*Y(Cr54)*state.rho;
    jac.set(Cr54, H1, scratch);

    scratch = -screened_rates(k_Cr54_He4_to_Fe58_reaclib)*Y(Cr54)*state.rho
              + screened_rates(k_V51_He4_to_p_Cr54_reaclib)*Y(V51)*state.rho;
    jac.set(Cr54, He4, scratch);

    jac.set(Cr54, C12, 0.0);

    jac.set(Cr54, O16, 0.0);

    jac.set(Cr54, Ne20, 0.0);

    jac.set(Cr54, Mg23, 0.0);

    jac.set(Cr54, Mg24, 0.0);

    jac.set(Cr54, Al27, 0.0);

    jac.set(Cr54, Si27, 0.0);

    jac.set(Cr54, Si28, 0.0);

    jac.set(Cr54, Si29, 0.0);

    jac.set(Cr54, Si30, 0.0);

    jac.set(Cr54, Si31, 0.0);

    jac.set(Cr54, P29, 0.0);

    jac.set(Cr54, P31, 0.0);

    jac.set(Cr54, P33, 0.0);

    jac.set(Cr54, S32, 0.0);

    jac.set(Cr54, S33, 0.0);

    jac.set(Cr54, S34, 0.0);

    jac.set(Cr54, S35, 0.0);

    jac.set(Cr54, Cl33, 0.0);

    jac.set(Cr54, Cl35, 0.0);

    jac.set(Cr54, Cl37, 0.0);

    jac.set(Cr54, Ar36, 0.0);

    jac.set(Cr54, Ar37, 0.0);

    jac.set(Cr54, Ar38, 0.0);

    jac.set(Cr54, K37, 0.0);

    jac.set(Cr54, K39, 0.0);

    jac.set(Cr54, Ca40, 0.0);

    jac.set(Cr54, Ca42, 0.0);

    jac.set(Cr54, Sc41, 0.0);

    jac.set(Cr54, Sc43, 0.0);

    jac.set(Cr54, Sc45, 0.0);

    jac.set(Cr54, Ti44, 0.0);

    jac.set(Cr54, Ti46, 0.0);

    jac.set(Cr54, Ti48, 0.0);

    jac.set(Cr54, V47, 0.0);

    jac.set(Cr54, V49, 0.0);

    scratch = screened_rates(k_V51_He4_to_p_Cr54_reaclib)*Y(He4)*state.rho;
    jac.set(Cr54, V51, scratch);

    jac.set(Cr54, Cr48, 0.0);

    jac.set(Cr54, Cr50, 0.0);

    scratch = 0.5*screened_rates(k_Cr52_n_n_to_Cr54_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(Cr54, Cr52, scratch);

    scratch = -screened_rates(k_Cr54_He4_to_Fe58_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Cr54_p_to_He4_V51_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Cr54_p_to_Mn55_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Cr54_to_Cr52_n_n_approx);
    jac.set(Cr54, Cr54, scratch);

    jac.set(Cr54, Mn51, 0.0);

    jac.set(Cr54, Mn53, 0.0);

    scratch = screened_rates(k_Mn55_to_p_Cr54_reaclib);
    jac.set(Cr54, Mn55, scratch);

    jac.set(Cr54, Fe52, 0.0);

    jac.set(Cr54, Fe54, 0.0);

    jac.set(Cr54, Fe56, 0.0);

    scratch = screened_rates(k_Fe58_to_He4_Cr54_reaclib);
    jac.set(Cr54, Fe58, scratch);

    jac.set(Cr54, Co55, 0.0);

    jac.set(Cr54, Co56, 0.0);

    jac.set(Cr54, Co57, 0.0);

    jac.set(Cr54, Co59, 0.0);

    jac.set(Cr54, Ni56, 0.0);

    jac.set(Cr54, Ni58, 0.0);

    jac.set(Cr54, Ni60, 0.0);

    jac.set(Cr54, Cu59, 0.0);

    jac.set(Cr54, Zn60, 0.0);


    // Mn51 row

    scratch = -1.0*screened_rates(k_Mn51_n_n_to_Mn53_approx)*Y(Mn51)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 0.5*Y(Mn51)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Mn51_n_n_to_Mn53_approx_dYN
              + Y(Mn53)*rate_eval.drate_Mn53_to_Mn51_n_n_approx_dYN;
    jac.set(Mn51, N, scratch);

    scratch = screened_rates(k_Cr50_p_to_Mn51_reaclib)*Y(Cr50)*state.rho
              + screened_rates(k_Fe54_p_to_He4_Mn51_reaclib)*Y(Fe54)*state.rho
              - screened_rates(k_Mn51_p_to_Fe52_reaclib)*Y(Mn51)*state.rho
              - screened_rates(k_Mn51_p_to_He4_Cr48_reaclib)*Y(Mn51)*state.rho;
    jac.set(Mn51, H1, scratch);

    scratch = screened_rates(k_Cr48_He4_to_p_Mn51_reaclib)*Y(Cr48)*state.rho
              - screened_rates(k_Mn51_He4_to_Co55_reaclib)*Y(Mn51)*state.rho
              - screened_rates(k_Mn51_He4_to_p_Fe54_reaclib)*Y(Mn51)*state.rho
              + screened_rates(k_V47_He4_to_Mn51_reaclib)*Y(V47)*state.rho;
    jac.set(Mn51, He4, scratch);

    jac.set(Mn51, C12, 0.0);

    jac.set(Mn51, O16, 0.0);

    jac.set(Mn51, Ne20, 0.0);

    jac.set(Mn51, Mg23, 0.0);

    jac.set(Mn51, Mg24, 0.0);

    jac.set(Mn51, Al27, 0.0);

    jac.set(Mn51, Si27, 0.0);

    jac.set(Mn51, Si28, 0.0);

    jac.set(Mn51, Si29, 0.0);

    jac.set(Mn51, Si30, 0.0);

    jac.set(Mn51, Si31, 0.0);

    jac.set(Mn51, P29, 0.0);

    jac.set(Mn51, P31, 0.0);

    jac.set(Mn51, P33, 0.0);

    jac.set(Mn51, S32, 0.0);

    jac.set(Mn51, S33, 0.0);

    jac.set(Mn51, S34, 0.0);

    jac.set(Mn51, S35, 0.0);

    jac.set(Mn51, Cl33, 0.0);

    jac.set(Mn51, Cl35, 0.0);

    jac.set(Mn51, Cl37, 0.0);

    jac.set(Mn51, Ar36, 0.0);

    jac.set(Mn51, Ar37, 0.0);

    jac.set(Mn51, Ar38, 0.0);

    jac.set(Mn51, K37, 0.0);

    jac.set(Mn51, K39, 0.0);

    jac.set(Mn51, Ca40, 0.0);

    jac.set(Mn51, Ca42, 0.0);

    jac.set(Mn51, Sc41, 0.0);

    jac.set(Mn51, Sc43, 0.0);

    jac.set(Mn51, Sc45, 0.0);

    jac.set(Mn51, Ti44, 0.0);

    jac.set(Mn51, Ti46, 0.0);

    jac.set(Mn51, Ti48, 0.0);

    scratch = screened_rates(k_V47_He4_to_Mn51_reaclib)*Y(He4)*state.rho;
    jac.set(Mn51, V47, scratch);

    jac.set(Mn51, V49, 0.0);

    jac.set(Mn51, V51, 0.0);

    scratch = screened_rates(k_Cr48_He4_to_p_Mn51_reaclib)*Y(He4)*state.rho;
    jac.set(Mn51, Cr48, scratch);

    scratch = screened_rates(k_Cr50_p_to_Mn51_reaclib)*Y(H1)*state.rho;
    jac.set(Mn51, Cr50, scratch);

    jac.set(Mn51, Cr52, 0.0);

    jac.set(Mn51, Cr54, 0.0);

    scratch = -screened_rates(k_Mn51_He4_to_Co55_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Mn51_He4_to_p_Fe54_reaclib)*Y(He4)*state.rho
              - 0.5*screened_rates(k_Mn51_n_n_to_Mn53_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_Mn51_p_to_Fe52_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Mn51_p_to_He4_Cr48_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Mn51_to_He4_V47_reaclib)
              - screened_rates(k_Mn51_to_p_Cr50_reaclib);
    jac.set(Mn51, Mn51, scratch);

    scratch = screened_rates(k_Mn53_to_Mn51_n_n_approx);
    jac.set(Mn51, Mn53, scratch);

    jac.set(Mn51, Mn55, 0.0);

    scratch = screened_rates(k_Fe52_to_p_Mn51_reaclib);
    jac.set(Mn51, Fe52, scratch);

    scratch = screened_rates(k_Fe54_p_to_He4_Mn51_reaclib)*Y(H1)*state.rho;
    jac.set(Mn51, Fe54, scratch);

    jac.set(Mn51, Fe56, 0.0);

    jac.set(Mn51, Fe58, 0.0);

    scratch = screened_rates(k_Co55_to_He4_Mn51_reaclib);
    jac.set(Mn51, Co55, scratch);

    jac.set(Mn51, Co56, 0.0);

    jac.set(Mn51, Co57, 0.0);

    jac.set(Mn51, Co59, 0.0);

    jac.set(Mn51, Ni56, 0.0);

    jac.set(Mn51, Ni58, 0.0);

    jac.set(Mn51, Ni60, 0.0);

    jac.set(Mn51, Cu59, 0.0);

    jac.set(Mn51, Zn60, 0.0);


    // Mn53 row

    scratch = screened_rates(k_Co56_n_to_He4_Mn53_reaclib)*Y(Co56)*state.rho
              + 1.0*screened_rates(k_Mn51_n_n_to_Mn53_approx)*Y(Mn51)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 1.0*screened_rates(k_Mn53_n_n_to_Mn55_approx)*Y(Mn53)*Y(N)*amrex::Math::powi<2>(state.rho)
              + 0.5*Y(Mn51)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Mn51_n_n_to_Mn53_approx_dYN
              - 0.5*Y(Mn53)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Mn53_n_n_to_Mn55_approx_dYN
              - Y(Mn53)*rate_eval.drate_Mn53_to_Mn51_n_n_approx_dYN
              + Y(Mn55)*rate_eval.drate_Mn55_to_Mn53_n_n_approx_dYN;
    jac.set(Mn53, N, scratch);

    scratch = screened_rates(k_Cr52_p_to_Mn53_reaclib)*Y(Cr52)*state.rho
              + screened_rates(k_Fe56_p_to_He4_Mn53_reaclib)*Y(Fe56)*state.rho
              - screened_rates(k_Mn53_p_to_Fe54_reaclib)*Y(Mn53)*state.rho
              - screened_rates(k_Mn53_p_to_He4_Cr50_reaclib)*Y(Mn53)*state.rho;
    jac.set(Mn53, H1, scratch);

    scratch = screened_rates(k_Cr50_He4_to_p_Mn53_reaclib)*Y(Cr50)*state.rho
              - screened_rates(k_Mn53_He4_to_Co57_reaclib)*Y(Mn53)*state.rho
              - screened_rates(k_Mn53_He4_to_n_Co56_reaclib)*Y(Mn53)*state.rho
              - screened_rates(k_Mn53_He4_to_p_Fe56_reaclib)*Y(Mn53)*state.rho
              + screened_rates(k_V49_He4_to_Mn53_reaclib)*Y(V49)*state.rho;
    jac.set(Mn53, He4, scratch);

    jac.set(Mn53, C12, 0.0);

    jac.set(Mn53, O16, 0.0);

    jac.set(Mn53, Ne20, 0.0);

    jac.set(Mn53, Mg23, 0.0);

    jac.set(Mn53, Mg24, 0.0);

    jac.set(Mn53, Al27, 0.0);

    jac.set(Mn53, Si27, 0.0);

    jac.set(Mn53, Si28, 0.0);

    jac.set(Mn53, Si29, 0.0);

    jac.set(Mn53, Si30, 0.0);

    jac.set(Mn53, Si31, 0.0);

    jac.set(Mn53, P29, 0.0);

    jac.set(Mn53, P31, 0.0);

    jac.set(Mn53, P33, 0.0);

    jac.set(Mn53, S32, 0.0);

    jac.set(Mn53, S33, 0.0);

    jac.set(Mn53, S34, 0.0);

    jac.set(Mn53, S35, 0.0);

    jac.set(Mn53, Cl33, 0.0);

    jac.set(Mn53, Cl35, 0.0);

    jac.set(Mn53, Cl37, 0.0);

    jac.set(Mn53, Ar36, 0.0);

    jac.set(Mn53, Ar37, 0.0);

    jac.set(Mn53, Ar38, 0.0);

    jac.set(Mn53, K37, 0.0);

    jac.set(Mn53, K39, 0.0);

    jac.set(Mn53, Ca40, 0.0);

    jac.set(Mn53, Ca42, 0.0);

    jac.set(Mn53, Sc41, 0.0);

    jac.set(Mn53, Sc43, 0.0);

    jac.set(Mn53, Sc45, 0.0);

    jac.set(Mn53, Ti44, 0.0);

    jac.set(Mn53, Ti46, 0.0);

    jac.set(Mn53, Ti48, 0.0);

    jac.set(Mn53, V47, 0.0);

    scratch = screened_rates(k_V49_He4_to_Mn53_reaclib)*Y(He4)*state.rho;
    jac.set(Mn53, V49, scratch);

    jac.set(Mn53, V51, 0.0);

    jac.set(Mn53, Cr48, 0.0);

    scratch = screened_rates(k_Cr50_He4_to_p_Mn53_reaclib)*Y(He4)*state.rho;
    jac.set(Mn53, Cr50, scratch);

    scratch = screened_rates(k_Cr52_p_to_Mn53_reaclib)*Y(H1)*state.rho;
    jac.set(Mn53, Cr52, scratch);

    jac.set(Mn53, Cr54, 0.0);

    scratch = 0.5*screened_rates(k_Mn51_n_n_to_Mn53_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(Mn53, Mn51, scratch);

    scratch = -screened_rates(k_Mn53_He4_to_Co57_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Mn53_He4_to_n_Co56_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Mn53_He4_to_p_Fe56_reaclib)*Y(He4)*state.rho
              - 0.5*screened_rates(k_Mn53_n_n_to_Mn55_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_Mn53_p_to_Fe54_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Mn53_p_to_He4_Cr50_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Mn53_to_He4_V49_reaclib)
              - screened_rates(k_Mn53_to_Mn51_n_n_approx)
              - screened_rates(k_Mn53_to_p_Cr52_reaclib);
    jac.set(Mn53, Mn53, scratch);

    scratch = screened_rates(k_Mn55_to_Mn53_n_n_approx);
    jac.set(Mn53, Mn55, scratch);

    jac.set(Mn53, Fe52, 0.0);

    scratch = screened_rates(k_Fe54_to_p_Mn53_reaclib);
    jac.set(Mn53, Fe54, scratch);

    scratch = screened_rates(k_Fe56_p_to_He4_Mn53_reaclib)*Y(H1)*state.rho;
    jac.set(Mn53, Fe56, scratch);

    jac.set(Mn53, Fe58, 0.0);

    jac.set(Mn53, Co55, 0.0);

    scratch = screened_rates(k_Co56_n_to_He4_Mn53_reaclib)*Y(N)*state.rho;
    jac.set(Mn53, Co56, scratch);

    scratch = screened_rates(k_Co57_to_He4_Mn53_reaclib);
    jac.set(Mn53, Co57, scratch);

    jac.set(Mn53, Co59, 0.0);

    jac.set(Mn53, Ni56, 0.0);

    jac.set(Mn53, Ni58, 0.0);

    jac.set(Mn53, Ni60, 0.0);

    jac.set(Mn53, Cu59, 0.0);

    jac.set(Mn53, Zn60, 0.0);


    // Mn55 row

    scratch = 1.0*screened_rates(k_Mn53_n_n_to_Mn55_approx)*Y(Mn53)*Y(N)*amrex::Math::powi<2>(state.rho)
              + 0.5*Y(Mn53)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Mn53_n_n_to_Mn55_approx_dYN
              - Y(Mn55)*rate_eval.drate_Mn55_to_Mn53_n_n_approx_dYN;
    jac.set(Mn55, N, scratch);

    scratch = screened_rates(k_Cr54_p_to_Mn55_reaclib)*Y(Cr54)*state.rho
              + screened_rates(k_Fe58_p_to_He4_Mn55_reaclib)*Y(Fe58)*state.rho
              - screened_rates(k_Mn55_p_to_Fe56_reaclib)*Y(Mn55)*state.rho
              - screened_rates(k_Mn55_p_to_He4_Cr52_reaclib)*Y(Mn55)*state.rho;
    jac.set(Mn55, H1, scratch);

    scratch = screened_rates(k_Cr52_He4_to_p_Mn55_reaclib)*Y(Cr52)*state.rho
              - screened_rates(k_Mn55_He4_to_Co59_reaclib)*Y(Mn55)*state.rho
              - screened_rates(k_Mn55_He4_to_p_Fe58_reaclib)*Y(Mn55)*state.rho
              + screened_rates(k_V51_He4_to_Mn55_reaclib)*Y(V51)*state.rho;
    jac.set(Mn55, He4, scratch);

    jac.set(Mn55, C12, 0.0);

    jac.set(Mn55, O16, 0.0);

    jac.set(Mn55, Ne20, 0.0);

    jac.set(Mn55, Mg23, 0.0);

    jac.set(Mn55, Mg24, 0.0);

    jac.set(Mn55, Al27, 0.0);

    jac.set(Mn55, Si27, 0.0);

    jac.set(Mn55, Si28, 0.0);

    jac.set(Mn55, Si29, 0.0);

    jac.set(Mn55, Si30, 0.0);

    jac.set(Mn55, Si31, 0.0);

    jac.set(Mn55, P29, 0.0);

    jac.set(Mn55, P31, 0.0);

    jac.set(Mn55, P33, 0.0);

    jac.set(Mn55, S32, 0.0);

    jac.set(Mn55, S33, 0.0);

    jac.set(Mn55, S34, 0.0);

    jac.set(Mn55, S35, 0.0);

    jac.set(Mn55, Cl33, 0.0);

    jac.set(Mn55, Cl35, 0.0);

    jac.set(Mn55, Cl37, 0.0);

    jac.set(Mn55, Ar36, 0.0);

    jac.set(Mn55, Ar37, 0.0);

    jac.set(Mn55, Ar38, 0.0);

    jac.set(Mn55, K37, 0.0);

    jac.set(Mn55, K39, 0.0);

    jac.set(Mn55, Ca40, 0.0);

    jac.set(Mn55, Ca42, 0.0);

    jac.set(Mn55, Sc41, 0.0);

    jac.set(Mn55, Sc43, 0.0);

    jac.set(Mn55, Sc45, 0.0);

    jac.set(Mn55, Ti44, 0.0);

    jac.set(Mn55, Ti46, 0.0);

    jac.set(Mn55, Ti48, 0.0);

    jac.set(Mn55, V47, 0.0);

    jac.set(Mn55, V49, 0.0);

    scratch = screened_rates(k_V51_He4_to_Mn55_reaclib)*Y(He4)*state.rho;
    jac.set(Mn55, V51, scratch);

    jac.set(Mn55, Cr48, 0.0);

    jac.set(Mn55, Cr50, 0.0);

    scratch = screened_rates(k_Cr52_He4_to_p_Mn55_reaclib)*Y(He4)*state.rho;
    jac.set(Mn55, Cr52, scratch);

    scratch = screened_rates(k_Cr54_p_to_Mn55_reaclib)*Y(H1)*state.rho;
    jac.set(Mn55, Cr54, scratch);

    jac.set(Mn55, Mn51, 0.0);

    scratch = 0.5*screened_rates(k_Mn53_n_n_to_Mn55_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(Mn55, Mn53, scratch);

    scratch = -screened_rates(k_Mn55_He4_to_Co59_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Mn55_He4_to_p_Fe58_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Mn55_p_to_Fe56_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Mn55_p_to_He4_Cr52_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Mn55_to_He4_V51_reaclib)
              - screened_rates(k_Mn55_to_Mn53_n_n_approx)
              - screened_rates(k_Mn55_to_p_Cr54_reaclib);
    jac.set(Mn55, Mn55, scratch);

    jac.set(Mn55, Fe52, 0.0);

    jac.set(Mn55, Fe54, 0.0);

    scratch = screened_rates(k_Fe56_to_p_Mn55_reaclib);
    jac.set(Mn55, Fe56, scratch);

    scratch = screened_rates(k_Fe58_p_to_He4_Mn55_reaclib)*Y(H1)*state.rho;
    jac.set(Mn55, Fe58, scratch);

    jac.set(Mn55, Co55, 0.0);

    jac.set(Mn55, Co56, 0.0);

    jac.set(Mn55, Co57, 0.0);

    scratch = screened_rates(k_Co59_to_He4_Mn55_reaclib);
    jac.set(Mn55, Co59, scratch);

    jac.set(Mn55, Ni56, 0.0);

    jac.set(Mn55, Ni58, 0.0);

    jac.set(Mn55, Ni60, 0.0);

    jac.set(Mn55, Cu59, 0.0);

    jac.set(Mn55, Zn60, 0.0);


    // Fe52 row

    scratch = -1.0*screened_rates(k_Fe52_n_n_to_Fe54_approx)*Y(Fe52)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 0.5*Y(Fe52)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Fe52_n_n_to_Fe54_approx_dYN
              + Y(Fe54)*rate_eval.drate_Fe54_to_Fe52_n_n_approx_dYN;
    jac.set(Fe52, N, scratch);

    scratch = screened_rates(k_Co55_p_to_He4_Fe52_reaclib)*Y(Co55)*state.rho
              + screened_rates(k_Mn51_p_to_Fe52_reaclib)*Y(Mn51)*state.rho;
    jac.set(Fe52, H1, scratch);

    scratch = screened_rates(k_Cr48_He4_to_Fe52_reaclib)*Y(Cr48)*state.rho
              - screened_rates(k_Fe52_He4_to_Ni56_reaclib)*Y(Fe52)*state.rho
              - screened_rates(k_Fe52_He4_to_p_Co55_reaclib)*Y(Fe52)*state.rho;
    jac.set(Fe52, He4, scratch);

    jac.set(Fe52, C12, 0.0);

    jac.set(Fe52, O16, 0.0);

    jac.set(Fe52, Ne20, 0.0);

    jac.set(Fe52, Mg23, 0.0);

    jac.set(Fe52, Mg24, 0.0);

    jac.set(Fe52, Al27, 0.0);

    jac.set(Fe52, Si27, 0.0);

    jac.set(Fe52, Si28, 0.0);

    jac.set(Fe52, Si29, 0.0);

    jac.set(Fe52, Si30, 0.0);

    jac.set(Fe52, Si31, 0.0);

    jac.set(Fe52, P29, 0.0);

    jac.set(Fe52, P31, 0.0);

    jac.set(Fe52, P33, 0.0);

    jac.set(Fe52, S32, 0.0);

    jac.set(Fe52, S33, 0.0);

    jac.set(Fe52, S34, 0.0);

    jac.set(Fe52, S35, 0.0);

    jac.set(Fe52, Cl33, 0.0);

    jac.set(Fe52, Cl35, 0.0);

    jac.set(Fe52, Cl37, 0.0);

    jac.set(Fe52, Ar36, 0.0);

    jac.set(Fe52, Ar37, 0.0);

    jac.set(Fe52, Ar38, 0.0);

    jac.set(Fe52, K37, 0.0);

    jac.set(Fe52, K39, 0.0);

    jac.set(Fe52, Ca40, 0.0);

    jac.set(Fe52, Ca42, 0.0);

    jac.set(Fe52, Sc41, 0.0);

    jac.set(Fe52, Sc43, 0.0);

    jac.set(Fe52, Sc45, 0.0);

    jac.set(Fe52, Ti44, 0.0);

    jac.set(Fe52, Ti46, 0.0);

    jac.set(Fe52, Ti48, 0.0);

    jac.set(Fe52, V47, 0.0);

    jac.set(Fe52, V49, 0.0);

    jac.set(Fe52, V51, 0.0);

    scratch = screened_rates(k_Cr48_He4_to_Fe52_reaclib)*Y(He4)*state.rho;
    jac.set(Fe52, Cr48, scratch);

    jac.set(Fe52, Cr50, 0.0);

    jac.set(Fe52, Cr52, 0.0);

    jac.set(Fe52, Cr54, 0.0);

    scratch = screened_rates(k_Mn51_p_to_Fe52_reaclib)*Y(H1)*state.rho;
    jac.set(Fe52, Mn51, scratch);

    jac.set(Fe52, Mn53, 0.0);

    jac.set(Fe52, Mn55, 0.0);

    scratch = -screened_rates(k_Fe52_He4_to_Ni56_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Fe52_He4_to_p_Co55_reaclib)*Y(He4)*state.rho
              - 0.5*screened_rates(k_Fe52_n_n_to_Fe54_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_Fe52_to_He4_Cr48_reaclib)
              - screened_rates(k_Fe52_to_p_Mn51_reaclib);
    jac.set(Fe52, Fe52, scratch);

    scratch = screened_rates(k_Fe54_to_Fe52_n_n_approx);
    jac.set(Fe52, Fe54, scratch);

    jac.set(Fe52, Fe56, 0.0);

    jac.set(Fe52, Fe58, 0.0);

    scratch = screened_rates(k_Co55_p_to_He4_Fe52_reaclib)*Y(H1)*state.rho;
    jac.set(Fe52, Co55, scratch);

    jac.set(Fe52, Co56, 0.0);

    jac.set(Fe52, Co57, 0.0);

    jac.set(Fe52, Co59, 0.0);

    scratch = screened_rates(k_Ni56_to_He4_Fe52_reaclib);
    jac.set(Fe52, Ni56, scratch);

    jac.set(Fe52, Ni58, 0.0);

    jac.set(Fe52, Ni60, 0.0);

    jac.set(Fe52, Cu59, 0.0);

    jac.set(Fe52, Zn60, 0.0);


    // Fe54 row

    scratch = 1.0*screened_rates(k_Fe52_n_n_to_Fe54_approx)*Y(Fe52)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 1.0*screened_rates(k_Fe54_n_n_to_Fe56_approx)*Y(Fe54)*Y(N)*amrex::Math::powi<2>(state.rho)
              + 0.5*Y(Fe52)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Fe52_n_n_to_Fe54_approx_dYN
              - 0.5*Y(Fe54)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Fe54_n_n_to_Fe56_approx_dYN
              - Y(Fe54)*rate_eval.drate_Fe54_to_Fe52_n_n_approx_dYN
              + Y(Fe56)*rate_eval.drate_Fe56_to_Fe54_n_n_approx_dYN;
    jac.set(Fe54, N, scratch);

    scratch = screened_rates(k_Co57_p_to_He4_Fe54_reaclib)*Y(Co57)*state.rho
              - screened_rates(k_Fe54_p_to_Co55_reaclib)*Y(Fe54)*state.rho
              - screened_rates(k_Fe54_p_to_He4_Mn51_reaclib)*Y(Fe54)*state.rho
              + screened_rates(k_Mn53_p_to_Fe54_reaclib)*Y(Mn53)*state.rho;
    jac.set(Fe54, H1, scratch);

    scratch = screened_rates(k_Cr50_He4_to_Fe54_reaclib)*Y(Cr50)*state.rho
              - screened_rates(k_Fe54_He4_to_Ni58_reaclib)*Y(Fe54)*state.rho
              - screened_rates(k_Fe54_He4_to_p_Co57_reaclib)*Y(Fe54)*state.rho
              + screened_rates(k_Mn51_He4_to_p_Fe54_reaclib)*Y(Mn51)*state.rho;
    jac.set(Fe54, He4, scratch);

    jac.set(Fe54, C12, 0.0);

    jac.set(Fe54, O16, 0.0);

    jac.set(Fe54, Ne20, 0.0);

    jac.set(Fe54, Mg23, 0.0);

    jac.set(Fe54, Mg24, 0.0);

    jac.set(Fe54, Al27, 0.0);

    jac.set(Fe54, Si27, 0.0);

    jac.set(Fe54, Si28, 0.0);

    jac.set(Fe54, Si29, 0.0);

    jac.set(Fe54, Si30, 0.0);

    jac.set(Fe54, Si31, 0.0);

    jac.set(Fe54, P29, 0.0);

    jac.set(Fe54, P31, 0.0);

    jac.set(Fe54, P33, 0.0);

    jac.set(Fe54, S32, 0.0);

    jac.set(Fe54, S33, 0.0);

    jac.set(Fe54, S34, 0.0);

    jac.set(Fe54, S35, 0.0);

    jac.set(Fe54, Cl33, 0.0);

    jac.set(Fe54, Cl35, 0.0);

    jac.set(Fe54, Cl37, 0.0);

    jac.set(Fe54, Ar36, 0.0);

    jac.set(Fe54, Ar37, 0.0);

    jac.set(Fe54, Ar38, 0.0);

    jac.set(Fe54, K37, 0.0);

    jac.set(Fe54, K39, 0.0);

    jac.set(Fe54, Ca40, 0.0);

    jac.set(Fe54, Ca42, 0.0);

    jac.set(Fe54, Sc41, 0.0);

    jac.set(Fe54, Sc43, 0.0);

    jac.set(Fe54, Sc45, 0.0);

    jac.set(Fe54, Ti44, 0.0);

    jac.set(Fe54, Ti46, 0.0);

    jac.set(Fe54, Ti48, 0.0);

    jac.set(Fe54, V47, 0.0);

    jac.set(Fe54, V49, 0.0);

    jac.set(Fe54, V51, 0.0);

    jac.set(Fe54, Cr48, 0.0);

    scratch = screened_rates(k_Cr50_He4_to_Fe54_reaclib)*Y(He4)*state.rho;
    jac.set(Fe54, Cr50, scratch);

    jac.set(Fe54, Cr52, 0.0);

    jac.set(Fe54, Cr54, 0.0);

    scratch = screened_rates(k_Mn51_He4_to_p_Fe54_reaclib)*Y(He4)*state.rho;
    jac.set(Fe54, Mn51, scratch);

    scratch = screened_rates(k_Mn53_p_to_Fe54_reaclib)*Y(H1)*state.rho;
    jac.set(Fe54, Mn53, scratch);

    jac.set(Fe54, Mn55, 0.0);

    scratch = 0.5*screened_rates(k_Fe52_n_n_to_Fe54_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(Fe54, Fe52, scratch);

    scratch = -screened_rates(k_Fe54_He4_to_Ni58_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Fe54_He4_to_p_Co57_reaclib)*Y(He4)*state.rho
              - 0.5*screened_rates(k_Fe54_n_n_to_Fe56_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_Fe54_p_to_Co55_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Fe54_p_to_He4_Mn51_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Fe54_to_Fe52_n_n_approx)
              - screened_rates(k_Fe54_to_He4_Cr50_reaclib)
              - screened_rates(k_Fe54_to_p_Mn53_reaclib);
    jac.set(Fe54, Fe54, scratch);

    scratch = screened_rates(k_Fe56_to_Fe54_n_n_approx);
    jac.set(Fe54, Fe56, scratch);

    jac.set(Fe54, Fe58, 0.0);

    scratch = screened_rates(k_Co55_to_p_Fe54_reaclib);
    jac.set(Fe54, Co55, scratch);

    jac.set(Fe54, Co56, 0.0);

    scratch = screened_rates(k_Co57_p_to_He4_Fe54_reaclib)*Y(H1)*state.rho;
    jac.set(Fe54, Co57, scratch);

    jac.set(Fe54, Co59, 0.0);

    jac.set(Fe54, Ni56, 0.0);

    scratch = screened_rates(k_Ni58_to_He4_Fe54_reaclib);
    jac.set(Fe54, Ni58, scratch);

    jac.set(Fe54, Ni60, 0.0);

    jac.set(Fe54, Cu59, 0.0);

    jac.set(Fe54, Zn60, 0.0);


    // Fe56 row

    scratch = screened_rates(k_Co56_n_to_p_Fe56_reaclib)*Y(Co56)*state.rho
              + 1.0*screened_rates(k_Fe54_n_n_to_Fe56_approx)*Y(Fe54)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 1.0*screened_rates(k_Fe56_n_n_to_Fe58_approx)*Y(Fe56)*Y(N)*amrex::Math::powi<2>(state.rho)
              + 0.5*Y(Fe54)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Fe54_n_n_to_Fe56_approx_dYN
              - 0.5*Y(Fe56)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Fe56_n_n_to_Fe58_approx_dYN
              - Y(Fe56)*rate_eval.drate_Fe56_to_Fe54_n_n_approx_dYN
              + Y(Fe58)*rate_eval.drate_Fe58_to_Fe56_n_n_approx_dYN;
    jac.set(Fe56, N, scratch);

    scratch = screened_rates(k_Co59_p_to_He4_Fe56_reaclib)*Y(Co59)*state.rho
              - screened_rates(k_Fe56_p_to_Co57_reaclib)*Y(Fe56)*state.rho
              - screened_rates(k_Fe56_p_to_He4_Mn53_reaclib)*Y(Fe56)*state.rho
              - screened_rates(k_Fe56_p_to_n_Co56_reaclib)*Y(Fe56)*state.rho
              + screened_rates(k_Mn55_p_to_Fe56_reaclib)*Y(Mn55)*state.rho;
    jac.set(Fe56, H1, scratch);

    scratch = screened_rates(k_Cr52_He4_to_Fe56_reaclib)*Y(Cr52)*state.rho
              - screened_rates(k_Fe56_He4_to_Ni60_reaclib)*Y(Fe56)*state.rho
              - screened_rates(k_Fe56_He4_to_p_Co59_reaclib)*Y(Fe56)*state.rho
              + screened_rates(k_Mn53_He4_to_p_Fe56_reaclib)*Y(Mn53)*state.rho;
    jac.set(Fe56, He4, scratch);

    jac.set(Fe56, C12, 0.0);

    jac.set(Fe56, O16, 0.0);

    jac.set(Fe56, Ne20, 0.0);

    jac.set(Fe56, Mg23, 0.0);

    jac.set(Fe56, Mg24, 0.0);

    jac.set(Fe56, Al27, 0.0);

    jac.set(Fe56, Si27, 0.0);

    jac.set(Fe56, Si28, 0.0);

    jac.set(Fe56, Si29, 0.0);

    jac.set(Fe56, Si30, 0.0);

    jac.set(Fe56, Si31, 0.0);

    jac.set(Fe56, P29, 0.0);

    jac.set(Fe56, P31, 0.0);

    jac.set(Fe56, P33, 0.0);

    jac.set(Fe56, S32, 0.0);

    jac.set(Fe56, S33, 0.0);

    jac.set(Fe56, S34, 0.0);

    jac.set(Fe56, S35, 0.0);

    jac.set(Fe56, Cl33, 0.0);

    jac.set(Fe56, Cl35, 0.0);

    jac.set(Fe56, Cl37, 0.0);

    jac.set(Fe56, Ar36, 0.0);

    jac.set(Fe56, Ar37, 0.0);

    jac.set(Fe56, Ar38, 0.0);

    jac.set(Fe56, K37, 0.0);

    jac.set(Fe56, K39, 0.0);

    jac.set(Fe56, Ca40, 0.0);

    jac.set(Fe56, Ca42, 0.0);

    jac.set(Fe56, Sc41, 0.0);

    jac.set(Fe56, Sc43, 0.0);

    jac.set(Fe56, Sc45, 0.0);

    jac.set(Fe56, Ti44, 0.0);

    jac.set(Fe56, Ti46, 0.0);

    jac.set(Fe56, Ti48, 0.0);

    jac.set(Fe56, V47, 0.0);

    jac.set(Fe56, V49, 0.0);

    jac.set(Fe56, V51, 0.0);

    jac.set(Fe56, Cr48, 0.0);

    jac.set(Fe56, Cr50, 0.0);

    scratch = screened_rates(k_Cr52_He4_to_Fe56_reaclib)*Y(He4)*state.rho;
    jac.set(Fe56, Cr52, scratch);

    jac.set(Fe56, Cr54, 0.0);

    jac.set(Fe56, Mn51, 0.0);

    scratch = screened_rates(k_Mn53_He4_to_p_Fe56_reaclib)*Y(He4)*state.rho;
    jac.set(Fe56, Mn53, scratch);

    scratch = screened_rates(k_Mn55_p_to_Fe56_reaclib)*Y(H1)*state.rho;
    jac.set(Fe56, Mn55, scratch);

    jac.set(Fe56, Fe52, 0.0);

    scratch = 0.5*screened_rates(k_Fe54_n_n_to_Fe56_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(Fe56, Fe54, scratch);

    scratch = -screened_rates(k_Fe56_He4_to_Ni60_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Fe56_He4_to_p_Co59_reaclib)*Y(He4)*state.rho
              - 0.5*screened_rates(k_Fe56_n_n_to_Fe58_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_Fe56_p_to_Co57_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Fe56_p_to_He4_Mn53_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Fe56_p_to_n_Co56_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Fe56_to_Co56_beta_neg_weaktab)
              - screened_rates(k_Fe56_to_Fe54_n_n_approx)
              - screened_rates(k_Fe56_to_He4_Cr52_reaclib)
              - screened_rates(k_Fe56_to_p_Mn55_reaclib);
    jac.set(Fe56, Fe56, scratch);

    scratch = screened_rates(k_Fe58_to_Fe56_n_n_approx);
    jac.set(Fe56, Fe58, scratch);

    jac.set(Fe56, Co55, 0.0);

    scratch = screened_rates(k_Co56_n_to_p_Fe56_reaclib)*Y(N)*state.rho
              + screened_rates(k_Co56_to_Fe56_electron_capture_weaktab);
    jac.set(Fe56, Co56, scratch);

    scratch = screened_rates(k_Co57_to_p_Fe56_reaclib);
    jac.set(Fe56, Co57, scratch);

    scratch = screened_rates(k_Co59_p_to_He4_Fe56_reaclib)*Y(H1)*state.rho;
    jac.set(Fe56, Co59, scratch);

    jac.set(Fe56, Ni56, 0.0);

    jac.set(Fe56, Ni58, 0.0);

    scratch = screened_rates(k_Ni60_to_He4_Fe56_reaclib);
    jac.set(Fe56, Ni60, scratch);

    jac.set(Fe56, Cu59, 0.0);

    jac.set(Fe56, Zn60, 0.0);


    // Fe58 row

    scratch = 1.0*screened_rates(k_Fe56_n_n_to_Fe58_approx)*Y(Fe56)*Y(N)*amrex::Math::powi<2>(state.rho)
              + 0.5*Y(Fe56)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Fe56_n_n_to_Fe58_approx_dYN
              - Y(Fe58)*rate_eval.drate_Fe58_to_Fe56_n_n_approx_dYN;
    jac.set(Fe58, N, scratch);

    scratch = -screened_rates(k_Fe58_p_to_Co59_reaclib)*Y(Fe58)*state.rho
              - screened_rates(k_Fe58_p_to_He4_Mn55_reaclib)*Y(Fe58)*state.rho;
    jac.set(Fe58, H1, scratch);

    scratch = screened_rates(k_Cr54_He4_to_Fe58_reaclib)*Y(Cr54)*state.rho
              + screened_rates(k_Mn55_He4_to_p_Fe58_reaclib)*Y(Mn55)*state.rho;
    jac.set(Fe58, He4, scratch);

    jac.set(Fe58, C12, 0.0);

    jac.set(Fe58, O16, 0.0);

    jac.set(Fe58, Ne20, 0.0);

    jac.set(Fe58, Mg23, 0.0);

    jac.set(Fe58, Mg24, 0.0);

    jac.set(Fe58, Al27, 0.0);

    jac.set(Fe58, Si27, 0.0);

    jac.set(Fe58, Si28, 0.0);

    jac.set(Fe58, Si29, 0.0);

    jac.set(Fe58, Si30, 0.0);

    jac.set(Fe58, Si31, 0.0);

    jac.set(Fe58, P29, 0.0);

    jac.set(Fe58, P31, 0.0);

    jac.set(Fe58, P33, 0.0);

    jac.set(Fe58, S32, 0.0);

    jac.set(Fe58, S33, 0.0);

    jac.set(Fe58, S34, 0.0);

    jac.set(Fe58, S35, 0.0);

    jac.set(Fe58, Cl33, 0.0);

    jac.set(Fe58, Cl35, 0.0);

    jac.set(Fe58, Cl37, 0.0);

    jac.set(Fe58, Ar36, 0.0);

    jac.set(Fe58, Ar37, 0.0);

    jac.set(Fe58, Ar38, 0.0);

    jac.set(Fe58, K37, 0.0);

    jac.set(Fe58, K39, 0.0);

    jac.set(Fe58, Ca40, 0.0);

    jac.set(Fe58, Ca42, 0.0);

    jac.set(Fe58, Sc41, 0.0);

    jac.set(Fe58, Sc43, 0.0);

    jac.set(Fe58, Sc45, 0.0);

    jac.set(Fe58, Ti44, 0.0);

    jac.set(Fe58, Ti46, 0.0);

    jac.set(Fe58, Ti48, 0.0);

    jac.set(Fe58, V47, 0.0);

    jac.set(Fe58, V49, 0.0);

    jac.set(Fe58, V51, 0.0);

    jac.set(Fe58, Cr48, 0.0);

    jac.set(Fe58, Cr50, 0.0);

    jac.set(Fe58, Cr52, 0.0);

    scratch = screened_rates(k_Cr54_He4_to_Fe58_reaclib)*Y(He4)*state.rho;
    jac.set(Fe58, Cr54, scratch);

    jac.set(Fe58, Mn51, 0.0);

    jac.set(Fe58, Mn53, 0.0);

    scratch = screened_rates(k_Mn55_He4_to_p_Fe58_reaclib)*Y(He4)*state.rho;
    jac.set(Fe58, Mn55, scratch);

    jac.set(Fe58, Fe52, 0.0);

    jac.set(Fe58, Fe54, 0.0);

    scratch = 0.5*screened_rates(k_Fe56_n_n_to_Fe58_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(Fe58, Fe56, scratch);

    scratch = -screened_rates(k_Fe58_p_to_Co59_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Fe58_p_to_He4_Mn55_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Fe58_to_Fe56_n_n_approx)
              - screened_rates(k_Fe58_to_He4_Cr54_reaclib);
    jac.set(Fe58, Fe58, scratch);

    jac.set(Fe58, Co55, 0.0);

    jac.set(Fe58, Co56, 0.0);

    jac.set(Fe58, Co57, 0.0);

    scratch = screened_rates(k_Co59_to_p_Fe58_reaclib);
    jac.set(Fe58, Co59, scratch);

    jac.set(Fe58, Ni56, 0.0);

    jac.set(Fe58, Ni58, 0.0);

    jac.set(Fe58, Ni60, 0.0);

    jac.set(Fe58, Cu59, 0.0);

    jac.set(Fe58, Zn60, 0.0);


    // Co55 row

    scratch = -screened_rates(k_Co55_n_to_Co56_reaclib)*Y(Co55)*state.rho;
    jac.set(Co55, N, scratch);

    scratch = -screened_rates(k_Co55_p_to_He4_Fe52_reaclib)*Y(Co55)*state.rho
              - screened_rates(k_Co55_p_to_Ni56_reaclib)*Y(Co55)*state.rho
              + screened_rates(k_Fe54_p_to_Co55_reaclib)*Y(Fe54)*state.rho
              + screened_rates(k_Ni58_p_to_He4_Co55_reaclib)*Y(Ni58)*state.rho;
    jac.set(Co55, H1, scratch);

    scratch = -screened_rates(k_Co55_He4_to_Cu59_reaclib)*Y(Co55)*state.rho
              - screened_rates(k_Co55_He4_to_p_Ni58_reaclib)*Y(Co55)*state.rho
              + screened_rates(k_Fe52_He4_to_p_Co55_reaclib)*Y(Fe52)*state.rho
              + screened_rates(k_Mn51_He4_to_Co55_reaclib)*Y(Mn51)*state.rho;
    jac.set(Co55, He4, scratch);

    jac.set(Co55, C12, 0.0);

    jac.set(Co55, O16, 0.0);

    jac.set(Co55, Ne20, 0.0);

    jac.set(Co55, Mg23, 0.0);

    jac.set(Co55, Mg24, 0.0);

    jac.set(Co55, Al27, 0.0);

    jac.set(Co55, Si27, 0.0);

    jac.set(Co55, Si28, 0.0);

    jac.set(Co55, Si29, 0.0);

    jac.set(Co55, Si30, 0.0);

    jac.set(Co55, Si31, 0.0);

    jac.set(Co55, P29, 0.0);

    jac.set(Co55, P31, 0.0);

    jac.set(Co55, P33, 0.0);

    jac.set(Co55, S32, 0.0);

    jac.set(Co55, S33, 0.0);

    jac.set(Co55, S34, 0.0);

    jac.set(Co55, S35, 0.0);

    jac.set(Co55, Cl33, 0.0);

    jac.set(Co55, Cl35, 0.0);

    jac.set(Co55, Cl37, 0.0);

    jac.set(Co55, Ar36, 0.0);

    jac.set(Co55, Ar37, 0.0);

    jac.set(Co55, Ar38, 0.0);

    jac.set(Co55, K37, 0.0);

    jac.set(Co55, K39, 0.0);

    jac.set(Co55, Ca40, 0.0);

    jac.set(Co55, Ca42, 0.0);

    jac.set(Co55, Sc41, 0.0);

    jac.set(Co55, Sc43, 0.0);

    jac.set(Co55, Sc45, 0.0);

    jac.set(Co55, Ti44, 0.0);

    jac.set(Co55, Ti46, 0.0);

    jac.set(Co55, Ti48, 0.0);

    jac.set(Co55, V47, 0.0);

    jac.set(Co55, V49, 0.0);

    jac.set(Co55, V51, 0.0);

    jac.set(Co55, Cr48, 0.0);

    jac.set(Co55, Cr50, 0.0);

    jac.set(Co55, Cr52, 0.0);

    jac.set(Co55, Cr54, 0.0);

    scratch = screened_rates(k_Mn51_He4_to_Co55_reaclib)*Y(He4)*state.rho;
    jac.set(Co55, Mn51, scratch);

    jac.set(Co55, Mn53, 0.0);

    jac.set(Co55, Mn55, 0.0);

    scratch = screened_rates(k_Fe52_He4_to_p_Co55_reaclib)*Y(He4)*state.rho;
    jac.set(Co55, Fe52, scratch);

    scratch = screened_rates(k_Fe54_p_to_Co55_reaclib)*Y(H1)*state.rho;
    jac.set(Co55, Fe54, scratch);

    jac.set(Co55, Fe56, 0.0);

    jac.set(Co55, Fe58, 0.0);

    scratch = -screened_rates(k_Co55_He4_to_Cu59_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Co55_He4_to_p_Ni58_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Co55_n_to_Co56_reaclib)*Y(N)*state.rho
              - screened_rates(k_Co55_p_to_He4_Fe52_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Co55_p_to_Ni56_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Co55_to_He4_Mn51_reaclib)
              - screened_rates(k_Co55_to_p_Fe54_reaclib);
    jac.set(Co55, Co55, scratch);

    scratch = screened_rates(k_Co56_to_n_Co55_reaclib);
    jac.set(Co55, Co56, scratch);

    jac.set(Co55, Co57, 0.0);

    jac.set(Co55, Co59, 0.0);

    scratch = screened_rates(k_Ni56_to_p_Co55_reaclib);
    jac.set(Co55, Ni56, scratch);

    scratch = screened_rates(k_Ni58_p_to_He4_Co55_reaclib)*Y(H1)*state.rho;
    jac.set(Co55, Ni58, scratch);

    jac.set(Co55, Ni60, 0.0);

    scratch = screened_rates(k_Cu59_to_He4_Co55_reaclib);
    jac.set(Co55, Cu59, scratch);

    jac.set(Co55, Zn60, 0.0);


    // Co56 row

    scratch = screened_rates(k_Co55_n_to_Co56_reaclib)*Y(Co55)*state.rho
              - screened_rates(k_Co56_n_to_Co57_reaclib)*Y(Co56)*state.rho
              - screened_rates(k_Co56_n_to_He4_Mn53_reaclib)*Y(Co56)*state.rho
              - screened_rates(k_Co56_n_to_p_Fe56_reaclib)*Y(Co56)*state.rho
              + screened_rates(k_Cu59_n_to_He4_Co56_reaclib)*Y(Cu59)*state.rho
              + screened_rates(k_Ni56_n_to_p_Co56_reaclib)*Y(Ni56)*state.rho;
    jac.set(Co56, N, scratch);

    scratch = -screened_rates(k_Co56_p_to_n_Ni56_reaclib)*Y(Co56)*state.rho
              + screened_rates(k_Fe56_p_to_n_Co56_reaclib)*Y(Fe56)*state.rho;
    jac.set(Co56, H1, scratch);

    scratch = -screened_rates(k_Co56_He4_to_n_Cu59_reaclib)*Y(Co56)*state.rho
              + screened_rates(k_Mn53_He4_to_n_Co56_reaclib)*Y(Mn53)*state.rho;
    jac.set(Co56, He4, scratch);

    jac.set(Co56, C12, 0.0);

    jac.set(Co56, O16, 0.0);

    jac.set(Co56, Ne20, 0.0);

    jac.set(Co56, Mg23, 0.0);

    jac.set(Co56, Mg24, 0.0);

    jac.set(Co56, Al27, 0.0);

    jac.set(Co56, Si27, 0.0);

    jac.set(Co56, Si28, 0.0);

    jac.set(Co56, Si29, 0.0);

    jac.set(Co56, Si30, 0.0);

    jac.set(Co56, Si31, 0.0);

    jac.set(Co56, P29, 0.0);

    jac.set(Co56, P31, 0.0);

    jac.set(Co56, P33, 0.0);

    jac.set(Co56, S32, 0.0);

    jac.set(Co56, S33, 0.0);

    jac.set(Co56, S34, 0.0);

    jac.set(Co56, S35, 0.0);

    jac.set(Co56, Cl33, 0.0);

    jac.set(Co56, Cl35, 0.0);

    jac.set(Co56, Cl37, 0.0);

    jac.set(Co56, Ar36, 0.0);

    jac.set(Co56, Ar37, 0.0);

    jac.set(Co56, Ar38, 0.0);

    jac.set(Co56, K37, 0.0);

    jac.set(Co56, K39, 0.0);

    jac.set(Co56, Ca40, 0.0);

    jac.set(Co56, Ca42, 0.0);

    jac.set(Co56, Sc41, 0.0);

    jac.set(Co56, Sc43, 0.0);

    jac.set(Co56, Sc45, 0.0);

    jac.set(Co56, Ti44, 0.0);

    jac.set(Co56, Ti46, 0.0);

    jac.set(Co56, Ti48, 0.0);

    jac.set(Co56, V47, 0.0);

    jac.set(Co56, V49, 0.0);

    jac.set(Co56, V51, 0.0);

    jac.set(Co56, Cr48, 0.0);

    jac.set(Co56, Cr50, 0.0);

    jac.set(Co56, Cr52, 0.0);

    jac.set(Co56, Cr54, 0.0);

    jac.set(Co56, Mn51, 0.0);

    scratch = screened_rates(k_Mn53_He4_to_n_Co56_reaclib)*Y(He4)*state.rho;
    jac.set(Co56, Mn53, scratch);

    jac.set(Co56, Mn55, 0.0);

    jac.set(Co56, Fe52, 0.0);

    jac.set(Co56, Fe54, 0.0);

    scratch = screened_rates(k_Fe56_p_to_n_Co56_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Fe56_to_Co56_beta_neg_weaktab);
    jac.set(Co56, Fe56, scratch);

    jac.set(Co56, Fe58, 0.0);

    scratch = screened_rates(k_Co55_n_to_Co56_reaclib)*Y(N)*state.rho;
    jac.set(Co56, Co55, scratch);

    scratch = -screened_rates(k_Co56_He4_to_n_Cu59_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Co56_n_to_Co57_reaclib)*Y(N)*state.rho
              - screened_rates(k_Co56_n_to_He4_Mn53_reaclib)*Y(N)*state.rho
              - screened_rates(k_Co56_n_to_p_Fe56_reaclib)*Y(N)*state.rho
              - screened_rates(k_Co56_p_to_n_Ni56_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Co56_to_Fe56_electron_capture_weaktab)
              - screened_rates(k_Co56_to_Ni56_beta_neg_weaktab)
              - screened_rates(k_Co56_to_n_Co55_reaclib);
    jac.set(Co56, Co56, scratch);

    scratch = screened_rates(k_Co57_to_n_Co56_reaclib);
    jac.set(Co56, Co57, scratch);

    jac.set(Co56, Co59, 0.0);

    scratch = screened_rates(k_Ni56_n_to_p_Co56_reaclib)*Y(N)*state.rho
              + screened_rates(k_Ni56_to_Co56_electron_capture_weaktab);
    jac.set(Co56, Ni56, scratch);

    jac.set(Co56, Ni58, 0.0);

    jac.set(Co56, Ni60, 0.0);

    scratch = screened_rates(k_Cu59_n_to_He4_Co56_reaclib)*Y(N)*state.rho;
    jac.set(Co56, Cu59, scratch);

    jac.set(Co56, Zn60, 0.0);


    // Co57 row

    scratch = screened_rates(k_Co56_n_to_Co57_reaclib)*Y(Co56)*state.rho
              - 1.0*screened_rates(k_Co57_n_n_to_Co59_approx)*Y(Co57)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 0.5*Y(Co57)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Co57_n_n_to_Co59_approx_dYN
              + Y(Co59)*rate_eval.drate_Co59_to_Co57_n_n_approx_dYN;
    jac.set(Co57, N, scratch);

    scratch = -screened_rates(k_Co57_p_to_He4_Fe54_reaclib)*Y(Co57)*state.rho
              - screened_rates(k_Co57_p_to_Ni58_reaclib)*Y(Co57)*state.rho
              + screened_rates(k_Fe56_p_to_Co57_reaclib)*Y(Fe56)*state.rho
              + screened_rates(k_Ni60_p_to_He4_Co57_reaclib)*Y(Ni60)*state.rho;
    jac.set(Co57, H1, scratch);

    scratch = -screened_rates(k_Co57_He4_to_p_Ni60_reaclib)*Y(Co57)*state.rho
              + screened_rates(k_Fe54_He4_to_p_Co57_reaclib)*Y(Fe54)*state.rho
              + screened_rates(k_Mn53_He4_to_Co57_reaclib)*Y(Mn53)*state.rho;
    jac.set(Co57, He4, scratch);

    jac.set(Co57, C12, 0.0);

    jac.set(Co57, O16, 0.0);

    jac.set(Co57, Ne20, 0.0);

    jac.set(Co57, Mg23, 0.0);

    jac.set(Co57, Mg24, 0.0);

    jac.set(Co57, Al27, 0.0);

    jac.set(Co57, Si27, 0.0);

    jac.set(Co57, Si28, 0.0);

    jac.set(Co57, Si29, 0.0);

    jac.set(Co57, Si30, 0.0);

    jac.set(Co57, Si31, 0.0);

    jac.set(Co57, P29, 0.0);

    jac.set(Co57, P31, 0.0);

    jac.set(Co57, P33, 0.0);

    jac.set(Co57, S32, 0.0);

    jac.set(Co57, S33, 0.0);

    jac.set(Co57, S34, 0.0);

    jac.set(Co57, S35, 0.0);

    jac.set(Co57, Cl33, 0.0);

    jac.set(Co57, Cl35, 0.0);

    jac.set(Co57, Cl37, 0.0);

    jac.set(Co57, Ar36, 0.0);

    jac.set(Co57, Ar37, 0.0);

    jac.set(Co57, Ar38, 0.0);

    jac.set(Co57, K37, 0.0);

    jac.set(Co57, K39, 0.0);

    jac.set(Co57, Ca40, 0.0);

    jac.set(Co57, Ca42, 0.0);

    jac.set(Co57, Sc41, 0.0);

    jac.set(Co57, Sc43, 0.0);

    jac.set(Co57, Sc45, 0.0);

    jac.set(Co57, Ti44, 0.0);

    jac.set(Co57, Ti46, 0.0);

    jac.set(Co57, Ti48, 0.0);

    jac.set(Co57, V47, 0.0);

    jac.set(Co57, V49, 0.0);

    jac.set(Co57, V51, 0.0);

    jac.set(Co57, Cr48, 0.0);

    jac.set(Co57, Cr50, 0.0);

    jac.set(Co57, Cr52, 0.0);

    jac.set(Co57, Cr54, 0.0);

    jac.set(Co57, Mn51, 0.0);

    scratch = screened_rates(k_Mn53_He4_to_Co57_reaclib)*Y(He4)*state.rho;
    jac.set(Co57, Mn53, scratch);

    jac.set(Co57, Mn55, 0.0);

    jac.set(Co57, Fe52, 0.0);

    scratch = screened_rates(k_Fe54_He4_to_p_Co57_reaclib)*Y(He4)*state.rho;
    jac.set(Co57, Fe54, scratch);

    scratch = screened_rates(k_Fe56_p_to_Co57_reaclib)*Y(H1)*state.rho;
    jac.set(Co57, Fe56, scratch);

    jac.set(Co57, Fe58, 0.0);

    jac.set(Co57, Co55, 0.0);

    scratch = screened_rates(k_Co56_n_to_Co57_reaclib)*Y(N)*state.rho;
    jac.set(Co57, Co56, scratch);

    scratch = -screened_rates(k_Co57_He4_to_p_Ni60_reaclib)*Y(He4)*state.rho
              - 0.5*screened_rates(k_Co57_n_n_to_Co59_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_Co57_p_to_He4_Fe54_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Co57_p_to_Ni58_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Co57_to_He4_Mn53_reaclib)
              - screened_rates(k_Co57_to_n_Co56_reaclib)
              - screened_rates(k_Co57_to_p_Fe56_reaclib);
    jac.set(Co57, Co57, scratch);

    scratch = screened_rates(k_Co59_to_Co57_n_n_approx);
    jac.set(Co57, Co59, scratch);

    jac.set(Co57, Ni56, 0.0);

    scratch = screened_rates(k_Ni58_to_p_Co57_reaclib);
    jac.set(Co57, Ni58, scratch);

    scratch = screened_rates(k_Ni60_p_to_He4_Co57_reaclib)*Y(H1)*state.rho;
    jac.set(Co57, Ni60, scratch);

    jac.set(Co57, Cu59, 0.0);

    jac.set(Co57, Zn60, 0.0);


    // Co59 row

    scratch = 1.0*screened_rates(k_Co57_n_n_to_Co59_approx)*Y(Co57)*Y(N)*amrex::Math::powi<2>(state.rho)
              + 0.5*Y(Co57)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Co57_n_n_to_Co59_approx_dYN
              - Y(Co59)*rate_eval.drate_Co59_to_Co57_n_n_approx_dYN;
    jac.set(Co59, N, scratch);

    scratch = -screened_rates(k_Co59_p_to_He4_Fe56_reaclib)*Y(Co59)*state.rho
              - screened_rates(k_Co59_p_to_Ni60_reaclib)*Y(Co59)*state.rho
              + screened_rates(k_Fe58_p_to_Co59_reaclib)*Y(Fe58)*state.rho;
    jac.set(Co59, H1, scratch);

    scratch = screened_rates(k_Fe56_He4_to_p_Co59_reaclib)*Y(Fe56)*state.rho
              + screened_rates(k_Mn55_He4_to_Co59_reaclib)*Y(Mn55)*state.rho;
    jac.set(Co59, He4, scratch);

    jac.set(Co59, C12, 0.0);

    jac.set(Co59, O16, 0.0);

    jac.set(Co59, Ne20, 0.0);

    jac.set(Co59, Mg23, 0.0);

    jac.set(Co59, Mg24, 0.0);

    jac.set(Co59, Al27, 0.0);

    jac.set(Co59, Si27, 0.0);

    jac.set(Co59, Si28, 0.0);

    jac.set(Co59, Si29, 0.0);

    jac.set(Co59, Si30, 0.0);

    jac.set(Co59, Si31, 0.0);

    jac.set(Co59, P29, 0.0);

    jac.set(Co59, P31, 0.0);

    jac.set(Co59, P33, 0.0);

    jac.set(Co59, S32, 0.0);

    jac.set(Co59, S33, 0.0);

    jac.set(Co59, S34, 0.0);

    jac.set(Co59, S35, 0.0);

    jac.set(Co59, Cl33, 0.0);

    jac.set(Co59, Cl35, 0.0);

    jac.set(Co59, Cl37, 0.0);

    jac.set(Co59, Ar36, 0.0);

    jac.set(Co59, Ar37, 0.0);

    jac.set(Co59, Ar38, 0.0);

    jac.set(Co59, K37, 0.0);

    jac.set(Co59, K39, 0.0);

    jac.set(Co59, Ca40, 0.0);

    jac.set(Co59, Ca42, 0.0);

    jac.set(Co59, Sc41, 0.0);

    jac.set(Co59, Sc43, 0.0);

    jac.set(Co59, Sc45, 0.0);

    jac.set(Co59, Ti44, 0.0);

    jac.set(Co59, Ti46, 0.0);

    jac.set(Co59, Ti48, 0.0);

    jac.set(Co59, V47, 0.0);

    jac.set(Co59, V49, 0.0);

    jac.set(Co59, V51, 0.0);

    jac.set(Co59, Cr48, 0.0);

    jac.set(Co59, Cr50, 0.0);

    jac.set(Co59, Cr52, 0.0);

    jac.set(Co59, Cr54, 0.0);

    jac.set(Co59, Mn51, 0.0);

    jac.set(Co59, Mn53, 0.0);

    scratch = screened_rates(k_Mn55_He4_to_Co59_reaclib)*Y(He4)*state.rho;
    jac.set(Co59, Mn55, scratch);

    jac.set(Co59, Fe52, 0.0);

    jac.set(Co59, Fe54, 0.0);

    scratch = screened_rates(k_Fe56_He4_to_p_Co59_reaclib)*Y(He4)*state.rho;
    jac.set(Co59, Fe56, scratch);

    scratch = screened_rates(k_Fe58_p_to_Co59_reaclib)*Y(H1)*state.rho;
    jac.set(Co59, Fe58, scratch);

    jac.set(Co59, Co55, 0.0);

    jac.set(Co59, Co56, 0.0);

    scratch = 0.5*screened_rates(k_Co57_n_n_to_Co59_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(Co59, Co57, scratch);

    scratch = -screened_rates(k_Co59_p_to_He4_Fe56_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Co59_p_to_Ni60_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Co59_to_Co57_n_n_approx)
              - screened_rates(k_Co59_to_He4_Mn55_reaclib)
              - screened_rates(k_Co59_to_p_Fe58_reaclib);
    jac.set(Co59, Co59, scratch);

    jac.set(Co59, Ni56, 0.0);

    jac.set(Co59, Ni58, 0.0);

    scratch = screened_rates(k_Ni60_to_p_Co59_reaclib);
    jac.set(Co59, Ni60, scratch);

    jac.set(Co59, Cu59, 0.0);

    jac.set(Co59, Zn60, 0.0);


    // Ni56 row

    scratch = -1.0*screened_rates(k_Ni56_n_n_to_Ni58_approx)*Y(Ni56)*Y(N)*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_Ni56_n_to_p_Co56_reaclib)*Y(Ni56)*state.rho
              - 0.5*Y(Ni56)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Ni56_n_n_to_Ni58_approx_dYN
              + Y(Ni58)*rate_eval.drate_Ni58_to_Ni56_n_n_approx_dYN;
    jac.set(Ni56, N, scratch);

    scratch = screened_rates(k_Co55_p_to_Ni56_reaclib)*Y(Co55)*state.rho
              + screened_rates(k_Co56_p_to_n_Ni56_reaclib)*Y(Co56)*state.rho
              + screened_rates(k_Cu59_p_to_He4_Ni56_reaclib)*Y(Cu59)*state.rho;
    jac.set(Ni56, H1, scratch);

    scratch = screened_rates(k_Fe52_He4_to_Ni56_reaclib)*Y(Fe52)*state.rho
              - screened_rates(k_Ni56_He4_to_Zn60_reaclib)*Y(Ni56)*state.rho
              - screened_rates(k_Ni56_He4_to_p_Cu59_reaclib)*Y(Ni56)*state.rho;
    jac.set(Ni56, He4, scratch);

    jac.set(Ni56, C12, 0.0);

    jac.set(Ni56, O16, 0.0);

    jac.set(Ni56, Ne20, 0.0);

    jac.set(Ni56, Mg23, 0.0);

    jac.set(Ni56, Mg24, 0.0);

    jac.set(Ni56, Al27, 0.0);

    jac.set(Ni56, Si27, 0.0);

    jac.set(Ni56, Si28, 0.0);

    jac.set(Ni56, Si29, 0.0);

    jac.set(Ni56, Si30, 0.0);

    jac.set(Ni56, Si31, 0.0);

    jac.set(Ni56, P29, 0.0);

    jac.set(Ni56, P31, 0.0);

    jac.set(Ni56, P33, 0.0);

    jac.set(Ni56, S32, 0.0);

    jac.set(Ni56, S33, 0.0);

    jac.set(Ni56, S34, 0.0);

    jac.set(Ni56, S35, 0.0);

    jac.set(Ni56, Cl33, 0.0);

    jac.set(Ni56, Cl35, 0.0);

    jac.set(Ni56, Cl37, 0.0);

    jac.set(Ni56, Ar36, 0.0);

    jac.set(Ni56, Ar37, 0.0);

    jac.set(Ni56, Ar38, 0.0);

    jac.set(Ni56, K37, 0.0);

    jac.set(Ni56, K39, 0.0);

    jac.set(Ni56, Ca40, 0.0);

    jac.set(Ni56, Ca42, 0.0);

    jac.set(Ni56, Sc41, 0.0);

    jac.set(Ni56, Sc43, 0.0);

    jac.set(Ni56, Sc45, 0.0);

    jac.set(Ni56, Ti44, 0.0);

    jac.set(Ni56, Ti46, 0.0);

    jac.set(Ni56, Ti48, 0.0);

    jac.set(Ni56, V47, 0.0);

    jac.set(Ni56, V49, 0.0);

    jac.set(Ni56, V51, 0.0);

    jac.set(Ni56, Cr48, 0.0);

    jac.set(Ni56, Cr50, 0.0);

    jac.set(Ni56, Cr52, 0.0);

    jac.set(Ni56, Cr54, 0.0);

    jac.set(Ni56, Mn51, 0.0);

    jac.set(Ni56, Mn53, 0.0);

    jac.set(Ni56, Mn55, 0.0);

    scratch = screened_rates(k_Fe52_He4_to_Ni56_reaclib)*Y(He4)*state.rho;
    jac.set(Ni56, Fe52, scratch);

    jac.set(Ni56, Fe54, 0.0);

    jac.set(Ni56, Fe56, 0.0);

    jac.set(Ni56, Fe58, 0.0);

    scratch = screened_rates(k_Co55_p_to_Ni56_reaclib)*Y(H1)*state.rho;
    jac.set(Ni56, Co55, scratch);

    scratch = screened_rates(k_Co56_p_to_n_Ni56_reaclib)*Y(H1)*state.rho
              + screened_rates(k_Co56_to_Ni56_beta_neg_weaktab);
    jac.set(Ni56, Co56, scratch);

    jac.set(Ni56, Co57, 0.0);

    jac.set(Ni56, Co59, 0.0);

    scratch = -screened_rates(k_Ni56_He4_to_Zn60_reaclib)*Y(He4)*state.rho
              - screened_rates(k_Ni56_He4_to_p_Cu59_reaclib)*Y(He4)*state.rho
              - 0.5*screened_rates(k_Ni56_n_n_to_Ni58_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_Ni56_n_to_p_Co56_reaclib)*Y(N)*state.rho
              - screened_rates(k_Ni56_to_Co56_electron_capture_weaktab)
              - screened_rates(k_Ni56_to_He4_Fe52_reaclib)
              - screened_rates(k_Ni56_to_p_Co55_reaclib);
    jac.set(Ni56, Ni56, scratch);

    scratch = screened_rates(k_Ni58_to_Ni56_n_n_approx);
    jac.set(Ni56, Ni58, scratch);

    jac.set(Ni56, Ni60, 0.0);

    scratch = screened_rates(k_Cu59_p_to_He4_Ni56_reaclib)*Y(H1)*state.rho;
    jac.set(Ni56, Cu59, scratch);

    scratch = screened_rates(k_Zn60_to_He4_Ni56_reaclib);
    jac.set(Ni56, Zn60, scratch);


    // Ni58 row

    scratch = 1.0*screened_rates(k_Ni56_n_n_to_Ni58_approx)*Y(Ni56)*Y(N)*amrex::Math::powi<2>(state.rho)
              - 1.0*screened_rates(k_Ni58_n_n_to_Ni60_approx)*Y(Ni58)*Y(N)*amrex::Math::powi<2>(state.rho)
              + 0.5*Y(Ni56)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Ni56_n_n_to_Ni58_approx_dYN
              - 0.5*Y(Ni58)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Ni58_n_n_to_Ni60_approx_dYN
              - Y(Ni58)*rate_eval.drate_Ni58_to_Ni56_n_n_approx_dYN
              + Y(Ni60)*rate_eval.drate_Ni60_to_Ni58_n_n_approx_dYN;
    jac.set(Ni58, N, scratch);

    scratch = screened_rates(k_Co57_p_to_Ni58_reaclib)*Y(Co57)*state.rho
              - screened_rates(k_Ni58_p_to_Cu59_reaclib)*Y(Ni58)*state.rho
              - screened_rates(k_Ni58_p_to_He4_Co55_reaclib)*Y(Ni58)*state.rho;
    jac.set(Ni58, H1, scratch);

    scratch = screened_rates(k_Co55_He4_to_p_Ni58_reaclib)*Y(Co55)*state.rho
              + screened_rates(k_Fe54_He4_to_Ni58_reaclib)*Y(Fe54)*state.rho;
    jac.set(Ni58, He4, scratch);

    jac.set(Ni58, C12, 0.0);

    jac.set(Ni58, O16, 0.0);

    jac.set(Ni58, Ne20, 0.0);

    jac.set(Ni58, Mg23, 0.0);

    jac.set(Ni58, Mg24, 0.0);

    jac.set(Ni58, Al27, 0.0);

    jac.set(Ni58, Si27, 0.0);

    jac.set(Ni58, Si28, 0.0);

    jac.set(Ni58, Si29, 0.0);

    jac.set(Ni58, Si30, 0.0);

    jac.set(Ni58, Si31, 0.0);

    jac.set(Ni58, P29, 0.0);

    jac.set(Ni58, P31, 0.0);

    jac.set(Ni58, P33, 0.0);

    jac.set(Ni58, S32, 0.0);

    jac.set(Ni58, S33, 0.0);

    jac.set(Ni58, S34, 0.0);

    jac.set(Ni58, S35, 0.0);

    jac.set(Ni58, Cl33, 0.0);

    jac.set(Ni58, Cl35, 0.0);

    jac.set(Ni58, Cl37, 0.0);

    jac.set(Ni58, Ar36, 0.0);

    jac.set(Ni58, Ar37, 0.0);

    jac.set(Ni58, Ar38, 0.0);

    jac.set(Ni58, K37, 0.0);

    jac.set(Ni58, K39, 0.0);

    jac.set(Ni58, Ca40, 0.0);

    jac.set(Ni58, Ca42, 0.0);

    jac.set(Ni58, Sc41, 0.0);

    jac.set(Ni58, Sc43, 0.0);

    jac.set(Ni58, Sc45, 0.0);

    jac.set(Ni58, Ti44, 0.0);

    jac.set(Ni58, Ti46, 0.0);

    jac.set(Ni58, Ti48, 0.0);

    jac.set(Ni58, V47, 0.0);

    jac.set(Ni58, V49, 0.0);

    jac.set(Ni58, V51, 0.0);

    jac.set(Ni58, Cr48, 0.0);

    jac.set(Ni58, Cr50, 0.0);

    jac.set(Ni58, Cr52, 0.0);

    jac.set(Ni58, Cr54, 0.0);

    jac.set(Ni58, Mn51, 0.0);

    jac.set(Ni58, Mn53, 0.0);

    jac.set(Ni58, Mn55, 0.0);

    jac.set(Ni58, Fe52, 0.0);

    scratch = screened_rates(k_Fe54_He4_to_Ni58_reaclib)*Y(He4)*state.rho;
    jac.set(Ni58, Fe54, scratch);

    jac.set(Ni58, Fe56, 0.0);

    jac.set(Ni58, Fe58, 0.0);

    scratch = screened_rates(k_Co55_He4_to_p_Ni58_reaclib)*Y(He4)*state.rho;
    jac.set(Ni58, Co55, scratch);

    jac.set(Ni58, Co56, 0.0);

    scratch = screened_rates(k_Co57_p_to_Ni58_reaclib)*Y(H1)*state.rho;
    jac.set(Ni58, Co57, scratch);

    jac.set(Ni58, Co59, 0.0);

    scratch = 0.5*screened_rates(k_Ni56_n_n_to_Ni58_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(Ni58, Ni56, scratch);

    scratch = -0.5*screened_rates(k_Ni58_n_n_to_Ni60_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)
              - screened_rates(k_Ni58_p_to_Cu59_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Ni58_p_to_He4_Co55_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Ni58_to_He4_Fe54_reaclib)
              - screened_rates(k_Ni58_to_Ni56_n_n_approx)
              - screened_rates(k_Ni58_to_p_Co57_reaclib);
    jac.set(Ni58, Ni58, scratch);

    scratch = screened_rates(k_Ni60_to_Ni58_n_n_approx);
    jac.set(Ni58, Ni60, scratch);

    scratch = screened_rates(k_Cu59_to_p_Ni58_reaclib);
    jac.set(Ni58, Cu59, scratch);

    jac.set(Ni58, Zn60, 0.0);


    // Ni60 row

    scratch = 1.0*screened_rates(k_Ni58_n_n_to_Ni60_approx)*Y(Ni58)*Y(N)*amrex::Math::powi<2>(state.rho)
              + 0.5*Y(Ni58)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho)*rate_eval.drate_Ni58_n_n_to_Ni60_approx_dYN
              - Y(Ni60)*rate_eval.drate_Ni60_to_Ni58_n_n_approx_dYN;
    jac.set(Ni60, N, scratch);

    scratch = screened_rates(k_Co59_p_to_Ni60_reaclib)*Y(Co59)*state.rho
              - screened_rates(k_Ni60_p_to_He4_Co57_reaclib)*Y(Ni60)*state.rho;
    jac.set(Ni60, H1, scratch);

    scratch = screened_rates(k_Co57_He4_to_p_Ni60_reaclib)*Y(Co57)*state.rho
              + screened_rates(k_Fe56_He4_to_Ni60_reaclib)*Y(Fe56)*state.rho;
    jac.set(Ni60, He4, scratch);

    jac.set(Ni60, C12, 0.0);

    jac.set(Ni60, O16, 0.0);

    jac.set(Ni60, Ne20, 0.0);

    jac.set(Ni60, Mg23, 0.0);

    jac.set(Ni60, Mg24, 0.0);

    jac.set(Ni60, Al27, 0.0);

    jac.set(Ni60, Si27, 0.0);

    jac.set(Ni60, Si28, 0.0);

    jac.set(Ni60, Si29, 0.0);

    jac.set(Ni60, Si30, 0.0);

    jac.set(Ni60, Si31, 0.0);

    jac.set(Ni60, P29, 0.0);

    jac.set(Ni60, P31, 0.0);

    jac.set(Ni60, P33, 0.0);

    jac.set(Ni60, S32, 0.0);

    jac.set(Ni60, S33, 0.0);

    jac.set(Ni60, S34, 0.0);

    jac.set(Ni60, S35, 0.0);

    jac.set(Ni60, Cl33, 0.0);

    jac.set(Ni60, Cl35, 0.0);

    jac.set(Ni60, Cl37, 0.0);

    jac.set(Ni60, Ar36, 0.0);

    jac.set(Ni60, Ar37, 0.0);

    jac.set(Ni60, Ar38, 0.0);

    jac.set(Ni60, K37, 0.0);

    jac.set(Ni60, K39, 0.0);

    jac.set(Ni60, Ca40, 0.0);

    jac.set(Ni60, Ca42, 0.0);

    jac.set(Ni60, Sc41, 0.0);

    jac.set(Ni60, Sc43, 0.0);

    jac.set(Ni60, Sc45, 0.0);

    jac.set(Ni60, Ti44, 0.0);

    jac.set(Ni60, Ti46, 0.0);

    jac.set(Ni60, Ti48, 0.0);

    jac.set(Ni60, V47, 0.0);

    jac.set(Ni60, V49, 0.0);

    jac.set(Ni60, V51, 0.0);

    jac.set(Ni60, Cr48, 0.0);

    jac.set(Ni60, Cr50, 0.0);

    jac.set(Ni60, Cr52, 0.0);

    jac.set(Ni60, Cr54, 0.0);

    jac.set(Ni60, Mn51, 0.0);

    jac.set(Ni60, Mn53, 0.0);

    jac.set(Ni60, Mn55, 0.0);

    jac.set(Ni60, Fe52, 0.0);

    jac.set(Ni60, Fe54, 0.0);

    scratch = screened_rates(k_Fe56_He4_to_Ni60_reaclib)*Y(He4)*state.rho;
    jac.set(Ni60, Fe56, scratch);

    jac.set(Ni60, Fe58, 0.0);

    jac.set(Ni60, Co55, 0.0);

    jac.set(Ni60, Co56, 0.0);

    scratch = screened_rates(k_Co57_He4_to_p_Ni60_reaclib)*Y(He4)*state.rho;
    jac.set(Ni60, Co57, scratch);

    scratch = screened_rates(k_Co59_p_to_Ni60_reaclib)*Y(H1)*state.rho;
    jac.set(Ni60, Co59, scratch);

    jac.set(Ni60, Ni56, 0.0);

    scratch = 0.5*screened_rates(k_Ni58_n_n_to_Ni60_approx)*amrex::Math::powi<2>(Y(N))*amrex::Math::powi<2>(state.rho);
    jac.set(Ni60, Ni58, scratch);

    scratch = -screened_rates(k_Ni60_p_to_He4_Co57_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Ni60_to_He4_Fe56_reaclib)
              - screened_rates(k_Ni60_to_Ni58_n_n_approx)
              - screened_rates(k_Ni60_to_p_Co59_reaclib);
    jac.set(Ni60, Ni60, scratch);

    jac.set(Ni60, Cu59, 0.0);

    jac.set(Ni60, Zn60, 0.0);


    // Cu59 row

    scratch = -screened_rates(k_Cu59_n_to_He4_Co56_reaclib)*Y(Cu59)*state.rho;
    jac.set(Cu59, N, scratch);

    scratch = -screened_rates(k_Cu59_p_to_He4_Ni56_reaclib)*Y(Cu59)*state.rho
              - screened_rates(k_Cu59_p_to_Zn60_reaclib)*Y(Cu59)*state.rho
              + screened_rates(k_Ni58_p_to_Cu59_reaclib)*Y(Ni58)*state.rho;
    jac.set(Cu59, H1, scratch);

    scratch = screened_rates(k_Co55_He4_to_Cu59_reaclib)*Y(Co55)*state.rho
              + screened_rates(k_Co56_He4_to_n_Cu59_reaclib)*Y(Co56)*state.rho
              + screened_rates(k_Ni56_He4_to_p_Cu59_reaclib)*Y(Ni56)*state.rho;
    jac.set(Cu59, He4, scratch);

    jac.set(Cu59, C12, 0.0);

    jac.set(Cu59, O16, 0.0);

    jac.set(Cu59, Ne20, 0.0);

    jac.set(Cu59, Mg23, 0.0);

    jac.set(Cu59, Mg24, 0.0);

    jac.set(Cu59, Al27, 0.0);

    jac.set(Cu59, Si27, 0.0);

    jac.set(Cu59, Si28, 0.0);

    jac.set(Cu59, Si29, 0.0);

    jac.set(Cu59, Si30, 0.0);

    jac.set(Cu59, Si31, 0.0);

    jac.set(Cu59, P29, 0.0);

    jac.set(Cu59, P31, 0.0);

    jac.set(Cu59, P33, 0.0);

    jac.set(Cu59, S32, 0.0);

    jac.set(Cu59, S33, 0.0);

    jac.set(Cu59, S34, 0.0);

    jac.set(Cu59, S35, 0.0);

    jac.set(Cu59, Cl33, 0.0);

    jac.set(Cu59, Cl35, 0.0);

    jac.set(Cu59, Cl37, 0.0);

    jac.set(Cu59, Ar36, 0.0);

    jac.set(Cu59, Ar37, 0.0);

    jac.set(Cu59, Ar38, 0.0);

    jac.set(Cu59, K37, 0.0);

    jac.set(Cu59, K39, 0.0);

    jac.set(Cu59, Ca40, 0.0);

    jac.set(Cu59, Ca42, 0.0);

    jac.set(Cu59, Sc41, 0.0);

    jac.set(Cu59, Sc43, 0.0);

    jac.set(Cu59, Sc45, 0.0);

    jac.set(Cu59, Ti44, 0.0);

    jac.set(Cu59, Ti46, 0.0);

    jac.set(Cu59, Ti48, 0.0);

    jac.set(Cu59, V47, 0.0);

    jac.set(Cu59, V49, 0.0);

    jac.set(Cu59, V51, 0.0);

    jac.set(Cu59, Cr48, 0.0);

    jac.set(Cu59, Cr50, 0.0);

    jac.set(Cu59, Cr52, 0.0);

    jac.set(Cu59, Cr54, 0.0);

    jac.set(Cu59, Mn51, 0.0);

    jac.set(Cu59, Mn53, 0.0);

    jac.set(Cu59, Mn55, 0.0);

    jac.set(Cu59, Fe52, 0.0);

    jac.set(Cu59, Fe54, 0.0);

    jac.set(Cu59, Fe56, 0.0);

    jac.set(Cu59, Fe58, 0.0);

    scratch = screened_rates(k_Co55_He4_to_Cu59_reaclib)*Y(He4)*state.rho;
    jac.set(Cu59, Co55, scratch);

    scratch = screened_rates(k_Co56_He4_to_n_Cu59_reaclib)*Y(He4)*state.rho;
    jac.set(Cu59, Co56, scratch);

    jac.set(Cu59, Co57, 0.0);

    jac.set(Cu59, Co59, 0.0);

    scratch = screened_rates(k_Ni56_He4_to_p_Cu59_reaclib)*Y(He4)*state.rho;
    jac.set(Cu59, Ni56, scratch);

    scratch = screened_rates(k_Ni58_p_to_Cu59_reaclib)*Y(H1)*state.rho;
    jac.set(Cu59, Ni58, scratch);

    jac.set(Cu59, Ni60, 0.0);

    scratch = -screened_rates(k_Cu59_n_to_He4_Co56_reaclib)*Y(N)*state.rho
              - screened_rates(k_Cu59_p_to_He4_Ni56_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Cu59_p_to_Zn60_reaclib)*Y(H1)*state.rho
              - screened_rates(k_Cu59_to_He4_Co55_reaclib)
              - screened_rates(k_Cu59_to_p_Ni58_reaclib);
    jac.set(Cu59, Cu59, scratch);

    scratch = screened_rates(k_Zn60_to_p_Cu59_reaclib);
    jac.set(Cu59, Zn60, scratch);


    // Zn60 row

    jac.set(Zn60, N, 0.0);

    scratch = screened_rates(k_Cu59_p_to_Zn60_reaclib)*Y(Cu59)*state.rho;
    jac.set(Zn60, H1, scratch);

    scratch = screened_rates(k_Ni56_He4_to_Zn60_reaclib)*Y(Ni56)*state.rho;
    jac.set(Zn60, He4, scratch);

    jac.set(Zn60, C12, 0.0);

    jac.set(Zn60, O16, 0.0);

    jac.set(Zn60, Ne20, 0.0);

    jac.set(Zn60, Mg23, 0.0);

    jac.set(Zn60, Mg24, 0.0);

    jac.set(Zn60, Al27, 0.0);

    jac.set(Zn60, Si27, 0.0);

    jac.set(Zn60, Si28, 0.0);

    jac.set(Zn60, Si29, 0.0);

    jac.set(Zn60, Si30, 0.0);

    jac.set(Zn60, Si31, 0.0);

    jac.set(Zn60, P29, 0.0);

    jac.set(Zn60, P31, 0.0);

    jac.set(Zn60, P33, 0.0);

    jac.set(Zn60, S32, 0.0);

    jac.set(Zn60, S33, 0.0);

    jac.set(Zn60, S34, 0.0);

    jac.set(Zn60, S35, 0.0);

    jac.set(Zn60, Cl33, 0.0);

    jac.set(Zn60, Cl35, 0.0);

    jac.set(Zn60, Cl37, 0.0);

    jac.set(Zn60, Ar36, 0.0);

    jac.set(Zn60, Ar37, 0.0);

    jac.set(Zn60, Ar38, 0.0);

    jac.set(Zn60, K37, 0.0);

    jac.set(Zn60, K39, 0.0);

    jac.set(Zn60, Ca40, 0.0);

    jac.set(Zn60, Ca42, 0.0);

    jac.set(Zn60, Sc41, 0.0);

    jac.set(Zn60, Sc43, 0.0);

    jac.set(Zn60, Sc45, 0.0);

    jac.set(Zn60, Ti44, 0.0);

    jac.set(Zn60, Ti46, 0.0);

    jac.set(Zn60, Ti48, 0.0);

    jac.set(Zn60, V47, 0.0);

    jac.set(Zn60, V49, 0.0);

    jac.set(Zn60, V51, 0.0);

    jac.set(Zn60, Cr48, 0.0);

    jac.set(Zn60, Cr50, 0.0);

    jac.set(Zn60, Cr52, 0.0);

    jac.set(Zn60, Cr54, 0.0);

    jac.set(Zn60, Mn51, 0.0);

    jac.set(Zn60, Mn53, 0.0);

    jac.set(Zn60, Mn55, 0.0);

    jac.set(Zn60, Fe52, 0.0);

    jac.set(Zn60, Fe54, 0.0);

    jac.set(Zn60, Fe56, 0.0);

    jac.set(Zn60, Fe58, 0.0);

    jac.set(Zn60, Co55, 0.0);

    jac.set(Zn60, Co56, 0.0);

    jac.set(Zn60, Co57, 0.0);

    jac.set(Zn60, Co59, 0.0);

    scratch = screened_rates(k_Ni56_He4_to_Zn60_reaclib)*Y(He4)*state.rho;
    jac.set(Zn60, Ni56, scratch);

    jac.set(Zn60, Ni58, 0.0);

    jac.set(Zn60, Ni60, 0.0);

    scratch = screened_rates(k_Cu59_p_to_Zn60_reaclib)*Y(H1)*state.rho;
    jac.set(Zn60, Cu59, scratch);

    scratch = -screened_rates(k_Zn60_to_He4_Ni56_reaclib)
              - screened_rates(k_Zn60_to_p_Cu59_reaclib);
    jac.set(Zn60, Zn60, scratch);



    // add ∂λ_weak / ∂Y_e terms now
    // this uses ∂Y_e/∂Y_i = Z_i
    for (int irow = 1; irow <= NumSpec; ++irow) {
        const auto dweak = rate_eval.dweak_ydot_dYe(irow);
        if (dweak == 0.0_rt) {
            continue;
        }
        for (int jcol = 1; jcol <= NumSpec; ++jcol) {
            jac.add(irow, jcol, dweak * NetworkProperties::zion(jcol));
        }
    }

}



template<class MatrixType>
AMREX_GPU_DEVICE AMREX_INLINE
void actual_jac(const burn_t& state, MatrixType& jac)
{

    // Set molar abundances
    amrex::Array1D<amrex::Real, 1, NumSpec> Y;
    for (int i = 1; i <= NumSpec; ++i) {
        Y(i) = state.xn[i-1] * aion_inv[i-1];
    }


    rate_derivs_t rate_eval;

    constexpr int do_T_derivatives = 1;

    evaluate_rates<do_T_derivatives, rate_derivs_t>(state, Y, rate_eval);

    // there are effectively 4 regions in the Jacobian:
    //
    //                    |
    //   1) ∂Ydot_i/∂Y_j  | 3) ∂Ydot/∂e
    //                    |
    //   -----------------+---------------
    //                    |
    //   2) ∂edot/∂Y_j    | 4) ∂edot/∂e
    //                    |


    // region 1

    // Species Jacobian elements with respect to other species
    // J_ij = ∂(Ydot_i)/∂Y_j

    jac_nuc(state, jac, Y, rate_eval);


    // region 2

    // Energy generation rate Jacobian elements with respect to species
    // J_{ienuc,j} = ∂(ε_nuc)/∂Y_j
    //
    // we compute this as ∂(ε_nuc)/∂Y_j = -N_A Σ_i ∂(Ydot_i)/∂Y_j m_i c²

    for (int j = 1; j <= NumSpec; ++j) {
        // create a function that gives us jac(i) for the column j
        auto jac_slice_2 = [&](int i) -> amrex::Real { return jac.get(i, j); };

        amrex::Real tmp{};
        ener_gener_rate(jac_slice_2, tmp);
        jac(net_ienuc,j) = tmp;
    }

    // Account for the weak rate neutrino losses
    // There are two contributions here:
    // ∂(ε_{ν,weak})/∂Y_j and ∂(ε_{ν,weak})/∂Y_e ∂Y_e/∂Y_j
    for (int jcol = 1; jcol <= NumSpec; ++jcol) {
        jac.add(net_ienuc, jcol,
                rate_eval.denuc_weak_dY(jcol) +
                rate_eval.denuc_weak_dYe * NetworkProperties::zion(jcol));
    }

    // Account for the thermal neutrino losses

    amrex::Real dsneutdt{};
#ifdef NEUTRINOS
    amrex::Real sneut, dsneutdd, dsnuda, dsnudz;
    constexpr int do_derivatives{1};
    neutrino_cooling<do_derivatives>(state.T, state.rho, state.abar, state.zbar,
                                     sneut, dsneutdt, dsneutdd, dsnuda, dsnudz);

    for (int j = 1; j <= NumSpec; ++j) {
       amrex::Real b1 = (-state.abar * state.abar * dsnuda + (zion[j-1] - state.zbar) * state.abar * dsnudz);
       jac.add(net_ienuc, j, -b1);
    }
#endif


    // region 3

    // J_{i,ienuc} = ∂(Ydot_i)/∂T for rows i = 1, ..., NumSpec

    amrex::Array1D<amrex::Real, 1, neqs> yderivs;

    // Evaluate the Jacobian elements with respect to energy by
    // calling the RHS using d(rate) / dT and then transform them
    // to our energy integration variable.

    rhs_nuc(state, yderivs, Y, rate_eval.dscreened_rates_dT);

    // convert the derivatives in this column from ∂/∂T to ∂/∂e
    for (int k = 1; k <= NumSpec; k++) {
        jac.set(k, net_ienuc, temperature_to_energy_jacobian(state, yderivs(k)));
    }


    // region 4

    // J_{ienuc, ienuc} = ∂(edot)/∂e.
    // We start by doing the derivatives as ∂/∂T.

    // first the contribution from the strong rates, ∂(ε_nuc)/∂T,
    // which we compute as ∂(ε_nuc)/∂T = -N_A Σ ∂(Ydot_i)/T m_i c²

    amrex::Real jac_e_T{};

    ener_gener_rate(yderivs, jac_e_T);

    // Account for the weak rate neutrino losses
    jac_e_T += rate_eval.denuc_weak_dT;

    // and the thermal neutrino contribution

    jac_e_T -= dsneutdt;

    // finally, convert d(de/dt)/dT to d(de/dt)/de

    jac.set(net_ienuc, net_ienuc, temperature_to_energy_jacobian(state, jac_e_T));

}


AMREX_INLINE
void actual_rhs_init () {


}


#endif
